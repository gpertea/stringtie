# StringTie merge-policy implementation worklog

## Objective and resulting contract

The implementation changes `stringtie --merge` from abundance-weighted,
permissive extension-assembly to structure-first meta-assembly with explicit
support, path, endpoint, guide, and bridge policies.

The resulting defaults are:

- do not filter by FPKM, TPM, or isoform fraction unless the user supplies
  `-F`, `-T`, or `-f`;
- count one structural vote per input file instead of pooling cross-sample
  normalized abundance;
- require one exact shared intron before joining overlapping MET paths;
- keep guide exon structures immutable;
- reject cross-guide-gene and same-gene coding-transcript bridges by default;
- preserve guide IDs, prevent duplicate generated structures and transcript
  IDs, and optionally write a deterministic merge decision report.

The approved implementation plan is recorded in
`worklogs/plan_26-08-09_02-30_merge-dev.md`.

## Option model and command-line handling

`rlink.h` now defines `MergePathMode`, `MergeGuidePolicy`, and `MergeOptions`.
`MergeOptions` is the single typed configuration object for support thresholds,
path mode, SET overlap, biological endpoint rules, guide policy, bridge
exceptions, explicit-filter flags, and the report path. Its constructor contains
the new defaults.

`processOptions()` in `stringtie.cpp` parses and validates every merge-only long
option. Integer support values must be positive where zero has no meaning;
`--min-introns` may be zero; slop and extension values are non-negative, with
`-1` reserved internally for no hard cap; fractions are constrained to `[0,1]`.
The same function rejects these options outside `--merge` and prevents `-o` and
`--report` from overwriting the same path.

The existing global `fpkm_thr`, `tpm_thr`, and `isofrac` values remain shared
with other StringTie modes, but merge mode now resets their unspecified values
to zero. `MergeOptions::fpkm_filter_set`, `tpm_filter_set`, and
`isofrac_filter_set` distinguish an explicit zero from an omitted option.

## Preserving input evidence

`TInputFiles::convert2BAM()` in `tmerge.cpp` still converts each annotation into
a temporary coordinate-sorted BAM stream, which lets the existing bundle engine
perform a multiway merge. The conversion now calls `GffReader::readAll()` with
close-exon merging disabled so input exon boundaries are not normalized before
policy evaluation.

Three BAM auxiliary tags carry merge metadata through the existing alignment
pipeline:

```text
Tag  Payload
---  -------------------------------------------------
ZS   cov|FPKM|TPM, with empty fields preserved
ZG   input gene_id
ZC   CDS start|end
```

`TAlnInfo` in `bundle.h` stores the decoded transcript ID, input-file index,
gene ID, optional abundance fields, CDS bounds, and `qpass`. `qpass` records an
explicit abundance-filter result without immediately deleting the structure;
this deferral is required for replicated rescue.

`processRead()` decodes these fields and marks explicit `-c`, `-F`, and `-T`
failures. Missing FPKM or TPM fails when that metric was explicitly requested,
which avoids treating an absent value as a valid zero or silently skipping the
filter.

## Exact-structure support and abundance filtering

`filter_merge_inputs()` in `rlink.cpp` is the bundle-level evidence gate. It
performs these stages in order:

1. `merge_read_crosses_guide_genes()` checks each same-strand input transcript
   against overlapping guide exons. A cross-gene input is removed before
   abundance or splice-graph processing unless fusions are allowed.
2. For explicit `-f`, a `GHash<double>` accumulates TPM, FPKM, or cov totals by
   input-file and gene ID. If `gene_id` is absent, the current bundle is the
   fallback scope. This prevents a sample's normalized values from being
   compared with another sample's values.
3. `merge_read_key()` serializes strand and the complete ordered exon list.
   `GHash<CMergeReadSupport*>` groups exact structures. Each
   `CMergeReadSupport` owns `GVec<int>` lists of member reads, all contributing
   files, and files containing a passing copy.
4. The group survives when it reaches `--min-support` and has a passing copy,
   or when every copy failed abundance but distinct-file support reaches
   `--rescue-support`. Each rejected or rescued member produces a report event.

Distinct file lists deliberately use gclib `GVec<int>` and the small
`merge_has_int()` helper. This avoids STL containers and guarantees one vote per
file even when a file repeats an exact transcript record.

## Splice graph, path joining, and endpoint evidence

The existing `build_merge()` driver still creates strand-specific splice graphs
and maps input structures to `CTransfrag` and `CMTransfrag` objects. Structural
abundance is now a per-file vote. Exact duplicates can reinforce a path, but
TPM/FPKM magnitudes do not increase graph weight.

`merge_transfrags()` remains the main path constructor. The new mechanics are:

- `merge_shared_introns()` compares graph edges and counts exact intron matches
  between a candidate transfrag and the current path contributors.
- `merge_can_join()` requires `--min-introns` evidence before it considers an
  extension. This also preserves the older intron-compatibility checks already
  enforced by `merge_topath()`.
- `merge_transfrag_support()` counts distinct input files for the candidate.
- `merge_extension_allowed()` first applies the biological-end hard cap, then
  accepts adequate support, and otherwise requires both base-pair and terminal
  exon fraction slop.
- `merge_can_join()` maps genomic left/right to biological 5-prime/3-prime by
  strand, so the independent endpoint options behave correctly on `-` models.

`store_merge_prediction()` creates the resulting `CPrediction`.
`set_merge_prediction_support()` derives these values from contributing
`CReadAln` objects:

- `merge_files`: sorted distinct files contributing any path evidence;
- `merge_path_support`: files carrying the complete ordered intron chain;
- `merge_piece_support`: the minimum distinct-file count over output introns;
- `merge_5p_support` and `merge_3p_support`: files matching each biological end.

The implementation keeps file indices rather than aggregate counts inside
`CPrediction` so later policies and output fractions use consistent evidence.

### Terminal SET/MET policy bypass found by testing

The focused regression fixture found that the legacy local-prediction cleanup
inside `merge_transfrags()` could still widen a MET through an overlapping SET
after `merge_can_join()` correctly rejected the join. There were two causes:

1. cleanup absorption did not call the new intron/endpoint policy;
2. a splice-graph node widened by the SET could give a one-transfrag MET the
   node's outer coordinate rather than its observed coordinate.

`merge_prediction_shared_introns()` and
`merge_prediction_extension_allowed()` now apply the same evidence rules before
cleanup absorbs and extends a prediction. `merge_set_observed_ends()` then reads
`CTransfrag::longstart` and `longend` from all path contributors, clamps the
first and last `GSeg` to the observed union, and recomputes transcript length.
With the default `--min-introns 1`, a MET and terminal SET remain separate;
`--min-introns 0` explicitly permits the legacy union.

## Guide identity, gene boundaries, and bridge barriers

Guide annotations are loaded without close-exon normalization in merge mode.
`registerMergeGuide()` creates a structural key from reference sequence, strand,
and ordered `GSeg` coordinates. Three gclib hashes in `stringtie.cpp` provide
global identity handling:

- `mergeGuideStructures` maps a structure to the lexicographically smallest
  guide ID;
- `mergeGuideAliases` stores other IDs for the same guide structure;
- `mergeGeneratedStructures` records every emitted generated structure.

The canonical guide is used only by `--dedup-guides`; without that option all
reference IDs remain present. Generated structures are always globally unique.

`merge_add_gene_span()` transitively coalesces overlapping guide spans when a
usable gene-level ID is unavailable or merely repeats the transcript ID. This
fallback treats overlapping reference transcripts as one inferred locus but
does not collapse non-overlapping genes into one boundary group.

`merge_generated_barriers()` evaluates a generated `CPrediction` against all
overlapping bundle guides:

1. Same-strand explicit gene IDs plus inferred gene spans determine how many
   guide-gene groups the path crosses. More than one is rejected unless
   `--allow-fusions` is set.
2. Opposite-strand overlaps are recorded as notices and do not block output.
3. For each coding guide, introns present in the prediction are collected.
   `merge_prediction_bridges_guides()` requires a guide-specific intron from
   each of two coding transcripts in the same gene before declaring a coding
   bridge. It is rejected unless `--allow-cds-bridges` is set.

The cross-gene check occurs both on input transcripts and final predictions, so
`--min-introns 0` cannot override guide-gene boundaries.

`merge_related_guide()` links only exact intron-chain containment or extension
variants. `--guide-policy both` retains guide and variant, `reference` retains
only the guide, and `best` replaces the guide only after complete-path support
and `--guide-replace-frac` pass. A changed model never reuses a guide ID.

## Final selection, IDs, and output evidence

`printMergeResults()` is the final policy stage. It sorts predictions, fills
support for exact structures that bypassed path synthesis, then applies guide
identity, bridge, minimum support, path mode, endpoint hard-cap, and guide-policy
checks.

Pairwise cleanup keeps immutable duplicate guides by default, canonicalizes
them when requested, and chooses one generated exact structure. SET/SET
deduplication requires reciprocal `--set-overlap` and selects one existing SET;
it never creates a transitive range union. A generated SET fully contained in a
retained MET exon is removed.

Before output, `registerMergeGeneratedStructure()` enforces process-wide
structural uniqueness. A local `GHash<bool> outputIds` independently prevents a
transcript ID from being printed twice. Output transcript records expose the
support counts and fractions, related guide, class code, and guide aliases.
Unbundled guides receive the same merge evidence schema with zero support.

## Deterministic rejection report

Workers call `mergeReportEvent()` to append `CMergeReportEvent` objects to the
owning `GPVec`; `mergeReportMutex` protects the append in threaded builds. Each
event stores reference, strand, range, action, reason, transcript ID, input file,
guide IDs, gene IDs, support, and details.

`writeMergeReport()` runs after bundle processing. `CMergeReportEvent::operator<`
defines a total key over every output field, so thread scheduling cannot change
row order. The function sorts the vector, counts reasons with `GHash<int>`,
sorts the reason names, writes summary comment lines, and then writes TSV rows.

## Supporting fix

The temporary BAM records built by `GSamRecord` gained more auxiliary data.
Under debug/Valgrind this exposed pre-existing uninitialized `_cigar` and `_read`
pointers in the build-from-scratch constructor. `gclib/GSam.cpp` now initializes
both debug members to null before any destructor path can inspect them.

## Tests

`test_merge_policy.sh` creates small self-contained GTF fixtures in a temporary
directory and checks:

- explicit low-FPKM filtering plus replicated structural rescue;
- minimum support, direct, replicated, and shared path behavior;
- default SET/MET terminal-overlap rejection and `--min-introns 0` relaxation;
- independent plus-strand 3-prime support and hard caps;
- independent guide 5-prime/3-prime caps on plus and minus strands;
- flexible hard-cap relaxation and `reference`, `both`, and `best` guide modes;
- cross-gene rejection even when the intron gate is disabled;
- explicit fusion and same-gene CDS-bridge authorization;
- identical-guide canonicalization and alias output;
- absence of duplicate transcript IDs.

The `Makefile` target `test-merge` builds StringTie and runs this suite. The
script also supports `STRINGTIE_BIN`, `VALGRIND=1`, and `KEEP_TMP=1`.

`compare_merge_outputs.pl` is a reusable structure-level comparison tool. It
parses records without relying on transcript-ID uniqueness, counts SETs, METs,
intron chains, exact structures, duplicate groups, exact/changed guide IDs, and
exact-source versus synthesized structures, then compares two output structure
sets.

## Full-data comparison

The 20 issue #247 GTFs were run with the same guide under a clean 3.0.3 build
from `dae0eff` and the new defaults. The measured report is in
`worklogs/merge-default-vs-3.0.3.md`.

The new default produced 192,173 transcript records versus 63,959 in 3.0.3, an
increase of 128,214. Exact source structures increased from 20,290 to 172,682,
while synthesized structures decreased from 43,577 to 19,461. All 26,059 guide
IDs retained their exact guide structure; 3.0.3 changed 13,958. Neither output
had duplicate transcript IDs, and neither had a duplicate group containing a
generated structure.

## Successful validation and analysis commands

```sh
make -j4 release
make test-merge
VALGRIND=1 ./test_merge_policy.sh
sh -n test_merge_policy.sh
perl -c compare_merge_outputs.pl
git diff --check
./compare_merge_outputs.pl "$GUIDE" v3.0.3.gtf new-default.gtf $INPUTS
```

A clean 3.0.3 release tree was built from `git archive HEAD`; both full-data
runs completed successfully. The code graph was then rebuilt in full mode as
`home-gpertea-src-stringtie-git` and reported 15,255 nodes, 86,887 edges, and
`status: indexed`.
