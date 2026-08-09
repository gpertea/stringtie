# StringTie merge mode

`stringtie --merge` combines transcript annotations from multiple GTF/GFF
inputs into a structurally non-redundant isoform set. It operates on assembled
transcripts, not raw RNA-seq reads.

## Workflow

```text
Input GTF/GFF files
  -> preserve tx ID, input file, gene ID, CDS, and available abundance fields
  -> convert transcripts to temporary alignment records
  -> coordinate-merge streams and form genomic bundles
  -> reject input bridges across distinct -G genes
  -> count exact structures by distinct input file
  -> apply optional abundance filters, with replicated rescue
  -> build strand-specific splice graphs and transfrags
  -> combine intron-compatible paths under the selected evidence policy
  -> enforce guide, CDS-bridge, endpoint, and path-support barriers
  -> deduplicate generated structures and write merged GTF/report
```

### 1. Input preparation and filtering

- `TInputFiles::start()` accepts individual files or a file containing their
  names. `TInputFiles::next()` performs the coordinate-ordered multiway merge.
- `TInputFiles::convert2BAM()` preserves the input exon structure and carries
  transcript ID, input-file index, gene ID, CDS bounds, cov, FPKM, and TPM in
  the temporary record.
- `processRead()` records whether an input passes explicit `-c`, `-F`, or `-T`
  filters. It defers removal until structural support is known.
- `filter_merge_inputs()` counts each exact structure once per input file.
  `--min-support` removes insufficiently replicated structures and
  `--rescue-support` can restore a structure rejected only by abundance filters.
- An explicit `-F` or `-T` requires that field to be present and pass. An
  explicit `-f` is evaluated within input file and input gene, or within the
  current locus when `gene_id` is absent. It prefers TPM, then FPKM, then cov.

The merge defaults are `-F 0`, `-T 0`, and `-f 0`. Cross-sample TPM/FPKM values
are not used as splice-graph weights; each input file contributes one structural
vote.

### 2. Genomic bundles and core merge

`main()` forms bundles and `processBundle()` calls `infer_transcripts()`, which
selects `build_merge()` in merge mode.

`build_merge()` is the main driver. It adds guides with `guide_to_read()`,
collapses exact duplicate input structures, groups overlapping transcripts,
creates strand-specific bundles, and sends graph-backed bundles through these
deeper workers:

- `create_graph()`: builds splice-graph nodes and edges.
- `construct_treepat()`: builds the graph-pattern lookup.
- `get_read_to_transfrag()`: maps each input tx to a graph path.
- `get_read_pattern()` and `update_abundance()`: encode paths and structural
  votes.
- `merge_transfrags()`: joins compatible paths subject to intron and endpoint
  evidence.
- `add_transfrag_to_path()`: extends a candidate graph path.
- `store_merge_prediction()`: creates a `CPrediction`.
- `merge_set_observed_ends()`: restores input terminal coordinates when graph
  nodes were widened by an overlapping SET.
- `set_merge_prediction_support()`: records distinct-file, whole-path,
  intron-piece, and biological 5-prime/3-prime support.

`process_merge_transfrags()` is present but not called by the active path.

### 3. Structural policies

The default shared-path policy requires at least one exact shared intron before
two MET paths can be combined. All overlapping introns must remain compatible.
SETs therefore do not extend METs merely through terminal-exon overlap. SETs
are deduplicated by reciprocal overlap without constructing a transitive union;
a non-guide SET fully contained in a MET exon is removed.

The intron and endpoint checks apply both while `merge_transfrags()` builds a
path and while its legacy local-prediction cleanup considers absorbing one model
into another. This closes a former path around the policy. Prediction ends are
then clamped to observed `CTransfrag::longstart` and `longend` coordinates so a
graph node widened by a neighboring SET cannot silently widen an unrelated MET.

`--min-introns 0` restores terminal-overlap path joining, but it never bypasses
guide-gene, CDS-bridge, hard endpoint, or support barriers. Guide-gene barriers
have precedence over shared-intron evidence.

Biological 5-prime and 3-prime ends are handled independently. A new extension
must meet its end-specific support threshold unless both its base-pair and
terminal-exon-fraction overhang are within the configured slop. Hard extension
caps apply regardless of support.

The path modes are:

```text
Mode        Requirement
----------  --------------------------------------------------------------
shared      Compatible paths share at least --min-introns exact introns
replicated  Every output intron has --min-piece-support input files
direct      The complete ordered intron chain has --min-support input files
```

### 4. Guide and bridge policies

`-G` guide exon structures are loaded without close-exon normalization and are
never rewritten. A changed model never reuses a guide transcript ID.

Barrier order is:

1. Reject a same-strand path crossing distinct guide genes unless
   `--allow-fusions` is set.
2. Reject a same-gene path combining unique introns from distinct coding guides
   unless `--allow-cds-bridges` is set.
3. Report opposite-strand guide overlaps without treating them as a barrier.
4. Apply biological 5-prime/3-prime hard extension caps.
5. Apply shared-intron, path, piece, and endpoint support policies.

When a reference annotation lacks a usable gene-level ID and repeats the
transcript ID as `gene_id`, overlapping same-strand guide spans are treated as
one inferred guide locus. Non-overlapping loci remain separate fusion barriers.

Guide output policies are:

```text
Policy     Result for a same-chain endpoint variant
---------  ---------------------------------------------------------------
both       Keep the unchanged guide and emit the variant under a new ID
reference  Keep the unchanged guide and suppress the endpoint variant
best       Replace the guide only with replicated, high-core-support evidence
```

Exact guide matches carry `reference_id` and `class_code "="`. Pure
same-chain containment/extensions carry `reference_id` and `class_code "k"`.
Identical guide structures keep every reference ID by default. With
`--dedup-guides`, the lexicographically smallest ID is kept and the others are
reported in `reference_aliases`. Generated structures are always deduplicated.

### 5. Output and report

`printMergeResults()` applies final barriers, structural deduplication, guide
identity rules, and SET containment before writing GTF. Merge transcript
records include:

```text
merge_support          distinct files contributing to the prediction
merge_path_support     files containing the complete intron chain
merge_piece_support    minimum distinct-file support over output introns
merge_5p_support       files supporting the biological 5-prime endpoint
merge_3p_support       files supporting the biological 3-prime endpoint
merge_locus_support    distinct files represented in the locus
merge_5p_fraction      5-prime support / merge_support
merge_3p_fraction      3-prime support / merge_support
reference_id           related guide, when present
class_code             "=" for exact guide, "k" for pure containment
reference_aliases      suppressed identical guide IDs, when requested
```

`--report FILE` writes a deterministic TSV after a comment summary. It records
filtering, replicated rescue, fusion and CDS-bridge decisions, opposite-strand
overlaps, containment, duplicate suppression, support, guide/gene IDs, input
transcript ID, and input file number.

## Merge-only controls

```text
Option                    Default
------------------------  ----------
--min-support             1
--rescue-support          2
--min-introns             1
--path-mode               shared
--min-piece-support       2
--set-overlap             0.80
--min-support-5p          2
--min-support-3p          2
--slop-bp-5p              100
--slop-bp-3p              100
--slop-frac-5p            0.50
--slop-frac-3p            0.50
--max-extension-5p        unlimited
--max-extension-3p        unlimited
--guide-policy            both
--guide-replace-frac      0.90
--dedup-guides            off
--allow-fusions           off
--allow-cds-bridges       off
--report                  unset
```

All controls above are rejected outside `--merge` mode.

### Support and assembly options

- `--min-support N` requires an exact input structure to occur in at least `N`
  distinct input files. It is also the complete-chain threshold for
  `--path-mode direct`.
- `--rescue-support N` restores an exact structure when every observed copy
  failed an explicitly requested abundance filter but the structure occurs in
  at least `N` files. It does not rescue a structural or guide barrier.
- `--min-introns N` requires `N` exact shared introns before compatible MET
  paths can be joined. The default `1` blocks SET-only terminal overlap;
  `0` deliberately restores terminal-overlap joining.
- `--path-mode shared` uses the shared-intron gate. `replicated` also requires
  every output intron to be replicated. `direct` requires replication of the
  complete ordered intron chain.
- `--min-piece-support N` is the per-intron distinct-file threshold used only
  by `--path-mode replicated`.
- `--set-overlap F` sets the reciprocal overlap fraction for choosing one of
  two generated SETs. It does not create a union of their ranges.

### Endpoint options

- `--min-support-5p N` and `--min-support-3p N` set independent distinct-file
  thresholds for an endpoint introduced by path joining.
- `--slop-bp-5p N` and `--slop-bp-3p N` set the maximum unsupported overhang in
  bases. An under-supported join must also pass the corresponding fractional
  slop.
- `--slop-frac-5p F` and `--slop-frac-3p F` set the maximum unsupported
  overhang as a fraction of the candidate terminal exon. Both the base-pair and
  fraction conditions must pass.
- `--max-extension-5p N` and `--max-extension-3p N` are absolute hard caps.
  They apply regardless of support; omission means unlimited.

Endpoint names are biological, not genomic: genomic right is 3-prime on `+`
and 5-prime on `-`. Support and slop govern synthesized joins. A directly
observed same-chain guide variant is controlled by the hard caps and guide
policy because its full endpoint is itself observed.

### Guide and barrier options

- `--guide-policy both` preserves an unchanged guide and emits an accepted
  same-chain endpoint variant under a generated ID.
- `--guide-policy reference` preserves the guide and suppresses its same-chain
  endpoint variant.
- `--guide-policy best` replaces the guide only when at least two files contain
  the variant's complete chain and complete-path support divided by total
  variant support reaches `--guide-replace-frac`.
- `--guide-replace-frac F` sets that `best` core-support fraction.
- `--dedup-guides` emits one canonical record for structurally identical guides.
  The lexicographically smallest transcript ID wins and aliases are retained in
  `reference_aliases`.
- `--allow-fusions` changes a distinct-guide-gene crossing from rejection to an
  accepted, reported fusion. This barrier otherwise precedes intron evidence.
- `--allow-cds-bridges` permits an intra-gene path that combines guide-specific
  introns from two coding guide transcripts.
- `--report FILE` writes the sorted decision events and per-reason counts to a
  TSV file after all worker threads have completed.

## Active call path

```text
main
  -> TInputFiles::start / convert2BAM / next
  -> BundleData::evalReadAln
  -> processRead
  -> processBundle
       -> infer_transcripts
            -> build_merge
                 -> filter_merge_inputs / guide_to_read
                 -> add_read_to_group / merge_read_to_group
                 -> create_bundle / add_group_to_bundle
                 -> create_graph / construct_treepat
                 -> get_read_to_transfrag
                      -> get_read_pattern / update_abundance
                 -> merge_transfrags
                      -> add_transfrag_to_path
                      -> store_merge_prediction
                           -> merge_set_observed_ends
                           -> set_merge_prediction_support
       -> printMergeResults
```

## Regression tests

```sh
make test-merge
VALGRIND=1 ./test_merge_policy.sh
```

`test_merge_policy.sh` builds temporary locus-level GTFs and checks support
rescue, path modes, SET/MET terminal overlap, independent biological endpoints,
flexible guide caps, guide policies, gene and CDS bridge barriers, guide dedup,
and transcript-ID uniqueness. `compare_merge_outputs.pl` compares two merge
GTFs at exact-structure, intron-chain, guide-identity, and source-versus-
synthesized levels.
