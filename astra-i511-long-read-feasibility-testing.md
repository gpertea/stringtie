# Issue 511 long-read bottleneck and test-data feasibility

Date: 2026-11-12

Repository: `/home/gpertea/src/stringtie`

Inspected revision: `53d964779c5782e0687ae9907dc0f31f6748d355` (`master`,
merge of PR #512 / merge-dev).

Status: investigation and feasibility assessment. No optimization, synthetic
dataset, new build, or StringTie benchmark was performed during this investigation.
The BAM examinations below were read-only. This document records both completed
inspection and proposed follow-up work; the latter is not an approved implementation
plan.

## Main findings

Ida's proposed hash-based acceleration addresses a real quadratic search in
`processRead()`. The reporter's experimental candidate index provides additional
evidence that this search is important on their data. A behavior-preserving fix is
nevertheless more subtle than hashing identical raw alignments: guide-aware junction
handling, strand assignment, and terminal-exon trimming affect the existing matching
semantics.

Local data is sufficient for a useful investigation without obtaining new reads or
rerunning an aligner. In particular, `issue_476/bundle_inf.bam` already contains a
natural same-start hotspot: 75,204 primary alignments at one coordinate, representing
7,147 distinct exon/strand patterns before StringTie processing. Start with this
natural input and controlled subsets. Synthetic inputs can then extend structural
diversity and isolate correctness cases.

Success requires both reduced producer work and unchanged behavior. Matching bundle
counts alone is insufficient: the chosen collapse representative, accumulated
evidence, transcript structures, and abundance estimates must also agree.

## Project context and source inspection

Read `AGENTS.md`, `README.md`, and `worklogs/handoff-stringtie-codebase.md`.
The handoff predates the merge-dev merge; its branch status and some code/test
descriptions are obsolete. Current source, rather than the handoff's line numbers,
was used for this investigation.

The important control flow is:

```text
main: read coordinate-sorted alignments
  -> BundleData::evalReadAln: guide overlaps, strand evidence, counts
  -> processRead: filtering, trimming, alignment collapse, junctions, pairing
  -> finish and queue the bundle
worker: processBundle
  -> "begins processing" log message
  -> infer_transcripts
  -> print results and clear bundle
```

The producer is a single thread. Idle assembly workers and a blocked upstream
`samtools` are consistent with producer-side computation, but do not identify a
specific function by themselves.

The code graph project `home-gpertea-src-stringtie` was discovered with
`list_projects`. `index_status` reported a ready, watched index rooted at this
repository, indexed at `2026-11-12T15:01:43Z`. Coverage checks reported matching
metadata and no recorded issue for `rlink.cpp`, `bundle.cpp`, `bundle.h`,
`stringtie.cpp`, `tablemaker.h`, `gclib/GList.hh`, and `gclib/GSam.cpp`.
Graph symbol searches and a direct caller/callee trace of `processRead` were
followed by exact source reads. Coverage is best-effort, not proof of completeness;
the graph also contains old source copies in debugging directories.

Relevant source locations at the inspected revision:

```text
Location                 Relevant behavior
-----------------------  --------------------------------------------------------
rlink.cpp:239            deljuncmatch: stored effective vs incoming deletion offsets
rlink.cpp:252            exonmatch: exact exon-coordinate comparison
rlink.cpp:663            processRead entry and long-read filtering
rlink.cpp:707            terminal polyA/T exon-trimming block
rlink.cpp:761            backward same-start candidate search
rlink.cpp:806            new stored alignment creation and polyA/T evidence
rlink.cpp:857            sorted junction insertion through AddIfNew
rlink.cpp:866            append new alignment to readlist
bundle.cpp:172          evalReadAln and guide-overlap processing
bundle.h:382            BundleData, including readlist and junction containers
bundle.h:464            bundle cleanup and readlist reset
stringtie.cpp:1153      evalReadAln followed by processRead
stringtie.cpp:1859      processBundle and processing-start log
tablemaker.h:329        findOvlExons: guide-exon overlap scan and cache
gclib/GSam.cpp:201      setupCoordinates: CIGAR to exons and deletion offsets
gclib/GSam.cpp:355      sequence: allocate and decode the entire BAM sequence
gclib/GList.hh:517      AddIfNew
gclib/GList.hh:619      sortInsert and pointer-array movement
```

The working directory already contained modified `.gitignore` and `AGENTS.md`,
many untracked data/build files, and multiple existing binaries. They were not
changed. `worklogs` is a symlink to `../stringtie-git/worklogs`; this document is
intentionally written through the requested path.

## Evidence from GitHub issue 511

Source: [issue #511](https://github.com/gpertea/stringtie/issues/511).

The reporter compared Conda StringTie 2.2.2 and 2.2.3 on the same complete
coordinate-sorted chr2 BAM, containing 927,114 alignment records, and the same GTF:

```bash
samtools view -h chr2.bam | \
  /path/to/stringtie -L -c 2 -v -p 4 \
    -G chr2.gtf -l test_chr2 -o output.gtf -
```

```text
Measurement                          2.2.2       2.2.3
-----------------------------------  ----------  ----------
StringTie wall time                  338.53 s    529.34 s
StringTie CPU time                   376.77 s    568.83 s
Gap before hotspot processing log    322 s       513 s
Hotspot processing start to done     about 2 s   about 1 s
```

The hotspot was `chr2:88,538,720-89,330,679`, with 634,090 alignments,
337,129 internal distinct structures, and 12,859 junctions. During the gap the main
thread used about one CPU, four assembly workers were mostly idle, and upstream
`samtools` was observed blocked in `pipe_write`. All 2,198 bundles matched in
coordinates and reported counts across the two versions.

A separate chr22 hotspot at `21,867,909-23,145,021` showed pre-start gaps of
91/151 seconds and processing times of roughly 2/3 seconds. Historical full-depth
times were much larger, but resource conditions were uncontrolled. The reporter
explicitly did not claim a proven version regression or function-level profile.

[Ida's reply](https://github.com/gpertea/stringtie/issues/511#issuecomment-6825744958)
identifies the same-start linear search, explains its quadratic worst case, and
proposes replacing it with a hash lookup. She attributes the loci to highly abundant
immunoglobulin transcripts with varied structures. That biological explanation was
not independently established from the reporter's reads in this investigation.

[The reporter's follow-up](https://github.com/gpertea/stringtie/issues/511#issuecomment-6826112393)
describes an experimental v2.2.2 candidate index retaining the original matching
checks:

```text
Dataset   Original   Experimental index
--------  ---------  ------------------
chr2      459.70 s   42.16 s
chr22      95.53 s   33.90 s
```

These are preliminary comparisons from different build environments/system loads.
Per-bundle input and potential-transcript counts matched; full output equivalence
was still pending. Neither the reporter's patch nor Ida's patch was inspected.

The bundle-wide count of 337,129 structures does not reveal how many share each
start. The distribution by start, and actual candidate-comparison counts, are
needed to quantify the quadratic explanation.

## Matching semantics that constrain a safe optimization

### The search window and ordering matter

In non-merge mode, `processRead` starts at the last stored alignment and walks
backward only while stored starts equal the incoming `brec.start`. It chooses the
first match encountered, which is the newest eligible representative. In normal
producer operation, `readlist` is appended to rather than maintained in sorted
order; workers can sort it later.

With K new distinct structures at one start, failed candidate checks can total
K*(K-1)/2. With N incoming reads and K retained structures, duplicate lookup can
also approach O(N*K), depending on representative order and which structures
recur. Exon/deletion comparisons add work proportional to the compared structure.

A hash index must preserve the eligible trailing run and newest-first choice if
the objective is a pure performance change. A bundle-wide lookup can find older
matches that the existing scan never considers.

### Equality is not raw CIGAR equality

Candidate matching checks:

- incoming vs stored strand;
- short/long-read status;
- unitig status when unitig handling is enabled;
- exact exon coordinates;
- junction-deletion compatibility in long-read or mixed mode.

Different CIGAR strings may map to the same exons. Conversely, identical exons may
need separate stored records because their junction-deletion evidence differs.
Changes must preserve updates to `nh`, read counts, polyA/T evidence, junction
counts, and pairing relationships, not just the number of stored records.

### Junction matching is asymmetric

`deljuncmatch` reconstructs deletion offsets from the stored exon and junction
coordinates, then compares those offsets to the incoming `brec.juncsdel`.
When a junction matches a guide, creation of the stored junction can suppress
the deletion-based coordinate adjustment. The stored effective offsets can
therefore differ from the original incoming offsets.

Hashing both sides' original deletion tags as a symmetric identity can change
matches. Keeping the original predicate as a final check is a conservative
approach, but the candidate key must also avoid excluding valid matches.

### Trimming and inferred strand affect indexing

Current master may remove a terminal polyA/T artifact exon before searching and
then update `brec.start` and `brec.end`. Coordinate-sorted input does not imply
monotonic post-trimming starts. Clearing an index merely because an incoming start
changes, or searching all earlier records with that start, needs careful comparison
against the existing trailing-run behavior.

For a newly stored long read, polyA/T evidence can assign a strand after lookup.
The index must reflect the stored matching state, which can differ from the
incoming strand used for the unsuccessful search. Unitig state is also populated
later in `processRead`.

The same-start loop is present in v2.2.3 and current master, as confirmed by
reading both. Their surrounding long-read processing is not identical, so a
2.x optimization cannot be assumed correct on master without these cases.

### Lifetime and other potential producer costs

Any index containing readlist indices or pointers needs a clear producer-only
lifetime and reset at bundle boundaries/reuse. It must not survive worker
reordering or freeing the reads. Merge mode bypasses this collapse loop and
should retain its current behavior.

Other source-level candidates for residual producer cost are:

- guide-exon overlap scans in `evalReadAln` / `findOvlExons`;
- sorted junction insertion, whose array insertion can move pointers;
- repeated calls to polyA/T helpers that allocate and decode the full read
  sequence through `GSamRecord::sequence()`.

These are plausible costs, not measured bottlenecks. Profile after accelerating
the candidate search before expanding the optimization scope.

## Local BAM inventory and measurement method

Discovery used `rg --files --hidden` under this repository, excluding `.git` and
irrelevant vendored directories. It found 36 BAM paths, one CRAM
(`eN_test/R6532.cram`), and one SAM (`mix/alns1.sam`). The numerical inventory below
covers BAMs only. Paths were not deduplicated by content or inode.

`samtools` was available at `/opt/sw/bin/samtools`. The default `python3` had no
`pysam`, so inspection used Python's standard library to consume `samtools view`.
No package installation was needed.

The initial pass retained at most 5,000 primary mapped alignments per file,
excluding flags `0x4`, `0x100`, and `0x800` (mask 2308). Query length was computed
from CIGAR operations consuming query sequence (`M`, `I`, `S`, `=`, `X`), including
soft clips and excluding hard clips. Median and P90 are sorted values at indices
`n//2` and `int(0.9*n)`. These are prefix samples in coordinate order, not random
samples or whole-file length distributions.

Every sampled primary record had sequence present. Length alone does not establish
sequencing technology; filename hints are also insufficient. For example,
`bundle2_long_c10A.bam` has a 150 bp median in its first 5,000 primary records.

Initial inventory (MiB means 1,048,576 bytes; rows with n < 5,000 reached EOF):

```text
Path                                             MiB        n  Median   P90    Max
--------------------------------------------  -------  ------  ------  ----  -----
usg_test/chrX.bam                               249.17    5000     101   101    101
e_test/missed_regions.bam                         0.02      69     100   100    100
bundles/bundle_sph_B.bam                        146.32    5000     101   101    101
bundles/bundle_sph_A.bam                        113.47    5000     101   101    101
ntest/chr1.bam                                   5.62    5000     100   100    100
issue_357/eLsel.bam                              0.83     495    1350  3024   6531
HPGL/chr1XY.bam                                252.36    5000     101   101    101
mix/c21_2bundles.bam                             1.00    5000      76    76     76
mix/chr21_short.bam                            667.35    5000      76    76     76
mix/chr1_lr.bam                                972.13    5000     282   405   8045
mix/chr21_long.bam                             105.92    5000     823  2130  11635
mix/chr1_sr.bam                               1286.17    5000      76    76     76
mix/chr21_2bundles.bam                            2.12    5000      76    76     76
issue_434/small_BAM_SR.bam                        0.02      44      37    37     37
issue_434/small_BAM_PE.bam                        0.02      86      37    37     37
eN_test/mix_short.bam                            0.64    5000      76    76     76
eN_test/missed_regions.bam                        0.02      69     100   100    100
issue_238/sel_short.bam                          1.08    5000     100   100    100
mixbug/chr1.long.bam                            32.12    5000     725  1568   6368
mixbug/chr1.short.bam                           81.59    5000      76    76     76
issue_465/bundle2_short_c10A.bam               7354.24    5000     101   101    101
issue_465/bundle2_long_c10A.bam                 240.16    5000     150   150   1301
debug_L320/bundle_c1.bam                      1561.74    5000     566   707   1937
issue_476/bundle_inf.bam                        91.92    5000    2308  3932   6532
issue_476/4801.heart.female.nosup.sorted.bam    2614.57    5000    1535  2334   5745
issue_465/shared_bam/iso.bundle.c19.238-244.bam    0.10     166    2386  4400   5719
issue_465/shared_bam/bundle.c19.238-244.bam       11.51    5000     110   151    151
readgroup_bug/iso.bundle.c19.238-244.bam           0.10     166    2386  4400   5719
readgroup_bug/bundle.c19.238-244.bam              11.51    5000     110   151    151
bugL357/eLsel.bam                                0.83     495    1350  3024   6531
bugL357/bundle2.bam                           1745.75    5000    1305  1826   5955
bugL357/long.bam                                 3.91    4005     960   997   3362
bugL357/r1.bam                                   0.08      51    1427  2209   3362
bugL357/b1_4893000_4898000.bam                    3.91    4005     960   997   3362
bugL357/bundle1.bam                            786.99    5000    1277  1960   6213
bugL357/bundle.bam                             782.20    5000     960  1429   5502
```

Several similarly named pairs have identical sample statistics and sizes. They may
be copies, but content identity was not established. Do not count them as independent
validation datasets without checking.

## Deeper examination of candidate datasets

For six selected inputs, exon blocks were reconstructed from CIGAR reference
consumption, split at `N`, and grouped by reference and alignment start. A shape was
the exon-coordinate tuple plus strand. Strand used `XS` when present, otherwise
`ts`, with reverse-alignment orientation applied to `ts`; absent tags gave `.`.

This exploratory parser is not a replacement for `GSamRecord::setupCoordinates`.
It does not reproduce special handling of unusual intron CIGARs, deletion-aware
matching, guide-derived strand, or polyA/T trimming. In particular, raw shape
counts must not be described as StringTie collapse counts. No `N-I-N` CIGAR pattern
was found in the complete issue #476 scan.

```text
Dataset                          Primary examined  Median   P90    Max  Shapes/start
-------------------------------  ----------------  ------  ----  -----  ------------
issue_476/bundle_inf.bam               173655, all    2322  4029   7013          7147
debug_L320/bundle_c1.bam                50000        552   771  75656          3308
mixbug/chr1.long.bam                    23536, all     748  1665   7662           389
bugL357/long.bam                         4005, all     960   997   3362            86
issue_357/eLsel.bam                       495, all    1350  3024   6531            13
readgroup_bug/iso.bundle.c19.238-244.bam   166, all    2386  4400   5719            22
```

### Issue 476: strongest natural performance candidate

Complete `issue_476/bundle_inf.bam` scan:

```text
Measurement                              Count
---------------------------------------  -------
Total records / primary mapped records   173655
Spliced primary alignments               172909
Soft-clipped primary alignments           41107
Junction-adjacent deletion alignments       2383
Alignments without XS/ts strand evidence       0
Raw exon/strand shapes across all starts   31437
```

Dominant starts, sorted by number of shapes:

```text
Reference  Start       Reads   Exon/strand shapes
---------  ----------  ------  ------------------
chr2       17815963     75204                7147
chr2       17843707     33370                3716
chr2       17843708     22827                3244
chr2       17815964     11131                2830
chr2       17843709     13118                2512
```

The first 50,000 primary records alone contained 43,754 reads and 5,558 shapes at
`chr2:17815963`; this was an intermediate sample, superseded by the complete counts
above. The first 5,000 records included exon counts from 1 through 39, providing
substantial multi-exon structure variety.

The BAM header declares coordinate sorting. Header history records extraction from
`4801.heart.female.nosup.sorted.bam` over `chr2:17677819-17894908`.
Do not assume that `chr2` denotes the same species/reference as the human issue #511
input; the local reference lengths differ from human GRCh38.

Existing `issue_476/log_bundle_inf_fix.log` records a historical StringTie 3.0.2
run with `./stringtie_ida_fix -v -L bundle_inf.bam`:

- 173,651 accepted alignments, 33,280 stored distinct entries, 2,570 junctions,
  and no guides;
- processing began at 16:36:50 and finished at 16:36:52;
- 121 processed potential transcripts;
- 173,655 long-read alignments examined, four discarded in singleton screening,
  and three exon-shortened.

The difference between the 31,437 raw shapes and 33,280 historical stored entries
is not an inconsistency: the measurements use different predicates and processing.
These old logs demonstrate a manageable inference case, but do not establish
current producer timing, current output, or which historical fix the binary held.

### debug_L320: substantial diversity with more complex alignment noise

For the first 50,000 primary records:

- 33,587 were spliced; all 50,000 were soft-clipped;
- 23,500 contained a deletion immediately adjacent to an intron;
- 29,339 distinct exon/strand shapes occurred across all starts;
- the leading starts were `1:6642` (8,051 reads, 3,308 shapes), `1:6648`
  (6,542 reads, 2,805 shapes), and `1:6645` (2,657 reads, 1,381 shapes).

The file is about 1.53 GiB. Its index reports 2,832,287 mapped alignments on reference
`1`. Header history identifies minimap2 2.24-r1122 with splice mode, `--secondary=no`,
`--MD`, and a junction BED. Local annotation/log names identify zebrafish reference
material. This is useful for deletion and trimming stress, but a bounded subset is
preferable to the full file initially.

Historical logs show different accepted counts: 2.2.2 recorded 2,832,287 alignments
and 1,816,971 distinct entries; a 3.0.2 run recorded 1,328,115 alignments and 899,345
distinct entries. Both reported zero guides in the bundle despite a supplied GTF.
Do not interpret these logs as a controlled timing comparison or assume the supplied
annotation actually exercises guide matching.

### Smaller correctness and mixed-mode inputs

`mixbug/chr1.long.bam` has 23,536 primary records, 16,519 spliced records,
3,340 with junction-adjacent deletions, and 7,182 without strand tags. All primary
records are soft-clipped. The top start, `chr1:1013497`, has 2,480 reads and 389
shapes. Its header identifies NA12878 cDNA alignment material. The directory also
has `chr1.short.bam`, `chr1.annotation.gtf`, and a recorded `--mix` command.

`bugL357/long.bam` has 4,005 primary records, 3,972 spliced, 739 with
junction-adjacent deletions, and 34 without strand tags. All are soft-clipped.
Its top start, `Chr1:4892536`, has 1,208 reads and 86 exon/strand shapes, despite
1,207 distinct CIGAR strings. This directly illustrates why unique CIGAR counts
overstate relevant structure diversity. Header history identifies Arabidopsis
Col-CEN alignments. The directory contains several guide/assembly GTF/GFF files;
verify their intended role before treating an assembly as reference annotation.

`issue_357/eLsel.bam` has 495 primary records, 372 spliced, 99 with
junction-adjacent deletions, and 123 without strand tags. It has a matching
`eLsel.gtf` and a recorded expression-only long-read command:

```bash
stringtie3 -e -L -o eLsel_strg3_out.gtf -G eLsel.gtf eLsel.bam
```

`readgroup_bug/iso.bundle.c19.238-244.bam` has 166 primary records, 140
spliced, 100 soft-clipped, and none with junction-adjacent deletions or missing
strand tags. It has a companion short-read BAM and apparent copies under
`issue_465/shared_bam`. It is useful as a small mixed-input check, not a convincing
performance stress case by itself.

## Feasible benchmark construction

### 1. Keep an unchanged natural baseline

Use the complete issue #476 bundle first, with its original sequences, CIGARs,
tags, and record order. It already provides thousands of shapes at common starts
and a historical example of short inference time. Measure the current binary built
from the inspected source before making assumptions from older local executables.

Keep source BAMs read-only. Place derived inputs, logs, and outputs in a dedicated
experiment directory. Record source revision, build flags, input identities,
commands, and the exact record-selection rule.

### 2. Derive controlled subsets from real records

Select records at the dominant start, varying the number of distinct structures
and the number of repeats independently. Preserve original relative order among
equal-coordinate records for the principal subset series. A separate ordering
experiment can intentionally group or interleave duplicate structures.

The 7,147 raw shapes at the dominant start support a useful initial series up to
several thousand shapes without inventing alignments. Count actual retained
structures after StringTie filtering; do not assume every selected raw shape survives.

Duplicating existing reads changes N but eventually leaves K fixed. This can test
duplicate lookup cost, particularly with interleaved structures, but cannot by
itself establish scaling as the number of novel structures grows. Test N and K
separately. Repetition under different names does not create structural diversity.

Use the smaller datasets for guide, expression-only, and mixed-mode equivalence.
Do not combine unrelated species/contigs into one apparent biological locus.

### 3. Extend with deliberately synthetic structures if needed

Use an exon scaffold informed by the real records, then generate controlled
combinations or endpoint variants to grow the number of distinct structures while
holding starts fixed. Separate this structural stress case from tests of realistic
sequence errors and polyA/T evidence.

Synthetic records must be valid SAM/BAM: query-consuming CIGAR length must agree
with sequence/quality lengths, reference-consuming coordinates must fit the header,
and tags must remain consistent with the generated alignment. Rebuild or remove
invalidated alignment-derived tags rather than copying them blindly. A coherent
synthetic reference and sequences are preferable if reference sequence is involved.

Changing only QNAME, sequence errors, or interior CIGAR details that leave exon
coordinates and relevant deletion offsets unchanged does not reliably create new
collapse candidates. Conversely, arbitrary splice-site proliferation can make
graph inference dominate; keep junction diversity controlled and measure producer
and worker time separately.

Synthetic data can prove algorithmic scaling and exercise semantic edge cases.
It cannot establish that the reporter's complete workload is fixed, nor validate
the biological explanation of the immunoglobulin hotspots.

## Proposed validation and profiling

### Measure the producer directly

Record elapsed and CPU time, peak memory, bundle handoff time, and worker inference
time. Attribute producer work among input decoding, guide evaluation, collapse
lookup, junction insertion, and polyA/T processing. Candidate-comparison counters
are useful because they expose algorithmic work independently of machine load.

Compare builds made with the same compiler and optimization settings, on identical
inputs and options, with repeated runs and controlled background load. Start with
direct BAM input; also reproduce the reporter's SAM pipe separately to distinguish
input-format overhead. Existing log gaps and old executable names are not substitutes
for these measurements.

The handoff warns that debug and release objects cannot safely be mixed and that
some header dependencies are incomplete. Preserve baseline binaries and use clean,
consistent builds for a later implementation experiment.

### Compare selected representatives, not only totals

A temporary diagnostic build could execute both the original scan and the indexed
lookup for each read, assert the same match/no-match and selected readlist index,
and apply accumulation only once. This would detect changed eligibility or ordering
before downstream differences become difficult to interpret. It is a correctness
tool, not a performance benchmark because it retains the slow scan.

Targeted cases should include:

- many same-start distinct structures and repeated structures in several orders;
- identical exons with different interior CIGAR operations;
- differing junction-deletion offsets, both with and without guide matches;
- incoming unknown strand and stored polyA/T-inferred strand;
- terminal-exon removal that changes the effective start, including returning to
  a start represented earlier in the bundle;
- mixed short/long-read records and unitig status where applicable;
- secondary/supplementary records, early filtering, and bundle reuse;
- changes in NH, count/evidence accumulation, and mate-pair linkage;
- hash collisions if the implementation uses a fingerprint plus exact checks;
- merge-mode bypass and multi-threaded bundle processing.

Final output checks must include transcript structures and abundance values, not
only counts. Normalize assembly output IDs/order when comparing thread counts, as
assembly-mode numbering can depend on completion order. Compare relevant guide,
gene-abundance, and count-table output when those options are exercised.

Run the existing assembly and merge tests after an implementation change, but do
not rely on their small fixtures as sufficient coverage for this bottleneck.

## Reproduction notes for the completed inspection

Representative read-only commands used:

```bash
rg --files --hidden -g '*.bam' -g '*.cram' -g '*.sam' \
  -g '!htslib/**' -g '!samtools-0.1.18/**' -g '!.git/**' \
  -g '!SuperReads_RNA/**' -g '!gclib/**' .
samtools view -H issue_476/bundle_inf.bam
samtools view -H debug_L320/bundle_c1.bam
samtools idxstats debug_L320/bundle_c1.bam
git show v2.2.3:rlink.cpp
```

The Python inspection scripts were executed inline and were not saved separately.
To reproduce the principal measurements, stream `samtools view FILE`, skip records
whose FLAG has any bit in 2308, and apply the length/shape definitions above. For
the inventory stop after 5,000 retained records; for `debug_L320` deeper sampling
stop at 50,000; for the other five deeper candidates scan to EOF. Group shapes by
`(RNAME, POS)`, not by POS alone. Preserve the distinction between a prefix sample
and a complete file scan in any subsequent report.

## Remaining questions

- How much current-master producer time does the natural issue #476 hotspot spend
  in candidate lookup, and how many candidates are actually compared?
- Which raw structures survive current filtering, trimming, and guide handling?
- Does Ida's implementation preserve the existing eligibility window, effective
  matching state, and newest-first representative choice?
- After lookup acceleration, which remaining producer costs become significant?
- Do the reporter's original inputs show the same improvement and full output
  equivalence under a controlled A/B comparison?

No additional data is needed to begin answering the first four questions locally.
The original issue #511 inputs or reporter collaboration are still needed for the
last one.
