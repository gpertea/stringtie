# StringTie `--merge` issues

## Problem statement

[Issue #247](https://github.com/gpertea/stringtie/issues/247) reports that
`stringtie --merge`, when combining sample assemblies with a guide annotation,
can change known reference transcript boundaries and join models into neighboring
gene regions. The immediate goal was to inspect the supplied example, identify
the main weakness behind that result, and determine whether it is an isolated
UTR-extension bug or evidence of a broader merge-mode problem.

## Main weakness

`--merge` looks like a conservative catalog union, but internally it performs
graph-based transcript reassembly. It does not guarantee that an input model is
preserved exactly, or that a transcript supplied with `-G` remains structurally
unchanged. Complete transcript models are projected onto splice graphs,
compatible paths are combined, endpoints may be extended, and reference IDs may
be attached to reconstructed paths.

This mismatch between expected and actual semantics is the central problem. The
[issue discussion](https://github.com/gpertea/stringtie/issues/247) and the
[summary comment](https://github.com/gpertea/stringtie/issues/247#issuecomment-5076112338)
describe it as a long-standing limitation rather than a single threshold bug.

## What the supplied example shows

For `Ec-04_004450.1`:

```text
Reference: chr_04:4570279-4600321, 31 exons
Merged:    chr_04:4570279-4601396, 30 exons
Change:    right endpoint extended by 1,075 bp
```

Two source models provide the longer downstream endpoint:

```text
SRR12897873.7658.2  4590729-4601396
SRR12897873.7658.3  4594733-4601396
```

They share the downstream splice chain and terminal exon ending at `4601396`.
Merge mode transfers that endpoint onto the model carrying the reference ID.
The result:

- overlaps the next reference gene, `Ec-04_004460.1`, which starts at `4601368`
  on the opposite strand;
- removes a tiny internal intron: reference exons `4584961-4587464` and
  `4587466-4588433` become one exon, `4584961-4588433`;
- retains transcript ID `Ec-04_004450.1` although its endpoint and exon chain no
  longer match the reference transcript.

The example is documented in the
[author's dataset comment](https://github.com/gpertea/stringtie/issues/247#issuecomment-5084881023).

## How large the problem is in this dataset

A full comparison of the supplied reference and merged output gave:

```text
Reference transcripts                         26,059
Reference IDs found in merged output           25,904
Reference IDs with changed outer coordinates   13,948  (53.8%)
Reference IDs with changed exon structures     13,958  (53.9%)
Extensions greater than 1 kb                    2,120
Largest extension                               8,929 bp
```

Every changed outer boundary is an extension of the right genomic coordinate.
That matches the asymmetric collapse in
[`build_merge()`](../rlink.cpp#L16022), which retains the earlier model and
expands its right endpoint to the maximum. Ten models retain their outer bounds
but have altered exon chains, mostly because very short introns are merged.
Guide loading enables this close-exon merging in
[`stringtie.cpp`](../stringtie.cpp#L301).

## Why this is systemic

1. Gene identity is discarded during reconstruction. GTF models are converted
   to synthetic alignment-like records, then grouped by overlap, strand, and
   junction graph rather than protected by their original gene boundaries.
2. Models with the same junction chain are collapsed asymmetrically. Extending
   the right coordinate is not strand-aware: it changes the 3-prime end on `+`
   transcripts and the 5-prime end on `-` transcripts.
3. Guides are graph evidence, not immutable records.
   [`merge_transfrags()`](../rlink.cpp#L7534) and
   [`store_merge_prediction()`](../rlink.cpp#L7142) reconstruct exons from graph
   nodes rather than copying the original `GffObj`. If a node is enlarged, the
   reference-associated prediction is enlarged too.
4. Compatibility can be transitive. If A overlaps B and B overlaps C, merge can
   reconstruct a larger path even when no input contained that complete model.
5. Final output performs another structural merge.
   [`printMergeResults()`](../rlink.cpp#L19613) removes weaker overlapping models
   and can join nearby single-exon models.

Increasing junction-support thresholds does not prevent the reported endpoint
extension because extending a terminal exon introduces no new junction.

## Additional behavior in StringTie 3.0.2

Rerunning the supplied inputs with the v3.0.2 binary reproduces the altered
model. The full run also emits 25,226 reference IDs twice: once as a graph-
reconstructed StringTie prediction and once as a later guide-source record. At
the example locus, the same transcript ID therefore appears with different
outer coordinates. Even the guide-source copy has already lost the short
internal intron during guide loading.

This duplicate-ID behavior is a separate identifier-integrity problem that was
not present in the older merged output supplied with the issue.

## Bottom line

The iceberg is not one overextended UTR. Users commonly expect union,
deduplication, and immutable guide models; the implementation instead performs
reassembly and synthesis. A conservative repair needs an explicit contract:

- never modify a `-G` transcript;
- never reuse a reference ID for a changed structure;
- preserve input exon chains unless a separate normalization option is enabled;
- prevent transitive synthesis across gene boundaries;
- make endpoint reconciliation optional rather than implicit.

## StringTie 3.0.3 follow-up

The supplied inputs were rerun with the newly built StringTie 3.0.3 binary,
using the original merge command and default thresholds. This version fixes the
duplicate transcript-ID bug, but it does not fix guide rewriting or all
duplicate-structure output.

```text
Metric                                           Result
------------------------------------------------ -------
Output transcript records                       63,959
Unique transcript IDs                            63,959
Duplicate transcript IDs                              0
Reference IDs emitted                            26,059 / 26,059
Reference IDs with changed exon structure         13,958
Reference IDs with changed outer coordinate       13,948
Largest endpoint extension                         8,929 bp
Exact duplicate output-structure groups                84
Extra records from duplicate structures                92
```

`Ec-04_004450.1` is still changed from 31 to 30 exons, and its right endpoint is
still extended from `4600321` to `4601396`.

The transcript-ID duplication was fixed by marking reconstructed guides as
stored in [`printMergeResults()`](../rlink.cpp#L19613), preventing
`writeUnbundledGuides()` from printing them again. Exact structure duplication
remains:

- all 84 duplicate groups contain different reference transcript IDs;
- 16 groups represent structures already duplicated in the reference;
- 68 groups are newly created when distinct reference isoforms are rewritten
  into the same output structure.

Thus v3.0.3 fixes duplicate IDs, but not duplicate structures caused by guide
mutation.

### Effects of `-F`, `-T`, and `-f`

Controlled reruns produced:

```text
Options                              Output transcripts
-----------------------------------  ------------------
default: -F 1 -T 1 -f 0.01                       63,959
-F 0                                               89,957
-T 0                                               63,959
-f 0                                               66,263
-F 0 -T 0                                        111,828
-F 0 -T 0 -f 0                                   122,763
```

All 927,737 input transcripts contain both FPKM and TPM attributes:

```text
Input transcripts                                 927,737
Rejected by default -F 1                          476,373
Rejected by default -T 1                          252,348
Rejected only by -T 1                                   0
Passing both defaults                             451,364
```

For this dataset, `-T 0` produces the same structure set as the default. The
FPKM threshold is the controlling absolute filter.

[`processRead()`](../rlink.cpp#L766) rejects an input transcript when any
available metric fails its threshold: coverage below `-c`, FPKM below `-F`, or
TPM below `-T`. Missing attributes retain the sentinel value `-1` and are not
filtered. Therefore `-F` and `-T` are meaningful mainly for input GTFs that
carry those attributes, normally StringTie output. When both attributes exist,
the transcript must pass both thresholds.

In the active merge path, `-f` is primarily applied by
[`printMergeResults()`](../rlink.cpp#L19613) as a final relative-abundance
filter. It:

- compares overlapping, same-strand predictions;
- removes a non-guide multi-exon prediction when its reconstructed coverage is
  less than `-f` times another prediction's coverage;
- protects guide-associated predictions;
- uses stricter logic for weaker overlapping single-exon predictions that is
  independent of `-f`.

This is a pairwise ratio, not a fraction of total locus abundance, and the
comparison is not restricted to transcripts from the same original gene. The
`good_merge_junc()` use of `-f` is inactive because its merge-mode call is
commented out.

### Suitability of normalized abundance for cross-sample merging

FPKM and TPM can be optional noise heuristics, but they are not statistically
calibrated evidence for cross-sample structural merging. When TPM is available,
[`build_merge()`](../rlink.cpp#L16022) uses it as transcript abundance ahead of
FPKM and coverage. [`store_merge_prediction()`](../rlink.cpp#L7142) combines
those values into reconstructed prediction coverage, which `-f` then compares.

Consequences include:

- TPM values from different samples are combined as if directly equivalent;
- library composition can change TPM without changing structural support;
- replicate count and sample-level support are not modeled;
- one highly expressed sample can dominate a structure found in other samples;
- a valid low-abundance splice chain can be discarded by `-f`;
- disabling the filters retains more structures but supplies more endpoint
  evidence to the graph, increasing changed reference structures from 13,958
  to 16,771 in this dataset.

There is also a mixed-input hazard: abundance availability is stored
bundle-wide. If any retained transcript supplies TPM, all non-guide transcripts
in that bundle are weighted through their TPM field; a transcript lacking TPM
retains the sentinel value `-1`.

The abundance filters are a major cause of disappearing input structures, but
disabling them does not solve merge correctness. It preserves more evidence
while increasing graph-driven endpoint extension and guide mutation.

## Merge-policy implementation and issue-dataset validation

The merge redesign now uses structural support by distinct input file. Merge
defaults are `-F 0`, `-T 0`, and `-f 0`; explicit abundance filters remain
available and a structure filtered in every input can be restored by replicated
support. The default path policy requires one exact shared intron, applies
independent biological 5-prime/3-prime extension evidence, rejects guide-gene
and coding-transcript bridges, and preserves every `-G` structure and ID.

The complete issue #247 dataset was rerun with all 20 input GTFs:

```text
Metric                                      New default
------------------------------------------  -----------
Reference transcripts                            26,059
Missing reference IDs                                 0
Changed reference structures                          0
Output transcripts                              144,781
Duplicate transcript IDs                              0
Duplicate generated structures                        0
Duplicate guide-structure groups                      21
```

The remaining duplicate structures under the default are reference models with
different reference IDs, retained intentionally. With `--dedup-guides`, the
lexicographically smallest reference ID is kept, aliases are reported, and the
same full run has zero duplicate structure groups. `Ec-04_004450.1` is emitted
with its original 31-exon structure and original `4570279-4600321` bounds.
