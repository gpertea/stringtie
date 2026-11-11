# New merge default compared with StringTie 3.0.3

Numbers updated 2026-11-11 for the final state of the `merge-dev` branch
(`3e5a15e`, merged locally with master `18669d0`; version 3.0.4). An earlier
version of this note described the branch at `a27c05d`, before intron-chain
support, end election and the default `--min-support` were added; its numbers
no longer apply.

## Scope

This comparison uses the 20 issue #247 input GTFs in
`merge-dev/merge-dev/example_dataset/stringtie_GTFs/` and the supplied
`Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gtf` guide.
Both runs used their own `--merge` defaults. With 20 input files the new
default requires an intron chain in 2 files (`--min-support 2`).

The baseline executable is a release build of version 3.0.3.

```sh
$BASELINE/stringtie --merge -p 8 -G "$GUIDE" -o v3.0.3.gtf $INPUTS
./stringtie --merge -p 8 -G "$GUIDE" -o new-default.gtf $INPUTS
./compare_merge_outputs.pl "$GUIDE" v3.0.3.gtf new-default.gtf $INPUTS
```

## Main result

The new default emits 7,325 more transcript records: 71,284 instead of 63,959
(+11.5%). The composition changes more than the total: structures that occur
in an input file or in the guide rise from 20,290 to 53,949, while structures
found in no input fall from 43,577 to 17,305.

```text
Metric                              3.0.3      New       Delta       Change
----------------------------------  ---------  --------  ----------  --------
Transcript records                    63,959    71,284      +7,325   +11.45%
Unique structures                     63,867    71,254      +7,387   +11.57%
SET records                             6,240    19,739     +13,499  +216.33%
MET records                            57,719    51,545      -6,174   -10.70%
Unique MET intron chains               56,428    50,245      -6,183   -10.96%
Generated-ID records                   37,900    45,225      +7,325   +19.33%
Exact source structures                20,290    53,949     +33,659  +165.89%
Synthesized structures                 43,577    17,305     -26,272   -60.29%
Guide IDs with exact guide structure   12,101    26,059     +13,958  +115.35%
Guide IDs with changed structure       13,958         0     -13,958  -100.00%
Duplicate transcript IDs                    0         0           0          -
Generated duplicate-structure groups       0         0           0          -
```

The union of exact structures in the input GTFs and guide contains 668,759
structures. The output retains 3.03% of that union in 3.0.3 and 8.07% under the
new default. This denominator includes many structures that should still be
removed as duplicates, contained models, weakly supported paths, or barrier
violations, so it is a retention comparison rather than a sensitivity claim.

Synthesized structures account for 68.23% of the 3.0.3 unique output and 24.29%
of the new unique output. The 3.0.3 count includes changed guide-ID structures,
which are correctly classified as different from their source guide structure.

## Guide and duplicate behavior

All 26,059 guide transcript IDs occur in both outputs. In 3.0.3, 13,958 of
those IDs carry a changed structure: 13,948 have changed outer bounds and 24
have a changed intron chain. Some guides have both changes. Under the new
default every guide ID retains its exact input exon structure.

Neither executable emits a duplicate transcript ID on this dataset. Exact
duplicate-structure groups fall from 84 to 21. Every duplicate group in both
outputs contains guide IDs only; no group contains a generated transcript. The
21 new-default groups are identical guide structures kept under distinct
reference IDs. `--dedup-guides` is the opt-in policy for collapsing those IDs.

## Structure-set overlap

```text
Comparison                       Count
-------------------------------  -------
Common unique structures          15,315
3.0.3-only unique structures      48,552
New-only unique structures        55,939
Jaccard overlap                    12.78%
```

The low overlap follows from changing both sides of the old behavior: the new
defaults stop abundance-based loss (`-F 0`, `-T 0`, `-f 0`) and keep observed
structures, and they stop many terminal-overlap unions and single-file
extensions. The SET increase is consistent with keeping terminal SET evidence
instead of using it to widen a MET.

## Human data

On 10 human brain samples per library type (assembled with 3.0.3, no `-G`),
merged with default settings and scored with gffcompare against GENCODE 49:

```text
Assemblies (10 samples)     Matching intron chains   Intron chain precision
                            3.0.3      New           3.0.3      New
--------------------------  ---------  ------------  ---------  ----------
polyA (unstranded)             20,360        33,021      49.5%       49.3%
RiboZero (--rf)                15,964        27,819      42.9%       44.6%
RiboZero (--rf -N)             15,612        28,076      48.5%       47.1%
```

## One-run resource observation

```text
Run          Wall time   Maximum RSS
-----------  ----------  -----------
3.0.3        15.53 s       87,328 KB
New default  16.31 s      110,272 KB
```

Alga dataset, `-p 4`, one run each. This is not a benchmark: process cache,
temporary I/O, and run order were not controlled.

## Reproducibility

`compare_merge_outputs.pl` parses transcript/exon structures, counts SETs and
METs, compares guide IDs to their original structures, classifies output
structures as exact source or synthesized, and compares the two structure sets.
Generated GTFs and temporary build trees were not added to the repository.
