# New merge default compared with StringTie 3.0.3

## Scope

This comparison uses the 20 issue #247 input GTFs in
`merge-dev/merge-dev/example_dataset/stringtie_GTFs/` and the supplied
`Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gtf` guide.
Both runs used four threads and otherwise used their own `--merge` defaults.

The baseline executable was a clean release build from commit `dae0eff`, whose
version is 3.0.3. The new executable was the release build from the `merge-dev`
working tree after the merge-policy implementation and terminal-boundary
regression fix. The source version string has not yet been changed from 3.0.3.

```sh
$BASELINE/stringtie --merge -p 4 -G "$GUIDE" -o v3.0.3.gtf $INPUTS
./stringtie --merge -p 4 -G "$GUIDE" -o new-default.gtf $INPUTS
./compare_merge_outputs.pl "$GUIDE" v3.0.3.gtf new-default.gtf $INPUTS
```

## Main result

The new default emits 128,214 more transcript records: 192,173 instead of
63,959, or 3.00 times the 3.0.3 total. This is expected from preserving more
observed structures and doing less synthetic extension-assembly. It is not an
increase in speculative synthesis: exact source structures increase sharply
while synthesized structures decrease.

```text
Metric                              3.0.3      New       Delta       Change
----------------------------------  ---------  --------  ----------  --------
Transcript records                    63,959   192,173    +128,214   +200.46%
Unique structures                     63,867   192,143    +128,276   +200.85%
SET records                             6,240    90,054     +83,814  +1343.17%
MET records                            57,719   102,119     +44,400    +76.92%
Unique MET intron chains               56,428    97,359     +40,931    +72.54%
Generated-ID records                   37,900   166,114    +128,214   +338.30%
Exact source structures                20,290   172,682    +152,392   +751.07%
Synthesized structures                 43,577    19,461     -24,116    -55.34%
Guide IDs with exact guide structure   12,101    26,059     +13,958   +115.35%
Guide IDs with changed structure       13,958         0     -13,958   -100.00%
Duplicate transcript IDs                    0         0           0          -
Generated duplicate-structure groups       0         0           0          -
```

The union of exact structures in the input GTFs and guide contains 668,759
structures. The output retains 3.03% of that union in 3.0.3 and 25.82% under the
new default. This denominator includes many structures that should still be
removed as duplicates, contained models, weakly supported paths, or barrier
violations, so it is a retention comparison rather than a sensitivity claim.

Synthesized structures account for 68.23% of the 3.0.3 unique output but 10.13%
of the new unique output. The 3.0.3 count includes changed guide-ID structures,
which are correctly classified as different from their source guide structure.

## Guide and duplicate behavior

All 26,059 guide transcript IDs occur in both outputs. In 3.0.3, 13,958 of
those IDs carry a changed structure: 13,948 have changed outer bounds and 24
have a changed intron chain. Some guides have both changes. Under the new
default every guide ID retains its exact input exon structure.

Neither executable emits a duplicate transcript ID on this dataset, confirming
that the same-ID reporting bug seen in 3.0.2 is already absent in 3.0.3. Exact
duplicate-structure groups fall from 84 to 21. Every duplicate group in both
outputs contains guide IDs only; no group contains a generated transcript. The
21 new-default groups represent genuinely identical guide structures retained
under distinct reference IDs. `--dedup-guides` is the opt-in policy for
canonicalizing those IDs.

## Structure-set overlap

```text
Comparison                       Count
-------------------------------  -------
Common unique structures          19,151
3.0.3-only unique structures      44,716
New-only unique structures       172,992
Jaccard overlap                     8.09%
```

The low overlap is a direct consequence of changing both sides of the old
behavior: the new defaults stop abundance-based loss (`-F 0`, `-T 0`, `-f 0`)
and also stop many terminal-overlap unions. The large SET increase is especially
consistent with preserving terminal SET evidence instead of using it to widen
a MET. The MET and unique-intron-chain increases show that the effect is not
limited to single-exon records.

## One-run resource observation

```text
Run          Wall time   Maximum RSS
-----------  ----------  -----------
3.0.3        15.37 s       87,380 KB
New default  16.51 s      181,460 KB
```

The new run was 7.4% slower and used 2.08 times the peak memory in this single
measurement. This is not a benchmark: process cache, temporary I/O, and run
order were not controlled. The memory increase is plausible because many more
candidate structures survive to final policy evaluation and output.

## Reproducibility

`compare_merge_outputs.pl` parses transcript/exon structures, counts SETs and
METs, compares guide IDs to their original structures, classifies output
structures as exact source or synthesized, and compares the two structure sets.
The raw comparison output for this run was saved temporarily as
`comparison.tsv`; generated GTFs and temporary build trees were not added to the
repository.
