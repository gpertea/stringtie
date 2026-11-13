# Issue 511 follow-up: leads worth keeping and the -falign-functions=64 option

Date: 2026-11-12. Author: Claude (Fable 5.2), forked session in
`/home/gpertea/src/stringtie`.

Context: this continues `worklogs/fable-i511-review-and-test.md`. After that
report, `ghash-clear` was merged into `i511-long-reads` on GitHub (`fec648b`). This
note records what came out of the later layout and timing experiments: three
things worth keeping, and the case for and against one compiler flag. All timings
are wall seconds on `mix/chr21_short.bam` (17.4M short reads), default threads,
Ryzen 9 5900X, `g++ 13.3 -g -O3`.

Data: `issue_511/run-index/layout_alignment_followup.tsv`,
`issue_511/run-index/layout_ghash_clear_timing.tsv`.

## 1. Worth keeping

### 1.1 Possible producer speedup: per-read tag lookups (not tried)

In a callgrind profile of master on the first 1.5M records of the file,
`bam_aux_get` accounts for 3.84e9 of 1.91e10 instructions, 20%. It is the largest
single function in the profile.

Each read triggers six or more separate lookups, and each one scans the read's
tag block from its start:

```text
Where                  Tags
---------------------  ----------------------------------------------
stringtie.cpp:859      XS, then ts if XS is absent (spliceStrand)
stringtie.cpp:915,917  NH, HI
rlink.cpp:754,756      NM, YK
rlink.cpp:911          YC
rlink.cpp:923          nM, when NM is 0
rlink.cpp:939          MD, for some short reads with mismatches
```

Tags that are usually absent (`YK`, `YC`, `nM`) cost a scan of the whole block.
Reading the needed tags in one pass over the block would remove most of this.

Limits of the evidence: this is an instruction share, not time, on one short-read
subset, and the profiled input was SAM text. The share of `bam_aux_get` does not
depend on the input format, but the total does. No change was written or timed.

### 1.2 Benchmarking rule for this machine

Any binary occasionally runs slow on this file: 28-37 s instead of 24-26 s. It
happens with address randomization on or off, and I did not find the cause. Seen
on master, on the run-index build and on aligned builds.

Consequences:

- A timing claim needs at least three runs per binary.
- Earlier single-run comparisons were affected. One claim in section 6 of the
  main report does not hold: "master relinked with 8 bytes of padding runs in
  28.0 s". The follow-up subsection there says so; the original sentence is still
  in the text and should be removed.
- The stack position is not a factor: one unchanged binary with the environment
  padded by 0 to 128 bytes, randomization off, stayed at 24.6-26.1 s.

### 1.3 Residual layout effect on the merged branch

With `ghash-clear` in and the default flags, 11 link layouts of the same objects
(padding of 0 to 56 bytes before `stringtie.o`) differ by up to 8%: the fastest
layout runs in 12.8-13.2 s and the slowest in 13.8-14.3 s, five runs each. Which
loop is responsible was not identified. Before `ghash-clear` the same effect was
about 25% and came from the `hashread.Clear()` loop.

Practical meaning: on this file a difference under about 1 s between two builds is
not evidence of anything.

### 1.4 Not tested

Profile-guided optimization, `-falign-loops`, and a non-PIE link.

## 2. The -falign-functions=64 option

### What it does here

The slowdown that made the run-index build look 25% slower than master was a
placement effect: `hashread.Clear()` is inlined in `main()`, its loop runs 7
instructions per empty bucket with two taken branches, and its speed changed with
where the linker put it. That build was slow on every run, with address
randomization on or off. Without `perf` I could not confirm the branch-prediction
mechanism directly.

`-falign-functions=64` starts every function on a 64-byte boundary, so code added
in another function no longer shifts the loop.

```text
Build                      Flags                  Wall
-------------------------  ---------------------  ---------------------------
master                     default                24.43 s mean, 16 runs
master                     -falign-functions=64   25.32 s mean, 16 runs
run index (before merge)   default                30.8-31.0 s, every run
run index (before merge)   -falign-functions=64   25.5-26.1 s, two layouts
ghash-clear                default                13.68 s mean, 21 runs
ghash-clear                -falign-functions=64   13.71 s mean, 15 runs
```

- The 3.7% cost on master is specific to master: the `Clear()` loop is pinned at
  a position that is slower than where the default layout happens to put it.
- On the `ghash-clear` code, which is what `i511-long-reads` now contains, the flag
  costs nothing measurable, and five layouts gave means of 13.4-13.9 s that cannot
  be separated from run-to-run noise (against the 8% of section 1.3).

### Downsides

Measured:

- Size: master's code section grows from 1,735,298 to 1,759,490 bytes, 1.4%.

Not measured, from how the flag works:

- It fixes the layout; it does not choose a good one. A hot loop can end up pinned
  at a slow position, as on master.
- It protects only against changes in other functions. Editing a function, or
  code inlined into it, still moves the loops inside it.
- htslib is built by its own script with its own flags, so BAM decoding code is
  not aligned by it.
- Code is spread out slightly more, which can cost some instruction-cache and TLB
  efficiency in call-heavy code. No sign of that on this file.
- The effect chased here depends on how this processor handles those branches;
  other CPUs may behave differently.

### Recommendation

Add it to the release flags. On the merged branch it cost nothing measurable, and
it keeps unrelated code changes from showing up as speed changes. It was tested on
one short-read file only; a long-read run and a `--mix` run with and without the
flag would settle it. No `Makefile` change has been made.
