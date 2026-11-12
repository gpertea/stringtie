# Issue 511: review of the same-start search, a behavior-preserving fix, and tests

Date: 2026-11-12. Author: Claude (Fable 5.2), session in `/home/gpertea/src/stringtie`.

Revision examined: `master` at `53d9647` (StringTie 3.0.4, after the merge-dev merge).
Baseline binary for all comparisons: a release build of that revision
(`g++ 13.3 -g -O3`), on a Ryzen 9 5900X.

Companion documents: `worklogs/astra-i511-long-read-feasibility-testing.md` (Codex
survey of the source and of local long-read data) and
`worklogs/fable-junction-deletions-issue-draft.md` (issue text for section 8).

## 1. Summary

- Ida's diagnosis of issue #511 is right: the duplicate search in `processRead()`
  (`rlink.cpp:761`) is quadratic in the number of distinct alignments sharing one
  start. On a synthetic locus it takes 85-97% of the run time.
- "Replace it with a hash lookup" is not a drop-in change. Which reads get
  collapsed affects the output values, the search only covers part of the reads
  with that start, and the stored read does not have the same key as the alignment
  that created it. A hash that collapses differently changes the output.
- Branch `i511-long-reads` (one commit, `9a398c5`, local, not pushed) indexes long
  same-start runs and returns exactly the entry the scan returns. On 19 runs over
  local data (up to 17.5M reads) a per-read check found no difference and every GTF
  is byte-identical to master at `-p 1`. Synthetic hot locus: 65.8 s to 1.8 s at
  200k reads; at 640k reads master needs 1,726 s and the branch 5 s.
- While building test cases I found four defects in the handling of deletions next
  to introns (long reads and `--mix`). They are separate from #511 and change
  results. Branch `lr-junction-deletions` (four commits, local) has proposed fixes;
  `issue_511/junction-deletions/` has reproducers and a regression script. The
  issue text is drafted but not filed (section 8).
- Hash functions (section 7): wyhash, libstdc++ `std::hash<std::string>`, FNV-1a,
  a boost-style `hash_combine` and a `std::unordered_map` version all give the
  same output and the same speed within noise. `hash_combine` produces full 64-bit
  collisions (6,050 on 340k keys), which is harmless here only because every hit
  is confirmed by the original comparison.
- A separate producer cost turned up while timing (section 9): the mate hash is
  cleared at every bundle and `GHash::Clear()` scans the whole table at its
  largest size. After one 14.6M-read bundle every later bundle costs 4 ms; on a
  chr21 short-read sample that is 14 of 24 s. Branch `ghash-clear` (one commit,
  `6e8ff08`, local) fixes `Clear()` and iteration in `gclib/GHashMap.hh`: 24.5 s
  to 13.6 s, identical output.

## 2. Context

### The issue

[#511](https://github.com/gpertea/stringtie/issues/511): with `-L -G` on a chr2
cDNA BAM (927,114 records), StringTie 2.2.2/2.2.3 spends 322 s/513 s at 100% of one
CPU before the bundle `chr2:88,538,720-89,330,679` (634,090 alignments, 337,129
distinct, 12,859 junctions) "begins processing"; the bundle itself then takes 1-2 s.
A similar gap precedes a chr22 bundle. The reporter suspected the same-start search
in `processRead()` and asked for a fix that keeps the assembly rules.

Ida confirmed the hypothesis (IGK/IGL loci: very many reads sharing a start, most of
them distinct structures) and said she is replacing the linear search with a hash
lookup. The reporter then reported a local "same-start candidate index retaining the
original matching checks": chr2 459.7 s to 42.2 s, chr22 95.5 s to 33.9 s, output
equivalence not yet checked.

The question for this session was whether the fix is as simple as that.

### The code

The main thread is the only producer: it reads sorted alignments, calls
`BundleData::evalReadAln()` and `processRead()` for each, and queues the finished
bundle for a worker. "begins processing" is printed by the worker, so a gap before
it is producer time.

`processRead()` collapses identical alignments into one `CReadAln`:

```cpp
int n=readlist.Count()-1;
if(!mergeMode) while(n>-1 && readlist[n]->start==brec.start) {
	if(strand==readlist[n]->strand && (readlist[n]->longread==longr) && (!isunitig || (unitig_cov>0) == readlist[n]->unitig)) {
		match=exonmatch(readlist[n]->segs,brec.exons);
		if(match && (longreads || mixedMode)) match=deljuncmatch(readlist[n],brec.juncsdel); //DEL AWARE
	}
	if(match) break;
	n--;
}
```

`readlist` is unsorted during production (`readlist(false,true)`, `bundle.h:413`);
new reads are appended. The loop walks back from the end while the stored start
equals the incoming start, and stops at the first match. With K distinct alignments
at one start, K*(K-1)/2 comparisons fail before the entries are all stored, and
every later duplicate scans until it meets its own entry.

## 3. Diagnosis, measured

### Synthetic reproducer

`issue_511/synthetic/gen.py N 0.7` writes N minus-strand reads that all start at
`chrT:100000`: a fixed first exon, then one of 5 "J" segments with a variable
boundary, one of 40 "V" segments, and a variable read end. About 53% of the reads
are distinct structures, close to the reporter's 337,129 of 634,090.

```text
Reads     Distinct   3.0.3 wall   3.0.4 wall   Scan loop     Comparisons
-------   --------   ----------   ----------   -----------   -----------
 25,000     13,363      0.83 s         -            -              -
 50,000     26,592      2.74 s       2.69 s     2.29 s (85%)     5.8e8
100,000     53,254     10.78 s      10.73 s     9.85 s (92%)     2.3e9
200,000    106,278     63.49 s      66.22 s    64.15 s (97%)     9.4e9
```

The scan-loop time and comparison counts come from a build of master with a timer
and counters around the loop (`perf` is not usable on this machine,
`perf_event_paranoid`). The cost per comparison grows from 3.9 ns to 6.8 ns as the
run outgrows the cache. The rest of the producer costs about 10 us per read.

The reporter's 322 s is consistent with this: at 7-10 ns per comparison it is
3-5e10 comparisons, about a quarter of what 634k reads would need if all 337k
distinct structures shared one start. I have not seen their data; the distribution
over starts is unknown.

### Local data

Same instrumented build, real data (all columns from master):

```text
Run            Options       Reads        Distinct    Longest run   Comparisons   Scan time
-------------  -----------  -----------  ----------  ------------  ------------  ---------
i476_L         -L               173,651      33,280        7,964         4.9e8     1.52 s
i476_full      -L             7,637,422   2,989,168        8,300         1.8e9     4.96 s
L320_LG        -L -c 2 -G     1,328,115     899,345        7,415         1.5e8     0.56 s
c1_L           -L               653,732     456,123        2,963         7.6e7     0.40 s
c1_LG          -L -G            663,430     469,461        4,245         5.1e7     0.35 s
b357_bundle2   -L             1,398,346     886,990        2,616         1.1e8     0.32 s
i465_L         -L             7,771,227   1,016,186          677         4.6e7     0.25 s
b357_bundle    -L               614,381     430,397        2,462         3.7e7     0.18 s
c21_mixG       --mix -G      17,466,391   1,043,713          925         3.2e8     0.83 s
c21_S          (short)       17,393,855     988,465           67         2.8e8     0.72 s
mixbug_mixG    --mix -G       1,477,889     131,053          567         1.4e7     0.06 s
```

Inputs: `issue_476/bundle_inf.bam` and `issue_476/4801.heart...sorted.bam`,
`debug_L320/bundle_c1.bam` with `chr.gtf`, `mix/chr1_lr.bam`, `mix/chr21_*.bam`,
`bugL357/bundle*.bam`, `issue_465/bundle2_long_c10A.bam`, `mixbug/chr1.*.bam`.

This answers the first open question of the Codex document: its natural hotspot
(`issue_476/bundle_inf.bam`, 75,204 reads at one start) costs master 4.9e8
comparisons and 1.5 s, with a longest run of 7,964 stored reads. None of the local
files comes near the reporter's case; the scan is a small part of a normal run.

I found no other quadratic step in the producer. After the fix, 640k synthetic
reads take 5 s in total.

## 4. Why a plain hash lookup changes behavior

### 4.1 Collapsing is not neutral for the output

Two variants of master, chr21 data, `-p 1`, compared with master:

```text
Data, options               Variant                    Structures     GTF lines differing
--------------------------  -------------------------  -------------  -------------------
chr21 long, -L              no collapsing              same                769 of 5,277
                            scan capped at 64 entries  same                635 of 5,277
chr21 long, -L -G           no collapsing              same                617 of 4,867
                            scan capped at 64 entries  same                460 of 4,867
chr21 short                 no collapsing              60 lines differ   5,917 of 6,051
                            scan capped at 64 entries  same                943 of 6,051
chr21 short + long, --mix   no collapsing              34 lines differ   2,955 of 6,863
                            scan capped at 64 entries  same                899 of 6,863
```

A control (master twice) differs on 0 lines. So a fix that collapses even slightly
more or less than the scan changes cov/FPKM/TPM values, and can change structures.
A fix intended as a pure speedup has to return the same entry as the scan for
every read.

### 4.2 The search covers only the tail of the list

The loop stops at the first stored read with a different start. Starts are not
monotonic: when a long read has its first exon removed as a polyT artifact,
`brec.start` is moved to the next exon (`rlink.cpp:748`) and the read is stored
under a later start. It then ends the run for the reads that follow at the earlier
position, and is itself not found when the input reaches its start.

An audit build (a bundle-wide exact lookup next to the scan) counted the reads for
which the scan found nothing although an identical stored read existed earlier in
the bundle: 1 on `mix/chr21_long.bam`, 2 on `mix/chr1_lr.bam`, 2 on
`issue_465/bundle2_long_c10A.bam`, 0 on `mixbug/chr1.long.bam`. Rare, but a
bundle-wide hash would collapse those reads and the output would no longer be
identical.

### 4.3 The stored read is not keyed like the alignment that created it

- Strand: an unstranded long read gets a strand from polyA/T evidence when it is
  stored (`rlink.cpp:825`). The next identical unstranded alignment compares strand
  0 with the stored 1 or -1 and does not match.
- Unitig flag: set after the read is stored (`rlink.cpp:915`).
- Deletion offsets: `deljuncmatch()` derives them from the stored junction, which
  is not adjusted when the junction matches a guide (`rlink.cpp:851`).
- More than one stored read can match (see the cases above); the scan returns the
  newest.

A hash keyed on the incoming alignment on both sides gets these wrong. The key of a
stored read has to be computed from the stored read.

### 4.4 Some identical alignments never match

Two kinds of alignments are never collapsed on master (section 8, defects 2 and 3):
a deletion right after an intron, and with `-G` a deletion next to a guide-matched
junction. They are 12% of the spliced alignments in `mix/chr1_lr.bam`. They make
runs longer, and a hash fix must reproduce the non-matches (or the output changes).

### 4.5 The common case is already cheap

On normal data the scan costs 0.2-0.8 s per file (table in section 3), and most
reads match one of the newest few entries. Hashing every read cannot gain anything
there and can lose; the index should only be used for long runs.

## 5. The run index (branch `i511-long-reads`)

One commit on top of master, `rlink.cpp` only, +177/-7 lines. No change to headers,
options or other modes; `gclib/wyhash.h` (already included through `GHashMap.hh`)
provides the hash.

- A "run" is the set of `readlist` entries at the end of the list having the start
  of the incoming alignment, i.e. exactly what the loop visits. A static
  `CRunIndex runidx` tracks its length. The run is reset when the last stored read
  has a different start (which is also what happens at the first read of a bundle).
- While the run has fewer than `RUNIDX_MIN` (64) entries, the original loop runs.
- At 64 entries the run is indexed: for each entry a 64-bit wyhash of the fields
  the comparison uses, taken from the stored read (strand, long/short, unitig flag
  if `-U`, exon coordinates, and in long-read or mixed mode the junction count and
  the deletion offsets as `deljuncmatch()` derives them). Entries are chained per
  bucket, newest first.
- A lookup hashes the same fields of the incoming alignment (its own strand, its
  `juncsdel` entries and their count), walks one bucket, and confirms every
  candidate with the original predicate (`sameReadAln()`, which is the old `if`
  body). The first confirmed candidate is the newest match, as in the scan. A hash
  collision costs one extra comparison and cannot change the result.
- A new read is added to the index at the end of `processRead()`, when its strand
  and unitig flag are final.
- The index holds positions within the run, not pointers, and is only touched by
  the producer thread, like the existing static `_id`. Merge mode does not use it.
- Memory: 12 bytes per run entry plus 4 bytes per bucket (at most 2x the entries),
  kept and reused. 640k reads at one start: 148 MB peak RSS in total, against
  51 MB for 200k on master (53 MB with the index).
- `-DRUNIDX_CHECK` builds a binary that also runs the plain search for every read
  and stops with an error if the result or the run length differs.

Why the key equality is exact: the incoming key and the stored key are equal if and
only if the predicate accepts the pair. The junction count is part of the key, so
an alignment with a doubled `juncsdel` entry (defect 2) hashes differently from
every stored read and misses, as it does on master.

## 6. Validation

### Per-read equivalence

`RUNIDX_CHECK` builds with thresholds 64 and 2 (the latter indexes almost every
run), `-p 4`. Zero mismatches in all of:

```text
Run            Options                    Reads         GTF vs master (-p 1)
-------------  -------------------------  ------------  --------------------
i476_L         -L                              173,651  identical
i476_full      -L                            7,637,422  identical
i357_eLG       -e -L -G                            422  identical
i357_L         -L                                  481  identical
rg_mix         --mix                           244,724  identical
mixbug_mixG    --mix -G                      1,477,889  identical
mixbug_LG      -L -G                            22,975  identical
b357_long      -L                                4,005  identical
b357_bundle    -L                              614,381  identical
b357_bundle2   -L                            1,398,346  identical
c21_L          -L                               71,708  identical
c21_LG         -L -G                            72,536  identical
c21_mixG       --mix -G                     17,466,391  identical
c21_S          (short reads)                17,393,855  identical
c1_L           -L                              653,732  identical
c1_LG          -L -G                           663,430  identical
i465_L         -L                            7,771,227  identical
L320_LG        -L -c 2 -G                    1,328,115  identical
syn200k        -L (synthetic)                  200,000  identical (no transcripts)
```

"identical" is `cmp` of the whole GTF except the `# ` header lines, both binaries
at default threads. The earlier prototype was also identical on a 100k synthetic
file that does produce transcripts (196 lines).

Negative control: a build in which one stored hash in seven is corrupted stops at
once with `RUNIDX_CHECK: search mismatch`, so the check can fail.

`make test` on the branch: 11 assembly tests OK, merge policy test OK.

The two deletion cases of section 4.4 were run through the check as well
(`t_ND.sam`, `t_DN.sam -G`): the non-matches are reproduced.

### Speed

Synthetic hot locus, `-L -p 2`:

```text
Reads     Master     Branch    Peak RSS master / branch
-------   --------   -------   ------------------------
100,000    11.0 s    0.84 s    30 MB / 31 MB
200,000    65.8 s    1.80 s    51 MB / 53 MB
640,000   1,726 s    4.97 s    140 MB / 148 MB
```

The 640k master run produces the same GTF as the indexed run, but that file has
no transcripts (the synthetic locus yields none at that depth), so it only shows
that both runs complete alike. The per-read check was run at 200k, not at 640k.

Real data, master against branch, `-p 4`, best of two alternating runs:

```text
Run          Master wall   Branch wall
-----------  -----------   -----------
i476_L           4.87 s        3.20 s
i476_full      116.43 s      112.24 s
c1_L             9.73 s        9.44 s
c1_LG            9.48 s        9.18 s
b357_bundle      7.73 s        7.59 s
L320_LG         22.00 s       21.53 s
i465_L          33.60 s       34.22 s
c21_S           24.75 s       31.17 s   <- not the index, see below
```

The long-read runs are equal or faster; `i465_L` (+1.8%) is within what I can
resolve with two runs.

### The short-read result is a code layout effect

`c21_S` (17.4M short reads) was reproducibly 25% slower with the branch. The index
is hardly used there (152,017 lookups, 2 index builds), so I traced it:

- The extra time is per bundle, not per read. After the bundle
  `chr21:7609049-8472236` (14.6M reads) master hands off about 240 bundles per
  second and the branch about 170, even for bundles of a few reads.
- Stack samples of the main thread in that phase (gdb, 30 samples): 28 are in
  `GHash<int>::Clear`, inlined in `main()` from `hashread.Clear()` at
  `stringtie.cpp:953`. That code is in `stringtie.o`, which is the same object
  file in both binaries.
- `main()` sits 112 bytes later in the branch binary than in master. Relinking the unchanged branch objects with padding in `.text.startup` before
  `stringtie.o` gives 29.2 s (no padding), 27.9 s (16 bytes), 25.2 s (32) and
  24.8 s (48). Master relinked with 8 bytes of padding in `.text` runs in 28.0 s
  instead of 24-25 s. So the speed of that one loop depends on where it lands,
  for master as well.

Two earlier explanations of mine were wrong and are recorded as such: a changed
inlining budget in `processRead()` (marking the index functions `noinline` changed
nothing), and the linker picking a different template instance of `Clear()` (it is
inlined; there is no such symbol).

With the cost of that loop removed (section 9), the index adds nothing measurable:

```text
Run, build of section 9 plus ...      Wall, three runs
------------------------------------  -----------------------
c21_S, master search, -p 1            14.00  14.00  13.99 s
c21_S, run index, -p 1                14.15  13.94  13.88 s
c21_mixG, master search               15.66  15.21  15.49 s
c21_mixG, run index                   15.61  16.31  16.15 s
c21_mixG, run index, threshold off    15.46  16.15  15.50 s
```

On `c21_mixG` the index is used for 401,524 of 17.5M reads; with the threshold set
so high that it is never used the times cover the same range, so the 0.5 s spread
there is not the index.

Consequence for review: on its own, this branch can show up as slower on deep
paired short-read data purely through layout. It should be timed together with
`ghash-clear`, or merged after it.

### Not covered

- The reporter's data. Per-bundle equivalence there should be checked with a
  `RUNIDX_CHECK` build, or by comparing GTFs at `-p 1`.
- `-U` (unitig) input: no local data; the flag is in the key and the predicate.
- Thread safety beyond the existing single-producer assumption.

## 7. Hash functions

The branch uses `wyhash` from `gclib/wyhash.h` over a buffer of 32-bit fields.
For comparison I built the same index with other hash functions, and one version
with the table replaced by `std::unordered_map<std::string,int>` (key bytes to
newest run position). Each ran on the 640k synthetic file, on
`issue_476/bundle_inf.bam` and on `mix/chr1_lr.bam`, with thresholds 64 and 2
(the latter sends nearly every lookup through the hash), default threads, best of
two runs.

```text
Variant                               syn640k        i476          chr1_lr
                                      thr64  thr2    thr64  thr2   thr64  thr2
------------------------------------  -----  -----   -----  -----  -----  -----
wyhash (branch)                        5.46   4.84    3.24   3.29  15.78  15.85
std::hash<std::string> (copies key)    4.87   4.78    3.32   3.25  15.76  15.79
libstdc++ std::_Hash_bytes, no copy    4.91   4.87    3.32   3.26  15.82  15.82
hash_combine over std::hash<uint32_t>  5.02   4.85    3.31   3.29  15.81  15.82
FNV-1a 64                              4.89   4.81    3.30   3.29  15.78  15.87
std::unordered_map<std::string,int>    5.18   5.51    3.25   3.25  15.83  15.84
```

Wall seconds. Table quality, threshold 64:

```text
Variant            Bucket steps per lookup      Equal 64-bit hash, different key
                   syn640k   i476   chr1_lr     syn640k   i476   chr1_lr
-----------------  -------   -----  -------     -------   ----   -------
wyhash              1.004    1.458   1.225          0       0        0
std::hash<string>   1.006    1.446   1.144          0       0        0
_Hash_bytes         1.006    1.446   1.144          0       0        0
hash_combine        1.013    1.406   1.177      6,050       0       68
FNV-1a              1.003    1.415   1.083          0       0        0
```

- Results: every variant produces a GTF identical to the reference on all three
  inputs (for the synthetic file the reference is the branch binary, see section
  6). `RUNIDX_CHECK` builds of the `hash_combine` and `unordered_map` variants
  report no mismatch on `bundle_inf.bam`, chr21 long `-G` and the `mixbug` mix.
- Speed: no difference I can measure. The spread between variants is the size of
  the spread between the two thresholds of one variant (wyhash 5.46 vs 4.84 on
  the synthetic file is run noise). Hashing 20-100 bytes per read is a small
  part of 5-10 us per read.
- Quality: `hash_combine` over libstdc++'s identity `std::hash<uint32_t>` gives
  6,050 pairs of different keys with the same 64-bit hash on the synthetic file
  and 68 on chr1. The others give none. The index tolerates this because each
  candidate is confirmed with the original comparison; a design that trusted
  the 64-bit hash alone would silently merge different reads with that function.
- `std::unordered_map`: same results, up to 14% slower on the synthetic file at
  threshold 2 (a string key is built for every lookup), and it brings STL
  containers into `rlink.cpp`. No reason to prefer it.
- The longest bucket chains (75-164 on chr1) are identical alignments that are
  never collapsed (section 8, defects 2 and 3): their stored keys are equal, so
  they share a bucket. A lookup passing through compares only the stored hashes.

wyhash stays: it is already in gclib, needs no key copy, and has no collisions
here.

## 8. Defects in the handling of deletions next to introns

Found while writing test cases for section 4. Details, reproducers and the proposed
fixes are in `worklogs/fable-junction-deletions-issue-draft.md`; in short:

```text
#  Where                         What happens                                  Visible effect
-  ----------------------------  --------------------------------------------  --------------------------
1  gclib/GSam.cpp:228,239        GSeg(prevdel,0) swaps its arguments: a         reads not rescued (longcov
                                 deletion before an intron is recorded as one   6 instead of 10), or exon
                                 after it                                       end off by the deletion
2  gclib/GSam.cpp:227-241        a deletion after an intron adds two juncsdel   identical reads never
                                 entries for one intron                         collapse; offsets of later
                                                                                junctions misapplied
3  rlink.cpp:239 (deljuncmatch)  with -G, offsets of a guide-matched junction   identical reads never
                                 compared with 0                                collapse
4  rlink.cpp:17617,17629         every read adds the whole shared junction      junction count 12 instead
                                 count to the rescued junction                  of 6; an error splice site
                                                                                reported as an isoform
```

The same code is in 3.0.3 and 2.2.3.

Reproducers: `issue_511/junction-deletions/` holds 8 small SAM files, `guide.gtf`,
the generator `mkreads.py`, and `check.sh STRINGTIE`, which prints ok/FAIL per case
and exits 1 on any failure. Master fails 7 of 9 cases; with fix 1, 4; with fixes
1-2, 2; with fixes 1-3, 1; with all four, 0. This script is the candidate
regression test.

Branch `lr-junction-deletions` (local): `9ccebd9` (1), `facc42f` (2), `c8e2046`
(3), `85b71ae` (4). `make test` on it: tests 6 and 7 (`--mix`) differ from the
expected GTFs in values only (largest cov change 0.23% and 3.5%); the rest pass.

Effect on real data, each fix added in turn, gffcompare against the annotation in
the same directory (IC = intron chain level):

```text
Run           Build      Distinct    Tx     Structure lines   IC sens   IC prec
                                            differing
------------  ---------  ----------  -----  ---------------   -------   -------
c1_L          master       456,123   4,161          0           20.3      40.5
              fix 1        456,471   4,153         54           20.3      40.6
              fixes 1-2    445,880   4,169         76           20.3      40.4
              fixes 1-3    445,880   4,169         76           20.3      40.4
              fixes 1-4    445,880   4,078        185           20.4      41.5
c1_LG         master       469,461   3,749          0           36.8      87.6
              fix 1        469,474   3,742         53           36.9      87.9
              fixes 1-2    469,345   3,743         62           36.8      87.8
              fixes 1-3    426,883   3,743         62           36.8      87.8
              fixes 1-4    426,883   3,743         62           36.8      87.8
c21_L         master        54,594     526          0           18.4      34.7
              fixes 1-2     53,960     533          7           18.4      34.3
              fixes 1-4     53,960     521         15           18.1      34.5
c21_LG        master        55,068     452          0           35.3      80.9
              fixes 1-4     52,713     451         11           35.1      80.6
mixbug_L      master        16,133     178          0            6.8      27.5
              fixes 1-4     16,094     162         24            6.8      30.3
mixbug_mixG   master       131,053     343          0           35.4      82.5
              fixes 1-4    130,445     345         12           35.4      81.9
c21_mix       master     1,043,071     812          0           22.9      30.7
              fixes 1-4  1,042,425     809         61           22.9      30.8
i476          master        33,280      20          0             -         -
              fixes 1-4     33,208      20          0             -         -
```

Reading: the fixes are corrections, not an accuracy gain. Changes are within about
1 point either way; fix 4 raises precision without guides (c1_L 40.4 to 41.5,
mixbug_L 28.4 to 30.3) and lowers sensitivity slightly on chr21. Fix 3 removes 9%
of the stored reads on c1_LG with no change in structures. Fix 4 changed nothing in
the `-G` runs (GTFs with fixes 1-3 and 1-4 are identical). The annotation is
RefSeq and the reads are NA12878 direct RNA aligned with minimap2 2.17, where 24%
of spliced alignments have such a deletion; on the `issue_465` long-read BAM it is
0.6%.

Whether the rescue is meant to work as I describe it is a question for Mihaela;
the issue draft asks it.

### Interaction with the run index

`CRunIndex::keyAln()` uses the raw `juncsdel` offsets of the incoming alignment,
matching today's `deljuncmatch()`. If fix 3 is merged, `keyAln()` has to use 0 for
guide-matched junctions as well, or reads that now match are missed by the index.
A `RUNIDX_CHECK` build shows this immediately.

## 9. The mate hash is scanned at every bundle (branch `ghash-clear`)

This came out of the timing in section 6 and is independent of issue #511, but it
is the same kind of problem: the single producer thread doing work that grows with
the largest bundle seen.

`hashread` (read name to `readlist` index, used to pair mates) is one
`GHash<int>` for the whole run, cleared at each new bundle
(`stringtie.cpp:953`). `GHash::Clear()` tests every bucket to free the key
strings, and the table never shrinks ("does not shrink!" in the source). Once a
very deep paired-read bundle has grown it, every later bundle pays for a scan of
the table at that size. On `mix/chr21_short.bam` one bundle has 14.6M of the 17.4M
reads; the roughly 3,000 bundles after it take 14 s, 4 to 5 ms each. A
whole-genome sample with a deep locus early on (rRNA, mitochondrial) would pay the
same on every later bundle. Paired short reads only: long reads do not use this
table.

Changes in `gclib/GHashMap.hh` (the request was to optimize GHash, `Clear()`
included):

- `Clear()` returns at once for a table without items (the used flags are already
  all zero: deletion clears them, there are no tombstones).
- The scan reads the used flags a 32-bucket word at a time and stops after the
  last item. The same helper (`gh_next_used()`) replaces the bucket-by-bucket
  loops in `Next()`/`NextData()`.
- A table above 65,536 buckets with fewer than 1/16 of them in use is freed by
  `Clear()`; the next insertion starts a small table. A table that is refilled to
  a similar size each time is kept.
- Bug fix: `GHashSet::Next()` returned `key(i_iter-1)`, the bucket before the one
  found. Nothing in stringtie calls it.

Measured with `issue_511/ghash/ghash_test bench` (old header, new header):

```text
Operation                                              Old          New
-----------------------------------------------------  -----------  ----------
Clear(), 8.4M buckets, 100 keys left (first time)        3,709 us    13,226 us
then 200 Add, 190 Remove, Clear, per cycle               3,701 us         7 us
Clear() of the emptied table                             3,667 us     0.001 us
iteration, 1,000 items in 4.2M buckets                   2,017 us        79 us
GIntHash: 1M Add then Clear, per cycle (Add / Clear)   12.9 / 0.002 ms   12.7 / 0.002 ms
```

The first `Clear()` is slower once because it returns about 200 MB to the system.

StringTie with the new header, default threads:

```text
Run                        Master     ghash-clear   Output
-------------------------  ---------  -----------   ---------
c21_S (short)              24.5 s       13.6 s      identical
c21_S -G -A                    -        13.9 s      identical
c21_mixG (--mix -G)        21.6 s       15.7 s      identical
c21_LG, mixbug_mixG, i476      -           -        identical
```

Checks: `ghash_test` runs random Add/Remove/Find/Clear/iteration sequences on
`GHash<int>`, `GHash<Obj*>` (owned values), `GHashMap<int,int>`, `GIntHash<int>`,
`GHashMap<int,Obj*>`, `GStrSet` and `GHashSet<int>` against `std::map`/`std::set`,
through phases that grow, empty and regrow the tables. Old and new headers give
the same checksums and pass every check, except that the old `GHashSet::Next()`
returns 85,435 wrong keys. Under valgrind the new build has 0 errors and no
leaks. `make test`: 12 OK.

One behavior to be aware of: the order in which `Next()` returns items depends on
the table size, so after a `Clear()` that freed a large table, items inserted
later come back in the order of a freshly grown table. The three places that
iterate (`stringtie.cpp:2134`, `rlink.cpp:16510`, `GFastaIndex.cpp:162`) do not
reuse a cleared large table as far as I can see, and all outputs above are
identical.

Not done: the hashing itself is under 2% of a short-read run in the callgrind
profile (`Add` 0.6%, `strlen` 0.7%), so I left the hash functions and the probe
code alone. `khashl.hh` `resize()` calls `memset` on a new array before checking
that the allocation succeeded; noted, not changed. `gclib/GHashMap.hh` here is
identical to the copy in `~/src/gclib`, so the change applies there as is.

An earlier version of this branch (commit `1e59117`, now replaced) worked around
the problem in `stringtie.cpp` by freeing `hashread` when it was large; it gave
the same times. The fix in `Clear()` makes that unnecessary.

## 10. Decisions and open items

Decisions from the maintainer (2026-11-12):

1. Success criterion: byte-identical GTF at `-p 1` against master on real long,
   short and mixed data, the per-read check, and no measurable slowdown on normal
   data. A small numeric drift would also be acceptable if accuracy is unaffected
   and the gain is large.
2. The fix goes to a separate branch for discussion with Ida, not a PR.
3. The collapse defects get their own issue and branch, with reproducers that can
   become regression tests.
4. Reporter data: can be requested; local data first.

Open:

- None of the three branches (`i511-long-reads`, `lr-junction-deletions`,
  `ghash-clear`) is pushed; `git push -u origin <branch>` when wanted. Each is
  based on master and independent of the others.
- `i511-long-reads` should be timed together with `ghash-clear` (section 6).
- The junction-deletion issue is drafted, not filed. I held it because the scope
  grew from the two collapse cases that were approved to four defects, two of which
  change assembly results, and filing notifies the repository watchers.
- Ida's patch has not been seen. The same battery and a `RUNIDX_CHECK`-style check
  apply to it; section 4 lists what to look for.
- A bundle-wide hash (drift accepted) would be simpler than the run index, but it
  only differs on the few reads of section 4.2 and saves no time over it, so there
  is little reason to prefer it.
- `rlink.cpp:17788` repeats the count pattern of defect 4 for unstranded junctions;
  not changed.
- `gclib/GSam.cpp` is vendored; fixes 1 and 2 belong in gclib as well.

## 11. Files

```text
Path                                         Content
-------------------------------------------  ---------------------------------------------
issue_511/synthetic/gen.py                   synthetic same-start locus (section 3)
issue_511/junction-deletions/                reproducers, mkreads.py, guide.gtf, check.sh
issue_511/tx.sh                              GTF -> sorted "ref strand exons" per transcript
issue_511/run-index/                         battery, timing, hash comparison and sampler
                                             scripts with their summary tables (copies)
issue_511/junction-deletions/eval.sh         effect of fixes 1-4; eval_summary.tsv
issue_511/work/i511/                         where those scripts ran: binaries, logs, GTFs
issue_511/work/i511/st_head, st_br           master and branch release binaries
issue_511/work/i511/st_brchk, st_brchk2      RUNIDX_CHECK builds, thresholds 64 and 2
issue_511/work/i511/st_time                  master with scan timer and counters
issue_511/work/i511/src/                     export of master used for instrumented builds
issue_511/work/i511/pmp.sh                   stack sampler (gdb as parent; perf and attach are blocked)
issue_511/work/i511/st_gh                    release binary of branch ghash-clear
issue_511/ghash/ghash_test.cpp               GHash test and benchmark (section 9)
issue_511/work/collapse/st_d0 .. st_d4       master and cumulative fix builds
issue_511/work/stringtie_i511, _lrdel        release binaries of the two branch tips
```

`issue_511/` is untracked, like the other `issue_*` directories. `work/` is about
3 GB, mostly SAM files, binaries and GTF outputs that the scripts regenerate; it
can be deleted once the numbers are no longer needed.

To rebuild a check binary from the branch:

```sh
git checkout i511-long-reads
make clean && CXXFLAGS="-g -O3 -DRUNIDX_CHECK" make -j8 release
```

## 12. Problems during this work

- I deleted the baseline binary copy with `rm ../st_h*` while cleaning hash
  variants. The dataset battery was running; its rows from `c21_S` on were invalid.
  I restored the binary from the unchanged object files, verified that it
  reproduces two earlier baseline GTFs, and reran those datasets.
- The first battery did not record run times (`/usr/bin/time` output redirected to
  `/dev/null`); timing was redone separately (section 6).
- The hash comparison began with master on the 640k synthetic file (the quadratic
  case, 29 minutes) without my having said so. I then stopped that run after 4.6
  minutes when asked about the delay. Stopping it was the mistake: the result was
  needed, and I had to run it again from the start.
- Twice I killed my own shell with `pkill -f PATTERN` because the pattern also
  matched the command line that contained it.
- The short-read slowdown (section 6) took two wrong explanations before stack
  samples and relinking with padding settled it.
- A binary that calls `GError()` from the producer thread (my negative control) did
  not exit: it stayed asleep after printing the message. I did not investigate; it
  may affect real error exits in threaded runs.
- An early table of the reproducers reported `longcov` as `cov`. The regression
  script exposed it; the numbers in section 8 and in the issue draft are from the
  corrected parsing.
- The session scratch directory was retired midway. Experiment files were moved to
  `issue_511/work/` and the temporary worktrees removed; both branches are intact.
- `.gitignore` and `AGENTS.md` in the working tree were modified at 10:01 by another
  agent (codebase-memory-mcp setup text). I did not touch them, and committed in
  separate worktrees so they are not part of either branch.
