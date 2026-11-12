# merge-dev review: what to fix before merging to master

Reviewed: branch `merge-dev` at `a27c05d` against `master` at `d1dc38d`
(StringTie 3.0.3). Review date 2026-11-09. Sections 1 to 9 are the plan plus the
evidence and reasoning behind it. Sections 10 to 12 (revision 4) record what
has since been implemented in the working tree, uncommitted, the results on a
second dataset, and where the implementation departs from the plan. Where
sections 4 and 8 disagree with section 10, section 10 is what the code does.

Revision 2 (2026-11-09). Incorporates the independent assessment in
[review-branch-response1.md](review-branch-response1.md), whose additional
failures all reproduce, and the maintainer decisions recorded in section 9.
Changes against revision 1: F1 to F4 are now implemented through one
restructuring of `merge_transfrags()` (section 4.1) under one option model
(section 4.2); F6 is broader and decided as pass-through; the `auto` default
proposal is withdrawn; two acceptance criteria were replaced; sections 3.5 to
3.8 are new evidence.

Revision 3 (2026-11-09). Every policy choice that was still open becomes a
command-line option with a stated default (section 4.3); a second dataset is
being prepared (section 5).

## 1. Verdict

Not ready to merge. The branch delivers the headline fix for issue #247 (a `-G`
transcript is never rewritten), but its own stated contracts do not hold yet:

```text
Contract stated in README/worklog                 Observed
------------------------------------------------  ---------------------------------
New 5p/3p ends need --min-support-5p/3p files     Not enforced for terminal-exon
or must fit the slop; hard caps always win        extension, the #247 mechanism
--guide-policy and caps govern guide variants     Bypassed by partial chains and by
                                                  variants shifted at both ends
One vote per file decides what is kept            Absorption ignores support: 17,800
                                                  replicated chains are lost under a
                                                  weaker path (section 3.5)
--path-mode restricts joins                       Applied after inputs are consumed;
                                                  direct mode can emit nothing
Every filtering decision is in --report           Absorbed inputs have no report row
```

In addition the options use two different definitions of "same structure"
(section 4.2), three behaviors regressed relative to 3.0.3 (F5, F6), and the
test script's negative assertions cannot fail (F7).

## 2. What already holds (verified, keep it)

- `make test-merge` passes; the 9 assembly tests in `run_tests.sh` pass; the
  merge fixtures are valgrind-clean.
- Issue #247 dataset (20 GTFs plus guide), new defaults: 192,173 transcripts,
  26,059 of 26,059 guide IDs keep their exact exon structure (3.0.3 changed
  13,958), 0 duplicate transcript IDs. Matches the worklog.
- Thread determinism: the output structure set is identical at `-p 1`, `4`, `8`
  and the report is byte-identical. GTF text differs only in `MSTRG` numbering
  and record order, as in 3.0.3.
- Printing and the generated-structure registry run under `printMutex`; report
  appends are under `mergeReportMutex`. No data race found by reading.

These were reproduced independently by the second reviewer. They are useful
results, not a broad correctness proof: the suite has no direct-mode test, no
support-value assertions, and most checks inspect outer bounds only.

## 3. Evidence base

All numbers below come from the issue #247 dataset in
`merge-dev/merge-dev/example_dataset/` (20 sample GTFs, 927,737 input
transcripts, guide with 26,059 transcripts), from the release build of
`a27c05d`, and from the 3.0.3 binary built from `master`. Small fixtures are
listed in section 8 so that every claim can be re-run.

Limits that apply throughout: one dataset, a brown alga with an incomplete
guide, so "novel" is not "false". The samples were assembled with `-G`, so
reference transcripts carry the guide's exact ends and their support is
inflated.

### 3.1 Support distribution in the inputs (independent of the merge code)

MET intron chains, counted by distinct input files carrying the chain:

```text
Files  Reference chains  Novel chains  Novel with max TPM >= 1
-----  ----------------  ------------  -----------------------
1                   473       124,086                   83,283
2                   370        25,857                   21,369
3                   329        11,324                   10,157
4                   330         6,295                    5,952
5+               14,049        22,160                   21,909
total            15,551       189,722
```

SET records, counted by distinct files holding a SET with reciprocal overlap of
0.80 or more:

```text
Files  SET records  With TPM >= 1  Mean length
-----  -----------  -------------  -----------
1          106,308         50,298          385
2           19,244         11,451          864
3           11,197          7,121        1,028
4            8,585          5,901        1,097
5+          86,358         58,254        1,030
```

### 3.2 Retention of input intron chains in the merged output

```text
Chain support   New, no guide   New, with guide   3.0.3, no guide
-------------   -------------   ---------------   ---------------
1 file                  47.7%             41.3%             11.6%
2-4 files               38.1%             35.9%             17.3%
5-9 files               45.6%             46.2%             32.8%
10+ files               57.7%             75.6%             55.0%
```

Retention need not be monotone in support (a lost chain can be a redundant
fragment), so this table is descriptive. Section 3.5 separates the two cases.

### 3.3 Runs at different settings (no guide unless noted)

```text
Run                          Transcripts     SET     MET  Ref chains  Novel chains
---------------------------  -----------  ------  ------  ----------  ------------
3.0.3 defaults                    44,428   3,908  40,520       8,318        32,202
3.0.3 -F 0 -T 0 -f 0             105,192  19,397  85,795       9,132        76,663
new default                      187,146  88,370  98,776      10,021        85,513
new --min-support 2               20,746   3,877  16,869      13,691         3,178
new -T 1                         118,608  43,518  75,090      10,019        62,564
new default, with -G             192,173  90,054 102,119      22,258             -
```

"Ref chains" is the number of reference intron chains present in the output,
out of 22,258 in the guide (15,551 of which occur in the inputs).

### 3.4 Composition of the new default output (with guide)

Of 166,114 generated (non-guide) transcripts, 128,206 (77%) have
`merge_support "1"`: 81,249 of 87,544 generated SETs and 46,957 of 78,570
generated METs. The growth from 63,959 (3.0.3) to 192,173 records is almost
exactly this single-file class (+128,214).

### 3.5 Why input chains are missing from the output (no guide, new default)

For each input chain absent from the output, its file support is compared with
`merge_piece_support` (weakest-intron support) of the best output path whose
chain contains it.

```text
Chain     Lost    Under a path with   Only under a   No containing
support           piece >= support    weaker path    output path
--------  ------  ------------------  -------------  -------------
1 file    65,186              63,106              0          2,080
2-4       27,545              20,494          6,448            603
5-9        7,898               3,787          3,942            169
10+        9,190               1,589          7,410            191
```

Most lost chains are fragments of an equally or better supported path, which
is legitimate. But 17,800 replicated chains, including 81% of the lost chains
seen in 10 or more files, sit only under a weaker path. 3,043 chains (963
replicated) have no containing output path; that loss is not traced.

### 3.6 Exact-structure support versus intron-chain support

Chains seen in 2 or more files that also have one exact structure (both ends
identical to the base) in 2 or more files:

```text
Class   Chain support  Chains  Exact structure in >= 2 files
------  -------------  ------  -----------------------------
ref     2+             15,078  15,078 (100%)
novel   2-4            43,476     574 (1.3%)
novel   5-9            12,631   1,036 (8.2%)
novel   10+             9,529   2,570 (27.0%)
novel   2+             65,636   4,180 (6.4%)
```

Reference chains reach 100% only because guided assembly copies the guide
ends. With the exact-structure key, `--min-support 2` therefore means "keep
reference transcripts, drop 94% of replicated novel chains" on this dataset.

### 3.7 End variability (novel chains seen in 3 or more files; 39,779 chains)

Cumulative share within a distance, in percent:

```text
Measure                              End   0   25   50  100  200  500  1000
-----------------------------------  ---  --  ---  ---  ---  ---  ---  ----
One file's end vs the chain median   5p   25   73   82   88   92   97    99
                                     3p   22   55   65   74   82   91    96
Most extended file beyond the        5p   14   64   74   82   89   95    98
second most extended file            3p    9   47   59   70   81   91    96
```

Partial-chain inputs anchored at a guide's terminal intron (input records, not
chains): 36.5% of 49,523 extend past the guide 5-prime end and 38.5% of 66,575
extend past the guide 3-prime end. Of the extending records:

```text
End   Records  <=25  <=50  <=100  <=200  <=500  <=1000
---   -------  ----  ----  -----  -----  -----  ------
5p     18,076   35%   47%    64%    80%    97%     99%
3p     25,631    6%   10%    16%    24%    44%     71%
```

Reading: 3-prime ends vary 2 to 3 times more than 5-prime ends. 100 bp covers
most 5-prime variation. The 3-prime distribution has no natural break, and 29%
of 3-prime extensions past a guide exceed 1 kb, so any 3-prime tolerance is a
policy choice, not a measured optimum. Full-chain matches to guides are
uninformative here (guided assembly copied the guide ends).

Single-file overhang (stretch between the second most extended and the most
extended file) for the 65,636 novel chains seen in 2 or more files, and how
often it reaches an annotated transcript on the same strand lying wholly
beyond the 2-file end:

```text
End  Overhang    Chains  Reaches another guide transcript
---  ---------   ------  --------------------------------
5p   1-100       44,215     31 (0.1%)
5p   101-500     10,299     64 (0.6%)
5p   501-1000     2,115     59 (2.8%)
5p   >1000        1,456     90 (6.2%)
3p   1-100       38,424     50 (0.1%)
3p   101-500     14,932    103 (0.7%)
3p   501-1000     3,614     95 (2.6%)
3p   >1000        3,221    183 (5.7%)
```

An unbounded 3-prime slop would let one file add more than 500 bp to 6,835
chains (10.4%) and reach a neighboring annotated transcript in 431, against
153 at a slop of 500 and 50 at 100. Collisions are rare in this genome; the
main effect is length. Total-RNA libraries with readthrough are the case to
check before choosing (section 5).

### 3.8 Effect of the number of input files (mean of 5 random subsets)

```text
N files  Ref >= 1  Ref >= 2  Novel >= 1  Novel >= 2  Novel >= 2 / >= 1
-------  --------  --------  ----------  ----------  -----------------
2          12,352     8,855      37,743       6,902               0.18
3          14,027    11,653      52,378      11,518               0.22
4          14,322    12,551      62,508      15,541               0.25
6          14,678    13,502      85,742      24,029               0.28
10         15,149    14,407     119,338      37,000               0.31
20         15,551    15,078     189,722      65,636               0.35
```

Reference chains saturate by about 4 files. Novel chains at support 1 grow
almost linearly (about 8,000 new chains per added file) with no sign of
saturation, which is how noise or sample-private models behave. At support 2
the count grows at about a third of that rate. At N = 2 a threshold of 2 would
remove 28% of the reference chains present.

## 4. Fixes

Sections 4.1 and 4.2 are the design agreed in revision 2. F1 to F4 describe
the individual defects and are fixed inside that design, not as separate
patches. Order of work is in section 6.

### 4.1 Restructure `merge_transfrags()`: observed chains first

Why. Every policy filter currently runs after destructive coassembly. A seed
consumes contained and joinable transfrags by zeroing their abundance
(`rlink.cpp:7789` to `7799`), one path per seed is stored, and the cleanup
loop (`rlink.cpp:7866` to `7915`) removes subset predictions and stretches the
survivor's ends. If the path is later rejected, or is weaker than what it
consumed, the observed inputs are already gone. F1, F3, F4 and the direct-mode
loss are instances of this one defect. Patching each site leaves the next one
waiting.

Design.

1. Candidates. Every guide is stored as today. Every non-guide intron chain
   that passes structure support (section 4.2, layer 1) becomes a candidate
   prediction with its own ends elected per F1. Transfrags are pooled by node
   pattern, and a terminal exon can span several nodes, so several transfrags
   can carry one chain: group them by chain before storing.
2. Synthesis. The existing seed-and-extend loop produces joined paths, but as
   additions: it no longer consumes candidates. A joined path is stored only
   if its chain is not already a candidate and it passes `--path-mode`,
   `--min-introns`, the endpoint policy and the gene/CDS barriers at build
   time. `--path-mode direct` skips synthesis altogether.
3. Redundancy. One explicit final pass replaces the cleanup loop. Candidate A
   is removed in favor of P only if all hold:
   - chain(A) is a contiguous sub-chain of chain(P);
   - A's outer ends fit P: not beyond P's ends by more than the slop of that
     end (section 4.2, layer 2);
   - an end of A that overhangs into an intron of P does so by at most
     `--max-intron-overhang` (default unlimited; 0 keeps every such fragment
     as a separate structure). An overhang spanning the whole intron is a
     retained intron and stays under the existing `-i` rule;
   - dominance: `merge_piece_support(P) >= support(A)`. A guide P is exempt
     unless `--guide-fragments dominance` is given.
   The survivor's ends are not stretched. Files of A are credited to P.
4. Accounting. Every chain removed in step 3 or rejected by any layer gets a
   report row with the reason, the surviving transcript if any, the number of
   input records and the number of files. Only then can report totals be
   reconciled with the 927,737 inputs.

Measured cost. Dominance releases about 17,800 chains on this dataset
(section 3.5), weighted toward the best-replicated ones. Accepted by the
maintainer. Runtime and memory effect is unknown until implemented; measure
on the full dataset.

Invariant to test. Every input chain that passes layer 1 is in the output, or
has a report row naming a transcript that contains it and satisfies
dominance. This replaces "monotone retention" and "at least 13,691 reference
chains" from revision 1.

### 4.2 One option model: layers that only restrict

Problem. The branch has three notions of "the same structure", and options
keyed on different ones contradict each other:

```text
Option or attribute                 Compares                       Site
----------------------------------  -----------------------------  ----------------
--min-support, --rescue-support     every exon start and end,      merge_read_key()
                                    including outer ends           rlink.cpp:16207
merge_path_support, --path-mode     intron chain only              merge_same_chain()
                                                                   rlink.cpp:7081
--min-support-5p/3p, slop, caps     graph-node bounds; skipped     merge_can_join()
                                    when no node is added          rlink.cpp:7592
--guide-policy, caps for variants   one span containing the other  merge_related_guide()
                                                                   rlink.cpp:20068
```

Consequences: fixture S (one chain in 3 files, ends differing by 10 bp) is
rejected entirely at `--min-support 2` while its `merge_path_support` would be
3; only 6.4% of replicated novel chains have an exact structure in 2 files
(section 3.6); the end options never see the ends they were designed for.

Model. Identity is defined once and each layer can only remove or trim what
the previous layer passed. No layer re-admits anything; `--rescue-support` is
the single explicit exception and uses the layer 1 identity.

```text
Layer            Question                         Options
---------------  -------------------------------  ---------------------------------
0 record         Is this input record usable?     -m -c -F -T -f
1 structure      Is the structure replicated?     --min-support, --rescue-support,
                                                  --set-overlap
2 ends           Where does each end go?          --min-support-5p/3p, --slop-bp-*,
                                                  --slop-frac-*, --max-extension-*
3 guide          Is it the guide, a variant, or   --guide-policy, --guide-replace-frac,
                 novel?                           layer 2 slop and caps
4 synthesis      May observed chains be joined?   --path-mode, --min-introns,
                                                  --min-piece-support, --allow-*
5 redundancy     Is it a fragment of a kept one?  dominance (no option)
```

- Identity. MET: strand plus ordered intron chain. SET: strand plus pairwise
  reciprocal overlap >= `--set-overlap` (no transitive closure). Structure
  support is the number of distinct files holding that identity.
- Layer 2 is the only place ends are judged, for every transcript, whether or
  not a graph node is added (F1).
- Layer 3 uses the same intron-chain identity and the same per-end slop and
  caps as layer 2, measured from the guide end (F3).

This is a contract change for `--min-support` (README today: "files required
for an exact structure"). An exact-ends mode is not kept: it cannot be made
consistent with layer 2, and no dataset result argues for it.

### 4.3 Policy switches (no choice hard-wired)

Maintainer requirement: each remaining policy choice is selectable on the
command line. Names are proposals. Every switch sits in one layer and, like
the existing options, can only restrict what earlier layers passed.

```text
Option                        Layer  Default      Effect
----------------------------  -----  -----------  ------------------------------
--end-variants MODE           2      single       single: one model per chain
                                                  with elected ends. replicated:
                                                  also emit an end variant when
                                                  a separate end cluster has
                                                  --min-support-5p/3p files
--slop-bp-5p N                2      100          unchanged
--slop-bp-3p N                2      500 (prov.)  was 100; accepts "unlimited"
--guide-slop-bp-5p/3p N       3      = slop-bp    fuzz allowed past a guide end
                                                  before an input is a variant
--guide-fragments MODE        3,5    absorb       absorb: a guide absorbs every
                                                  contained fragment. dominance:
                                                  guides follow the dominance
                                                  rule like any other transcript
--max-intron-overhang N       5      unlimited    largest overhang into an
                                                  intron of the containing
                                                  transcript for a fragment to
                                                  be absorbed; 0 = always keep
                                                  such fragments separate
```

Notes.

- `--end-variants replicated` needs end clusters (ends within the slop of each
  other), not reach, because every file reaching a long end also reaches a
  shorter one. It is the one switch with real design risk; implement it last
  and keep `single` as default until it has fixtures and a dataset measurement.
- `--slop-bp-3p unlimited` makes `--min-support-3p` inactive unless
  `--slop-frac-3p` or `--max-extension-3p` binds. That is a legitimate setting
  but as a default it would switch the 3-prime end layer off. 500 is
  provisional: it keeps the most extended 3-prime end for 91% of replicated
  novel chains (section 3.7). Final value after the survey of section 5.
- `--max-intron-overhang unlimited` follows the maintainer's reading that such
  overhangs are usually incompletely processed RNA. Dominance still applies,
  so a fragment better supported than its container is kept either way.
- All defaults in this table are to be confirmed on the second dataset, which
  has matched polyA and total-RNA libraries and so tests exactly the intron
  overhang and 3-prime readthrough cases.

### F1. Apply the endpoint policy to terminal-exon extensions

Problem. `--min-support-5p/3p`, `--slop-*` and `--max-extension-*` only act
when a join adds new graph nodes. A longer first or last exon adds no node, so
it is never evaluated. This is the exact mechanism of issue #247 ("extending a
terminal exon introduces no new junction", `worklogs/merge-issues.md`).

Evidence.

- Fixture G: chain 1000-4300 in 2 files plus one partial-chain file ending at
  5300 gives a single output 1000-5300 with `merge_3p_support "1"`; the 2-file
  end is gone.
- Fixture H/H2: same with a full-chain single-file extension; identical result
  and `--max-extension-3p 100` does not stop the 1,000 bp extension, although
  the README says "Hard caps always win".
- Fixture K: one file extended on the left, one on the right, gives 300-5300,
  a pair of ends never observed together.
- Fixture I versus H: a left extension yields two outputs, a right extension
  yields one. The left/right asymmetry of 3.0.3 is still present, so on the
  minus strand it is the 5-prime end that is silently extended.
- Dataset: of 27,208 generated METs whose exact intron chain occurs in 2 or
  more files, 12,851 (47%) have an end more than 100 bp beyond the point that
  2 files carrying that chain reach (5-prime 4,824; 3-prime 9,509). This is an
  upper bound because other chains sharing the terminal exon may support the
  end.

Root cause.

1. `merge_can_join()` (`rlink.cpp:7592`) measures extension between graph-node
   bounds. Terminal-exon variants share one terminal node, so the extension is
   computed as 0.
2. Inputs with the same node pattern are pooled into one transfrag whose
   `longstart`/`longend` hold the minimum start and maximum end
   (`rlink.cpp:4880` area). `merge_set_observed_ends()` (`rlink.cpp:7094`) then
   restores that union over all path contributors. The cleanup loop stretches
   ends once more (`rlink.cpp:7895` to `7904`).
3. The support passed to `merge_extension_allowed()` (`rlink.cpp:7571`) is the
   candidate transfrag's total file count, not the number of files observing
   that particular end.

Plan.

1. Stop taking the union. In `store_merge_prediction()` elect each end
   independently from the contributing non-guide reads whose terminal exon
   ends at the path's terminal splice site. `mt->read` already links every
   transfrag to its reads, so `rd.start`, `rd.end` and `tinfo->fileidx` are
   available without new storage.
2. Election rule for one end, with k = `--min-support-5p` or `-3p` for the
   biological end in question:
   - support(x) = number of distinct files whose end reaches x or beyond;
   - elected end = the most extended x with support(x) >= k, that is the k-th
     most extended end over distinct files;
   - if fewer than k files contribute, use k = number of contributing files
     (a single-file structure keeps its own end; it is governed by
     `--min-support`, not by the endpoint policy);
   - the emitted end may exceed the elected end up to the furthest observed
     end only if that distance is within `--slop-bp-*` and within
     `--slop-frac-*` of the elected terminal exon length;
   - `--max-extension-*` bounds the distance between the emitted end and the
     elected end (and the guide end for guide-related variants), always.
3. In `merge_can_join()` keep the new-node check, but compute distances from
   observed ends and pass the number of files observing the candidate's end.
4. Redefine `merge_5p_support` and `merge_3p_support` as "files reaching the
   emitted end within slop" and document it. Today they count exact coordinate
   matches, which is 1 for almost every novel transcript.
5. Write a report row (`trim`, `unsupported_end_extension`, file, distance)
   for each input whose end was not adopted.

Reasoning.

- The per-end knobs were designed for exactly this case (`planning.txt`:
  SET or UTR overhang, polymerase run-off, "biology is fuzzy" hence slop). As
  implemented they govern the rare case (new exon) and skip the common one.
- "Reach" semantics rather than exact coordinates: only 14% (5-prime) and 9%
  (3-prime) of replicated novel chains have their two most extended files end
  at the same base (section 3.7).
- The k-th most extended end is symmetric in strand and in left/right,
  independent of file order, deterministic, monotone in k, and k = 1
  reproduces the legacy union, so there is a clean opt-out.
- Alternatives rejected: the modal end is ill-defined when ends do not
  coincide; the median ignores the user's threshold; clustering ends into
  windows adds parameters without a clear default.
- One model per intron chain is the default; `--end-variants replicated`
  (section 4.3) emits further end variants.

Expected effect of the default slops (100 bp 5-prime, 500 bp 3-prime): the
most extended end is kept for 82% (5-prime) and 91% (3-prime) of replicated
novel chains and trimmed to the 2-file end for the rest (section 3.7).

Acceptance. Fixtures G, H, H2, K, I and their minus-strand mirrors.
Self-consistency on the dataset: every emitted end either reaches its
`--min-support-5p/3p` threshold by the definition in step 2, or lies within
the recorded slop of an end that does; 0 exceptions. The same-chain statistic
(12,851 of 27,208) is reported before and after but has no target, because a
correct fix may use evidence from other compatible chains.

### F2. Count `--min-support` by intron chain (MET) and by overlap (SET)

Problem. Support is keyed on the complete exon list including the outer ends,
so two files agree only if their transcript ends are identical to the base.
The README documents this ("exact structure"), so it is not a bug against the
stated contract; it is a contract that contradicts the rest of the option set
(section 4.2).

Evidence.

- Fixture S: one intron chain in 3 files with ends differing by 10 bp;
  `--min-support 2` rejects all three (`minimum_support`, support 1) and the
  output is empty.
- Dataset: 4,180 of 65,636 replicated novel chains (6.4%) have an exact
  structure in 2 or more files (section 3.6); the `--min-support 2` run keeps
  3,178 novel chains.

Root cause. `merge_read_key()` (`rlink.cpp:16207`) serializes strand plus every
exon start and end; `filter_merge_inputs()` (`rlink.cpp:16264`) applies
`--min-support` and `--rescue-support` to those groups.

Plan.

1. MET key = strand plus the ordered intron chain (drop the first exon start
   and the last exon end).
2. SET support = distinct files holding a SET on the same strand with
   reciprocal overlap >= `--set-overlap` with this one. Pairwise, no transitive
   closure, so a chain of partially overlapping SETs cannot vote for each
   other. This is the same criterion the output SET deduplication already uses.
3. Use the same grouping for `--rescue-support` and for the passing-peer
   logic, so that abundance rescue and structural support agree.
4. `merge_support` on an output transcript is the layer 1 count, so the
   attribute and the filter cannot disagree.

Acceptance. Fixture S keeps one model at `--min-support 2` and `3`, none at
`4`; the dataset run at `--min-support 2` keeps a number of novel chains of the
same order as the 65,636 eligible ones; README and `-h` text state the new
definition.

### F3. Guide relation: fragments, variants, and their accounting

Problem. The relation between an input and a guide is decided by span
containment, so partial chains and shifted variants escape the guide options,
and absorbed evidence is discarded without a trace.

Evidence.

- Fixture A: guide with 4 exons ending at 4300; one input with the last two
  introns and a last exon ending at 5300. Output is the guide only, with
  `merge_support "0"`, and the report is empty.
- Fixture D: same with `--min-support-3p 1`. Same result.
- Fixture E: the input present in 3 files. Same result.
- Fixture V (second reviewer): guide `1000-1199,1300-1499`, input
  `900-1199,1300-1399`. Same intron chain, longer at one end and shorter at
  the other. It is emitted as `MSTRG.1.1` with no reference relationship under
  the default, under `--guide-policy reference`, and under
  `--max-extension-5p 0`.
- Dataset: the two source models of the #247 example,
  `SRR12897873.7658.2` (cov 17.6) and `.3` (cov 95.2, TPM 11.9), both ending at
  4601396, are in neither the output nor the report. 7,478 guides show
  `merge_support "0"` even where many samples carry fragments.

Root cause. With a guide as path seed, the loop at `rlink.cpp:7789` adds any
transfrag whose pattern is contained in the path pattern (`rlink.cpp:7792`);
the extension is measured as 0 (F1), the seed's node range did not change, so
nothing is stored (`rlink.cpp:7811`). `merge_related_guide()`
(`rlink.cpp:20068`) requires one complete span to contain the other.

Plan (layer 3 of section 4.2).

1. Related to guide G: the input chain equals G's chain, or is a contiguous
   sub-chain of it.
2. Each end is judged separately, and only where the input shares G's terminal
   intron on that side (an end falling inside G is not an end variant):
   - not beyond G's end, or beyond it by at most `--guide-slop-bp-*` of that
     end (default: the value of `--slop-bp-*`): the input is G at that end
     (fuzzy TSS/TES);
   - beyond by more than the slop: an extension, which needs
     `--min-support-5p/3p` files reaching it and is bounded by
     `--max-extension-*` measured from G's end.
3. If no end is an accepted extension the input is a fragment of G: absorbed,
   credited to G (`merge_contained_support`, files), reported. With
   `--guide-fragments dominance` it is absorbed only if G's piece support is
   at least its own.
4. If an end is an accepted extension the input is a guide variant: emitted
   under a generated ID with `reference_id` under `--guide-policy both`,
   dropped with a report row under `reference`, handled as today under `best`.
   An extension that fails support or the cap is trimmed to G's end, which
   makes the input a fragment (item 3) with a `trim` row.

Tolerance default. `--guide-slop-bp-5p/3p` default to `--slop-bp-5p/3p`, so
there is one definition of fuzziness per end unless the user separates them.
On this dataset 100 bp treats 64% of 5-prime past-guide extensions as fuzz,
and 500 bp treats 44% of 3-prime ones as fuzz (16% at 100 bp; section 3.7).
The 3-prime extensions are long and smoothly distributed, plausibly
incomplete UTR annotation in this genome; check on the second dataset.

Reasoning. Keeping the guide unchanged is right and is what #247 asked for.
Discarding the evidence is a separate decision that the planning notes did
not want hard-wired: reference annotations "are never perfect", and an
extension seen in several samples should be something the user can keep or
reject.

Acceptance. A stays guide-only by default with a `trim` row and credited
support; D and E additionally emit the extended variant with a `reference_id`;
`--guide-policy reference` suppresses it. V at the default is absorbed into
the guide (100 bp is within slop); with `--guide-slop-bp-5p 50` it is a
variant under `both`, dropped under `reference`, and dropped by
`--max-extension-5p 0`, each with a report row. The two #247 source models appear in the report.

### F4. Absorption must respect support, credit it, and report it

Problem. A contained model is absorbed by a longer one regardless of how many
files support each, the absorbed files are not credited, and nothing is
reported. With `--path-mode direct` a rejected join also destroys the observed
inputs it consumed.

Evidence.

- Fixture Q (issue #370 shape): a 3-exon model in 1 file and a 2-exon isoform
  with a different 5-prime start in 3 files. Output is the 3-exon model only,
  with `merge_support "1"`, empty report. 3.0.3 does the same.
- Fixture W (second reviewer): inputs `1000-1100,2000-2100,3000-3300` and
  `2000-2100,3000-3300,4000-4300` with `--path-mode direct` and endpoint
  thresholds 1. The only result is a report row
  `1000 4300 reject direct_path_support 0`; neither original is emitted.
- Dataset: section 3.5. 17,800 replicated chains are lost only under a weaker
  path. (Revision 1 argued this from aggregate retention and from the
  10,021 to 13,691 difference; both were confounded and are withdrawn.)

Root cause.

1. Path construction absorbs every contained-pattern transfrag
   (`rlink.cpp:7792`) and joins others, consuming them.
2. The prediction cleanup removes any prediction whose introns are a subset of
   another's (`rlink.cpp:7866` to `7915`), compares only reconstructed `cov`,
   and does not merge `merge_files` into the survivor.
3. `--path-mode` is checked on the finished path in `printMergeResults()`
   (`rlink.cpp:20220`), after the inputs are gone.
4. Neither absorption site calls `mergeReportEvent()`.

Plan. Section 4.1, steps 2 to 4. Piece support is the distinct-file count of
the path's weakest intron, already computed by
`set_merge_prediction_support()` (`rlink.cpp:7110`).

Reasoning.

- The dominance rule needs no new option and uses a quantity the branch
  already computes. A join of two fragments each seen in 10 files has piece
  support 10 and still absorbs them, so legitimate coassembly is unaffected. A
  single-file long model has piece support 1 and can no longer erase an
  isoform seen in 10 files.
- Lost transcripts are the largest class of open merge complaints (#370, #267,
  #402, #487, #461, and the closed #491 and #66).

Acceptance. Fixture Q emits both models with supports 1 and 3; Q2 (equal
support) emits one model by default and two with `--max-intron-overhang 0`.
Fixture W
emits both originals under `direct` and the joined path under `shared`. The
invariant of section 4.1 holds on the dataset with 0 exceptions. The 3,043
chains with no containing output path are either emitted or carry a report
row.

### F5. `contained_set` must respect strand (regression)

Problem. A generated SET lying inside an exon of any retained MET is rejected,
including a MET on the opposite strand.

Evidence. Fixture L: plus-strand MET with exon 2000-2600 and a minus-strand SET
2100-2500. New code rejects the SET (`contained_set`); 3.0.3 keeps it.

Root cause. The loop at `rlink.cpp:20294` has no strand test. The removed
`equal_strand()` helper used to provide it.

Plan. Require `pred[m]->strand==pred[n]->strand`, or the SET being unstranded
(`.`), before declaring containment.

Reasoning. Antisense single-exon transcripts are real and common in compact
and bacterial genomes; issues #105 and #468 are about wrong cross-strand
handling. One-line fix, independent of section 4.1.

### F6. Unstranded inputs disappear (regression)

Problem. Unstranded inputs still go through the legacy path, which splits a
MET into exons and unions overlapping SETs. The new support lookup then finds
no input with that exact structure and rejects the result with support 0.

Evidence.

- Fixture N: unstranded MET `8000-8100,9000-9100` gives no output and two
  report rows (`8000-8100` and `9000-9100`, `reject minimum_support`, support
  0). 3.0.3 emits the two exons as unstranded SETs. (Revision 1 said "no
  report row"; that run lacked `--report`.)
- Fixture U (second reviewer): unstranded SETs 1000-1300 and 1200-1500 give
  no output and one row `1000-1500 reject minimum_support 0`. 3.0.3 emits the
  union 1000-1500.

Root cause. Unstranded bundle output at `rlink.cpp:16765`;
`merge_fill_prediction_support()` at `rlink.cpp:20046`.

Plan (decided: pass-through). Unstranded inputs bypass graph construction.
Each is emitted with its own structure and strand `.`, deduplicated by the
identity of section 4.2 (intron chain, or reciprocal overlap for SETs),
counted and filtered by `--min-support` like stranded ones. No splitting and
no union.

Consequence to document. Fixture U yields two SETs, not the 3.0.3 union,
because their overlap is below `--set-overlap`. That is the same policy the
branch applies to stranded SETs. Related: #224.

### F7. Test suite

1. Negative assertions are no-ops. In POSIX sh, `set -e` is ignored for a
   pipeline that begins with `!`, so `! has_tx ...` never aborts the script.
   All 9 such lines in `test_merge_policy.sh` are affected. They happen to
   hold today, but a deliberately false one still prints `OK`. Replace them
   with a helper:

   ```sh
   no_tx() { if has_tx "$@"; then echo "unexpected transcript: $*" >&2; exit 1; fi; }
   ```

2. `make test` does not run the merge tests, and `run_tests.sh` never exits
   non-zero on a mismatch (`#exit 1` is commented out). Fix both.
3. Assert complete exon chains and support attribute values, not only outer
   bounds.
4. Add the fixtures of section 8, each with a minus-strand mirror where an end
   is involved. Fixtures for defects not yet fixed go in as expected failures
   so the suite stays green and each fix flips one.
5. Add coverage for what the plan listed but the script lacks: `-T`, within-
   file `-f`, `--set-overlap`, `--slop-bp-*`, `--slop-frac-*`,
   `--rescue-support`, `--path-mode direct` and `replicated`, missing-abundance
   notices, option-validation failures (assert a non-zero exit), and thread
   determinism (sorted structure set and report identical at `-p 1` and
   `-p 4`).
6. Keep the dataset checks of section 7 as scripts next to
   `compare_merge_outputs.pl`; the dataset is too large for the repository, so
   this remains a manual pre-release step.

### F8. Smaller items

- Numeric options are parsed with `asInt()`/`asDouble()`. `--slop-frac-3p x`
  and `--min-introns abc` become 0 (the latter disables the intron gate);
  `--min-support 2x`, `--set-overlap nan` and `--slop-frac-3p nan` are also
  accepted. Parse strictly (whole string, finite, in range) and name the
  offending option in the error.
- A guide transcript named like a generated ID (`MSTRG.1.1`) collides with the
  generated one. Also present on master. Skip IDs used by guides.
- `VERSION` is still 3.0.3 (`stringtie.cpp:14`) although merge defaults and
  the output schema changed. Bump it and add a short "behavior changes" note
  with the options that approximate the old behavior.
- `worklogs/merge-issues.md` still reports 144,781 output transcripts; the
  current figure is 192,173. Correct the README statements that are not true
  today ("Hard caps always win", "at most one vote to an exact structure") and
  rewrite the option section around the layers of section 4.2.
- Output order: 7,487 out-of-order transcript lines at `-p 1` versus 1,337 in
  3.0.3, because changed guides are appended by `writeUnbundledGuides()`. Not
  blocking; relevant to #177.
- `writeUnbundledGuides()` rescans all guides of a reference sequence twice
  per guide under `--dedup-guides`; use the `mergeGuideStructures` and
  `mergeGuideAliases` hashes that `registerMergeGuide()` already fills.
- `mergeGeneratedStructures` is process-wide but a structure cannot recur
  across bundles; a bundle-local hash is enough and does not grow with the
  genome.
- A 588 MB `core` file (from codex, not StringTie) sits untracked in the
  repository root.

## 5. `--min-support` default

Decision. The `auto` proposal of revision 1 (1 below 4 files, else 2) is
withdrawn: with it, adding an unrelated fourth file would change filtering
everywhere. The branch keeps 1 while sections 4.1 and 4.2 are implemented.
The default is then evaluated as a separate step, and the release that changes
merge defaults is held until that evaluation, so that defaults do not change
in two consecutive releases. Merging to master and releasing are separate.

Why keeping 1 is not neutral. Together with `-F 0 -T 0 -f 0` it is a default
change from 63,959 to 192,173 records on this dataset, 77% of the generated
ones single-file (section 3.4), validated on one dataset only.

What this dataset says (not sufficient to decide):

1. Noise is concentrated in single-file structures: 65% of novel chains and
   46% of SET records occur in one file; single-file SETs average 385 bp
   against about 1,000 bp for replicated ones.
2. Novel chains at support 1 grow linearly with the number of files, about
   8,000 per file, while reference chains saturate by 4 files (section 3.8).
3. The sensitivity cost of 2 by the only truth proxy available is 473 of
   15,551 reference chains (3.0%) at N = 20, but 28% at N = 2. The proxy is
   biased upward by guided assembly.
4. Abundance is not a substitute: `-T 1` removes 68,538 records without
   recovering one additional reference chain (10,019 versus 10,021).

Surveys needed before choosing:

```text
Id  Survey                                         Decides
--  ---------------------------------------------  ------------------------------
S1  Subsampling ladder on the alga dataset         Done (section 3.8)
S2  Second dataset, no -G at assembly: merge       Precision and sensitivity of
    without a guide and score intron chains        min-support 1, 2, 3 with and
    against GENCODE (gffcompare); compare with     without -T 1, against 3.0.3
    3.0.3
S3  Same dataset, merge with -G from which about   Recovery of real transcripts
    20% of transcripts were removed; score         the guide lacks, and how many
    recovery of the removed ones                   guide variants are emitted
S4  Repeat sections 3.1 and 3.5 to 3.8 on the      End variability and slop
    second dataset, polyA and total RNA apart      defaults where ends are not
                                                   copied from a guide; intron
                                                   overhang; 3-prime readthrough
S5  Small merges: subsets of 2 and 3 files, and    Small-N behavior of any
    reference plus 1 or 2 GTFs (#247, #276)        default
S6  Re-run sections 3.2 to 3.5 on the fixed code   Actual cost of dominance and of
                                                   chain-keyed support
S7  Long-read merge of 2-3 files (#461 shape)      No data yet
```

Second dataset (in preparation). Human brain, 10 samples each sequenced as
polyA and as RiboZero total RNA (20 CRAMs, hg38). Assembled with StringTie
3.0.3 without `-G`: polyA with defaults, RiboZero with `-N`, and RiboZero
with defaults as a control. Output in `merge-dev/brainseq_eval/` (ignored by
git). S2 to S5 can run on the current branch for the "before" numbers and
again after section 4 is implemented. Caveats to carry: 10 individuals of one
tissue, so sample-private isoforms are under-represented compared with a
multi-tissue design; GENCODE is the truth proxy, so real unannotated
transcripts count as false.

Candidates to compare in S2: 1; 2; 2 for SETs with 1 for METs; any of these
with `-T 1`. A fraction of N is not a candidate: a tissue-specific isoform
appears in a small fixed number of samples whatever the cohort size (#487).

## 6. Order of work

1. Tests first: F7 items 1 to 4 (failure propagation, exact-structure and
   support assertions, fixtures as expected failures).
2. Independent fixes: F5 (strand), F8 option validation, F6 (unstranded
   pass-through).
3. Layer 1 identity (F2), needed by everything after it.
4. Restructure `merge_transfrags()` (section 4.1) with end election (F1),
   guide relation (F3), dominance and accounting (F4).
5. Policy switches of section 4.3, `--end-variants replicated` last.
6. Surveys S2 to S6 on the fixed code (S2 to S5 also on the current branch as
   the baseline); remaining F7 and F8 items; README rewritten around sections
   4.2 and 4.3.
7. Separately: choice of defaults from the surveys, version bump, release.

## 7. Acceptance criteria

On the issue #247 dataset at the branch defaults:

```text
Check                                                Now            Required
---------------------------------------------------  -------------  ---------------
Guide IDs with exact guide structure                 26,059/26,059  unchanged
Duplicate transcript IDs                             0              0
Emitted ends consistent with the endpoint policy     not checkable  0 exceptions
(F1 acceptance)
Input chains passing layer 1 that are neither in     at least       0
the output nor in a report row satisfying dominance  20,843
Input records accounted for by output plus report    not possible   927,737/927,737
member counts
--min-support 2: replicated novel chains kept or     3,178 kept of  65,636/65,636
reported as absorbed                                 65,636
Structure set and report identical at -p 1, 4, 8     yes            yes
```

"At least 20,843" is 17,800 replicated chains under a weaker path plus 3,043
with no containing path (no-guide run; the with-guide figure is to be measured
in S6).

Plus: all fixtures of section 8 pass with effective negative assertions,
valgrind-clean; `run_tests.sh` 9 of 9 and exiting non-zero on a mismatch;
`make test` runs the merge suite; runtime and peak memory on the full dataset
reported against the current branch.

## 8. Fixtures

All on one reference sequence, plus strand unless noted. "x3" means the same
model supplied in 3 input files. `full` is 1000-1100, 2000-2100, 3000-3100,
4000-4300.

```text
Id  Inputs                                          Now                 Expected
--  ----------------------------------------------  ------------------  ------------------------
A   guide=full; 2050-2100,3000-3100,4000-5300 x1    guide only,         guide, credited support,
                                                    no report           trim row
D   A with --min-support-3p 1                       guide only          guide + extended variant
E   A with the input x3                             guide only          guide + extended variant
F   guide=full; full with last exon to 5300 x1      guide + variant     same
G   full x2; partial of A x1                        1000-5300           1000-4300, trim row
H   full x2; full with last exon to 5300 x1         1000-5300           1000-4300, trim row
H2  H with --max-extension-3p 100                   1000-5300           1000-4300
I   full x2; full with first exon from 300 x1       300-4300 and        1000-4300, trim row
                                                    1000-4300
K   first exon from 300 x1; last exon to 5300 x1    300-5300            one model, ends per F1
S   1000-1100,2000-2100,3000-3300 with ends         empty at            one model at
    shifted by -10, 0, +10 bp, one per file         --min-support 2     --min-support 2 and 3
Q   1000-1100,2000-2100,3000-3300 x1;               3-exon model only,  both models,
    1950-2100,3000-3300 x3                          support 1           supports 1 and 3
Q2  Q with the 3-exon model x3                      3-exon model only   default: 3-exon model,
                                                                        absorb row, support 3;
                                                                        --max-intron-overhang 0:
                                                                        both models
L   + MET 1000-1100,2000-2600,3000-3100;            SET rejected        both kept
    - SET 2100-2500
N   unstranded 8000-8100,9000-9100                  nothing; 2 rows     unstranded MET, unchanged
                                                    support 0
U   unstranded SETs 1000-1300 and 1200-1500         nothing; 1 row      both SETs, unchanged
                                                    support 0
V   guide 1000-1199,1300-1499;                      MSTRG.1.1, no       absorbed at default;
    input 900-1199,1300-1399                        reference relation  variant, dropped, or
                                                    under any option    capped per F3 otherwise
W   1000-1100,2000-2100,3000-3300 x1;               direct: nothing     direct: both originals
    2000-2100,3000-3300,4000-4300 x1;               shared: 1000-4300   shared: 1000-4300
    --min-support-5p 1 --min-support-3p 1
X   --min-support 2x, --set-overlap nan,            exit 0              non-zero exit naming
    --slop-frac-3p nan or x, --min-introns abc                          the option
```

For G, H and I the expected end assumes the default `--min-support-5p/3p 2`
and an extension beyond the slop.

The scripts used for sections 3 and 7 (input support statistics, chain
retention, end check, absorber classification, end survey, subsampling ladder)
exist only in the review session's scratch directory and are not part of the
repository. They are short Perl scripts over transcript and exon records;
recreate them alongside `compare_merge_outputs.pl` as part of F7 item 6. The
second reviewer's artifacts are listed in `review-branch-response1.md` and are
likewise temporary.

## 9. Decisions

Decided by the maintainer (2026-11-09):

1. Dominance rule accepted, including about 17,800 additional transcripts on
   this dataset.
2. `merge_transfrags()` is restructured to "observed chains first"
   (section 4.1) rather than patched site by site.
3. Options must be consistent and compound: intron-chain comparison is the
   base, boundaries are layered on top, each option only restricts
   (section 4.2).
4. `--min-support` default: not `auto`; decide after more surveys (section 5).
5. Unstranded inputs pass through unchanged (F6).
6. Fuzzy ends relative to a guide are tolerated up to a user-tunable threshold,
   working value 100 bp. A partial chain extending past a guide end by less
   than that counts as the guide. This applies to guide end coordinates only,
   not to intron coordinates.

Decided by the maintainer (2026-11-09, second round):

7. No policy choice is hard-wired: end variants, guide tolerance, guide
   fragment handling and intron overhang are command-line switches
   (section 4.3).
8. Intron overhang: absorbed by default, with an option to keep such
   fragments strictly separate.
9. The 3-prime slop default is larger than the 5-prime one. Unbounded was
   suggested; 500 is provisional until the second dataset is surveyed.
10. Second dataset: the brain polyA and RiboZero samples, assembled without
    `-G` (section 5).

Still open:

1. Final defaults for `--slop-bp-3p` (500, another value, or unlimited) and
   `--max-intron-overhang`, from survey S4.
2. `--min-support` default, from surveys S2, S3 and S5.
3. Option names in section 4.3.
4. A long-read dataset for S7, if that case should gate the release.
5. Version number for the release that changes merge defaults.

## 10. Implementation status (revision 4, 2026-11-10)

Everything below is in the working tree on top of `a27c05d`, uncommitted. A
copy of the diff is in `merge-dev/brainseq_eval/survey/working_tree_checkpoint.patch`
(taken before the fixes of section 12; `git diff` is authoritative).

```text
Item                                   State
-------------------------------------  ----------------------------------------
F7 tests                               Done: no_tx, exact exon-chain and support
                                       assertions, fixtures of section 8 plus
                                       review fixtures; run_tests.sh exits
                                       non-zero on a mismatch; make test runs
                                       the merge suite
F5 strand in contained_set             Done
F8 strict numeric options              Done (whole value, finite, in range)
F6 unstranded pass-through             Done
F2 intron-chain support                Done (MET chain key; SET pairwise overlap
                                       in the filter and in output support)
4.1 observed chains first              Done, see "Departures" 1
F1 end election                        Done, see "Departures" 2
F3 guide relation                      Done, see "Departures" 3
F4 dominance, credit, absorb rows      Done
4.3 --guide-slop-bp-5p/3p,             Done
    --guide-fragments,
    --max-intron-overhang,
    "unlimited" for --slop-bp-*
4.3 --slop-bp-3p default 500           Done (provisional)
4.3 --end-variants                     Not implemented
F8 remaining items                     Not done: ID collision with a guide named
                                       like a generated ID, version bump,
                                       writeUnbundledGuides scans, stale number
                                       in merge-issues.md
Accounting of all input records        Partial: absorb rows carry record and file
                                       counts; no summed reconciliation yet
README.md, MERGE_MODE.md, -h text      Updated to the implemented behavior
```

Verification.

- `make test`: 9 of 9 assembly tests and the merge suite pass; the merge suite
  has no known failures left and is valgrind-clean.
- Alga dataset with guide: 193,102 transcripts; 26,059 of 26,059 guide IDs with
  their exact structure; 0 duplicate IDs; structure set and report identical at
  `-p 1` and `-p 8`; 25 s at `-p 1` (the unmodified branch takes 34 s at
  `-p 8`). Guide variants emitted: 13 with the guide's chain, 1,595 with part
  of it; 12,543 guides carry `merge_contained_support`. Out-of-order transcript
  lines at `-p 1`: 8,587 (was 7,487); not addressed.
- Human polyA, no guide: input chains lost only under a weaker path fell from
  11,629 to 58 (section 11.3).

Departures from the plan, decided while implementing.

1. `merge_transfrags()` itself was not rewritten. `build_merge()` records every
   observed stranded intron chain before reads are collapsed, and
   `printMergeResults()` makes one candidate per chain, keeps graph paths only
   when their chain was not observed (joins), and runs one redundancy pass
   after all per-prediction rejections. The legacy cleanup loop still runs but
   can no longer remove an observed chain. Same design, smaller change.
2. Ends are elected from the files that carry the chain itself. A first
   version pooled every input sharing the terminal intron; on the alga data
   that produced about 16,000 guide variants whose end combination was never
   observed, so it was dropped. Consequence: a chain seen in one file keeps its
   own ends (it is governed by `--min-support`), and a chain seen in several
   files gets the end that `--min-support-5p/3p` files reach. For a chain
   without support (a join) the hard cap bounds the emitted end against the
   elected end, not against any outside reference.
3. An end that passes a guide end, or the end of a containing transcript, by
   more than the slop is a real extension only if `--min-support-5p/3p` files
   of that chain reach that far. Otherwise the input is a fragment and is
   absorbed (credited, reported). This follows F3 plan item 2 but changes a
   behavior the branch tested: a single-file endpoint variant of a guide is no
   longer emitted by default; it needs `--min-support-5p 1` or `-3p 1`.
   Fixture F of section 8 therefore now gives the guide only. Existing tests
   that expected such variants were given explicit options.
4. Fuzz within `--guide-slop-bp-*` (default 100 / 500) is credited to the
   guide's `merge_support`; fragments go to `merge_contained_support`. Guide
   support is now derived from the observed chains, not from graph pooling.
5. Guides whose ends the graph had moved are restored and emitted in place.
6. With `--min-introns 0` a terminal exon may still be extended by a
   single-exon input (the legacy behavior that option asks for); only ends
   beyond every observed end of the chain are taken from the graph.
7. A SET that overlaps a single-exon guide within the guide slop is absorbed
   into the guide (second code review, finding 6).
8. `-i` keeps retained-intron variants; without it they are absorbed under the
   same dominance rule as other fragments.

Side effect to know about: a crash during development (mixed debug and release
objects after a header change) made the kernel overwrite the stray 588 MB
`core` file in the repository root with a 1 MB one. Core dumps were disabled
for all later runs. Build with `make release` after header changes; plain
`make` produces debug objects that do not link safely with the release gclib
objects.

## 11. Second dataset: human brain, polyA and RiboZero

Ten individuals, each sequenced as polyA (unstranded) and RiboZero
(fr-firststrand), hg38. Assembled with StringTie 3.0.3 without `-G`: polyA with
defaults, RiboZero with `--rf -N`, RiboZero with `--rf` only as a control. The
reference was verified against every CRAM header (names, lengths, M5), and
strandedness was established from the data. Scored against GENCODE 49 (475,887
distinct intron chains) with gffcompare at the intron-chain level. Files are in
`merge-dev/brainseq_eval/` (git-ignored); scripts and outputs in `survey/`,
results summary in `survey/RESULTS_polyA.txt`.

Caveats: one tissue, 10 individuals; GENCODE as truth counts real unannotated
transcripts as false; absolute sensitivity is low because one tissue expresses
a small part of GENCODE.

### 11.1 Input chains by support (independent of the merge code)

Share of chains at each support level that match a GENCODE chain:

```text
Files   polyA (10)   RiboZero -N (10)   all 20
-----   ----------   ----------------   ------
1             9.3%               7.4%     5.9%
2            23.7%              21.2%    15.5%
3            35.3%              33.0%    23.5%
5            53.8%              51.2%        -
10           93.7%              94.9%        -
```

Cumulative, polyA: 54,316 matching of 271,077 chains at 1 or more files
(20.0%); 36,453 of 78,478 at 2 or more (46.4%); 29,206 of 47,869 at 3 or more
(61.0%). Requiring 2 files removes 33% of the matching chains, not 3% as the
guided alga data suggested (section 5 item 3): that figure was an artifact of
guided assembly.

Chains seen in 2 or more files that have one base-identical structure in 2
files: 5.2% of matching and 1.2% of novel chains (polyA). With the old
exact-ends key `--min-support 2` kept 2,487 transcripts out of 78,478 eligible
chains. This settles F2.

Subsampling: novel chains at 1 file grow linearly with the number of files
(polyA: 52,634 at N = 2, 216,761 at N = 10); precision of "2 or more files"
falls with N (78.8% at N = 2, 46.4% at N = 10, 34.3% at N = 20 mixed), so a
fixed threshold is not equivalent across cohort sizes.

### 11.2 No-guide merges scored against GENCODE

Transcripts / matching intron chains / sensitivity / precision.

```text
Run                        polyA                       RiboZero -N
-------------------------  --------------------------  --------------------------
3.0.3 default               42,210  20,360  4.3  49.5   33,279  15,612  3.3  48.5
3.0.3 -F 0 -T 0 -f 0       183,515  41,180  8.6  22.7  198,970  34,121  7.2  17.5
a27c05d default            212,004  43,132  9.0  21.3  236,550  36,504  7.7  16.5
a27c05d --min-support 2      2,487   1,875  0.4  79.8    1,984   1,448  0.3  80.2
a27c05d -T 1               119,513  33,213  7.0  30.0  132,045  27,879  5.8  23.7
now, --min-support 1       205,608  45,127  9.5  23.1  230,984  39,167  8.2  18.2
now, --min-support 2        70,218  33,021  6.9  49.3   62,908  28,076  5.9  47.1
now, --min-support 3        45,264  27,303  5.7  62.8   38,628  23,002  4.8  62.1
```

```text
Run                        RiboZero default            all 20 (polyA + RiboZero -N)
-------------------------  --------------------------  --------------------------
3.0.3 default               38,303  15,964  3.3  42.9   51,005  20,438  4.3  41.1
a27c05d default            277,667  36,260  7.6  14.0  370,024  49,075 10.3  14.1
now, --min-support 1       270,803  38,752  8.1  15.5  354,629  52,348 11.0  15.8
now, --min-support 2        66,236  27,819  5.8  44.6  105,924  39,195  8.2  39.2
now, --min-support 3        39,458  22,872  4.8  60.6   67,215  33,042  6.9  51.6
```

Reading.

- In all four sets the implemented code at `--min-support 3` is better than the
  3.0.3 default on both axes (polyA: 34% more matching chains, precision 62.8
  against 49.5). At `--min-support 2` it finds 62 to 92% more matching chains at
  about the same precision as 3.0.3.
- At `--min-support 1` precision is 16 to 23%: four of five output transcripts
  do not match GENCODE, and single-file chains match it 6 to 9% of the time.
- Observed-chains-first with dominance improved both axes at every setting
  over chain-keyed support alone (polyA at 2: 31,991 to 33,021 matches,
  precision 47.3 to 49.3).
- `-T 1` helps here (precision 21 to 30 on polyA) but less than replication.
- `-N` for RiboZero assembly gives fewer and cleaner inputs than the default
  (precision at 2: 47.1 against 44.6).

### 11.3 Lost chains and hold-out

Input chains missing from the no-guide polyA output whose best containing
output path is weaker than they are: 7,295 / 3,072 / 1,262 (support 2-4, 5-9,
10) on `a27c05d`; 50 / 7 / 1 now.

Hold-out, polyA merged with `-G` from which 20% of GENCODE transcripts were
removed (94,819 hidden chains, 10,704 of them present in the inputs):

```text
Run                               Hidden recovered   Novel chains   Hidden share
--------------------------------  ----------------   ------------   ------------
3.0.3 default                     3,689 (34.5%)            25,111          14.7%
a27c05d default                   7,406 (69.2%)           129,740           5.7%
now, default                      7,148 (66.8%)           121,325           5.9%
now, --guide-fragments dominance  7,890 (73.7%)           132,811           5.9%
now, --min-support 2              5,087 (47.5%)            27,984          18.2%
now, --min-support 3              4,189 (39.1%)            15,409          27.2%
```

At `--min-support 2` the merge recovers 38% more hidden real transcripts than
3.0.3 with a similar number of novel models; at 3 it recovers 14% more with
40% fewer novel models. The default loses 258 hidden transcripts against
`a27c05d` because guides absorb their fragments regardless of support;
`--guide-fragments dominance` recovers them and more, at the cost of 11,500
additional novel chains.

### 11.4 Ends and intron overhang

Share of replicated novel chains (3 or more files) whose most extended file is
within a distance of the second most extended file:

```text
Set            End  100 bp  500 bp  1000 bp
-------------  ---  ------  ------  -------
polyA          5p      73%     90%      95%
polyA          3p      69%     85%      91%
RiboZero -N    5p      72%     88%      93%
RiboZero -N    3p      60%     78%      85%
RiboZero def.  3p      57%     75%      83%
```

Assembled ends against GENCODE ends for the same chain (polyA): 5-prime within
100 bp 70%, within 500 bp 92%; 3-prime within 100 bp 53%, within 500 bp 69%;
13% of records pass the annotated 3-prime end by more than 1 kb (15% in
RiboZero -N).

Fragment ends next to an internal exon of a GENCODE transcript: 54 to 57%
(polyA) and 57 to 67% (RiboZero) overhang into the intron or span it; the
overhang exceeds 1 kb in 5 to 11% (polyA) and 11 to 18% (RiboZero) of cases.

Reading for the defaults.

- 3-prime ends are less reproducible than 5-prime ends, and clearly less so in
  RiboZero. A 500 bp 3-prime slop keeps the longest end for 78 to 85% of
  replicated chains. An unbounded default would let a single file set the
  3-prime end of the remaining 15 to 22%, most of them more than 1 kb long in
  RiboZero. I recommend keeping 500 as the default and offering `unlimited`.
- Intron overhang is the norm for fragments, as the maintainer expected.
  Absorbing by default (`--max-intron-overhang unlimited`, subject to
  dominance) is consistent with the data.

## 12. Second code review (independent agent, of the implementation)

An independent review of the working-tree diff confirmed, with fixtures, the
defects below; all are fixed and each has a regression fixture.

```text
Finding                                               Fix
----------------------------------------------------  ---------------------------
Unstranded pass-through used read indices that go      Keep pointers; fixture R1
stale when a guide is in the bundle (regression)
--min-introns 0 bypassed end election for every        Only ends beyond all
observed chain                                         observed ends count; R5
Guide piece support overwritten with 0, which broke    Do not refill support that
--guide-fragments dominance                            was already derived; R8
SET support missed identical records collapsed by      Credit collapsed records;
build_merge                                            R2
Single-exon inputs near a single-exon guide were not   Absorb within guide slop;
treated as the guide                                   R3
```

Reported and left as they are:

- End election is not monotone in the number of files: one file keeps its own
  end; adding a second file that ends earlier shortens the model. This is the
  stated rule (a single-file chain is governed by `--min-support`).
- A join is always absorbable by a longer container and is never a guide
  variant of the partial-chain kind, whatever supports its ends.
- If `--guide-policy best` or a later duplicate check removes a guide, the
  fragments it absorbed lose their credit (the structure is still emitted).
- With `--min-introns 0`, `trim` rows can describe an end as not adopted when a
  single-exon join later extends the transcript further.

The review found no memory, thread-safety or determinism problem in the new
code; valgrind is clean on the suite and on its extra fixtures.

## 13. Decisions for the maintainer after revision 4

Item 1 is settled (2 for 4 or more files, otherwise 1) and implemented; see
section 14, which also has new measurements for items 3 and 4.

1. `--min-support` default. The data now argue against 1: on 10 human samples
   it gives 16 to 23% precision, and 2 or 3 beat 3.0.3 on the measures
   available. Against a fixed 2: with 2 or 3 input files it removes most real
   transcripts (at N = 2, 46% of matching chains). My recommendation: default 2
   when 4 or more files are given and 1 otherwise is still the best fit to the
   data, but the objection that the threshold then depends on N stands
   (section 5). The alternative that avoids it is to keep 1 and print a
   prominent note when 4 or more files are merged. This needs your call.
2. Departure 3 of section 10: single-file endpoint variants of a guide are no
   longer emitted by default. Confirm, or ask for the old behavior as default.
3. `--guide-fragments`: `absorb` (default, fewer transcripts) or `dominance`
   (more recovered real transcripts, section 11.3).
4. `--slop-bp-3p` 500 as the default, or `unlimited`.
5. `--end-variants` (not implemented): still wanted?
6. Whether to commit the working tree as it is, split into commits per item,
   or revise first.

## 14. Follow-up after the maintainer's answer on decision 1 (2026-11-10)

Decision 1 is settled and implemented: without `--min-support` the threshold
is 2 when 4 or more files are merged and 1 otherwise. The value used is in the
`# min_support` line of the report and in the `-v` log. Fixtures Sd3, Sd4,
Sd41. Alga dataset with guide under the new default: 70,972 transcripts
(193,102 at `--min-support 1`), 26,059 of 26,059 guides exact, 0 duplicate
IDs, structures and report identical at `-p 1` and `-p 8`.

### 14.1 Decision 4: 3-prime slop (no guide, `--min-support 2`)

3-prime end of each output transcript against the MANE Select transcript with
the same intron chain (11,300 to 12,300 chains per set). "Over" and "short"
are more than 1 kb longer or shorter than MANE. Script
`scripts/end3ref.pl`, runs in `*_slop3/`, table in `slop3.summary.txt`.

```text
                       polyA                 RiboZero -N           all 20
Setting                <=500  over  short    <=500  over  short    <=500  over  short
---------------------  -----  ----  -----    -----  ----  -----    -----  ----  -----
--slop-bp-3p 100       76.7%  7.1%  10.3%    72.2%  9.8%  11.2%    72.7% 12.3%   8.2%
--slop-bp-3p 500       76.4%  7.4%  10.2%    71.4% 10.2%  11.1%    71.7% 12.7%   8.1%
--slop-bp-3p 2000      75.3%  8.3%  10.0%    69.9% 11.7%  10.8%    69.9% 14.4%   7.9%
--slop-bp-3p unlimited 75.3%  8.3%   9.9%    69.8% 11.9%  10.7%    69.9% 14.5%   7.8%
 + --slop-frac-3p 1    75.1%  9.1%   9.4%    69.1% 13.3%  10.0%    68.9% 15.9%   7.2%
--min-support-3p 1     74.9% 11.5%   7.3%    68.4% 15.9%   8.0%    67.5% 18.8%   5.7%
```

Base-level precision and number of loci (gffcompare, all GENCODE):

```text
Setting                polyA            RiboZero -N      all 20
---------------------  ---------------  ---------------  ---------------
--slop-bp-3p 500       76.8%  27,727    65.2%  27,369    61.6%  33,466
--slop-bp-3p unlimited 76.2%  27,686    64.4%  27,328    61.0%  33,389
--min-support-3p 1     72.6%  27,294    59.0%  26,906    55.5%  32,613
```

Reading.

- `--slop-bp-3p unlimited` is not "the longest end wins": `--slop-frac-3p`
  (0.50, at most 1) still limits the extension to half of the elected terminal
  exon. The real "longest end wins" is `--min-support-3p 1`.
- 500 to unlimited costs about 1 to 2 points of ends within 500 bp of MANE,
  adds 1 to 2 points of ends more than 1 kb too long, and removes 0.3 points
  of truncated ends. Truncation relative to MANE is mostly present in every
  input file, so a wider slop does not repair it.
- Letting one file set the end (`--min-support-3p 1`) costs 4 to 6 points of
  base-level precision and joins 430 to 850 loci.
- Overshoot grows with the number of files under every setting (20 files
  against 10), faster under the wider ones.

### 14.2 Decision 3: `--guide-fragments` at `--min-support 2` (hold-out)

```text
Set          Mode       Hidden recovered  Novel chains  Hidden share  5+ files
-----------  ---------  ----------------  ------------  ------------  --------
polyA        absorb     5,087 (47.5%)           27,984         18.2%     75.9%
polyA        dominance  5,681 (53.1%)           33,982         16.7%     84.1%
RiboZero -N  absorb     4,129 (44.0%)           24,563         16.8%     74.5%
RiboZero -N  dominance  4,763 (50.8%)           31,252         15.2%     84.4%
```

`dominance` recovers 12 to 15% more hidden transcripts and emits 21 to 27%
more novel chains; of the additional chains 9.5 to 9.9% are hidden reference
transcripts, about half the rate of the other novel chains. `absorb` loses
about 8 to 10% of hidden transcripts that 5 or more files support. Runs in
`gfrag_ms2/`.

## 15. Remaining items closed (2026-11-10, version 3.0.4)

Decisions 2 to 5 of section 13 were confirmed as implemented (single-file guide
end variants absorbed, `--guide-fragments absorb`, `--slop-bp-3p 500`, no
`--end-variants`). Closed afterwards:

- `VERSION` is 3.0.4; README has a "Changes in merge mode since version 3.0.3"
  table with the closest option for each change.
- Output order: merge records go to a temporary file with an in-memory index
  and are written sorted by reference, start, end. References keep the order
  of the input files (#177). Unbundled guides are written in place. Human
  polyA with guide: 50,590 out-of-order transcript lines before, 0 now; alga:
  0.
- Generated IDs are numbered in output order: output is byte-identical at any
  `-p` (alga `-p 1`/`-p 8`, human `-p 1`/`-p 6`).
- A guide ID or gene ID of the form `<label>.N[.M]` reserves gene number N.
- `writeUnbundledGuides()` uses the guide hashes instead of rescanning;
  the generated-structure hash is per bundle.
- `worklogs/merge-issues.md`: 70,972 (default) and 193,102 (`--min-support 1`).
- The stray `core` file was removed.

Structure sets and reports are unchanged by these edits on both datasets.

Input record reconciliation (added after the push, header only): the report
header has `# input_records`, one `# fate <name> <count>` line per final fate,
and `# fate unaccounted`. Unaccounted is 0 on the alga run (927,737 records)
and on the human polyA and RiboZero -N runs with and without guide; every
fixture of the merge suite asserts 0. Multi-exon records are traced through
the decision for their chain; stranded single-exon records are classified by
the output transcript that covers them. Observation from the tally: with the
GENCODE guide, 78,690 of 572,057 polyA input records (13.8%) are rejected as
`cross_gene_bridge`.

## 16. The cross-gene barrier on GENCODE (2026-11-10, not changed)

`filter_merge_inputs()` rejects an input record whose exons overlap exons of
guide transcripts from two or more gene IDs on the same strand
(`cross_gene_bridge`), unless `--allow-fusions`. The record is discarded: it is
not emitted and supports nothing. Unstranded records never match a stranded
guide and are exempt. Human polyA, 10 files, `-G guide80.gtf`
(`scripts/xgene.pl`, `fate_check/xgene.txt`):

```text
Rejected records                                   78,690  (13.8% of 572,057)
  the guide genes overlap each other               67,945  86.3%
  the guide genes are at separate loci             10,745  13.7%
  a GENCODE readthrough gene is involved           34,338  43.6%
  own intron chain is a guide transcript           24,888  31.6%
  own intron chain is a hidden GENCODE transcript   5,889   7.5%
  separate loci and chain not annotated             8,955  11.4%
```

RiboZero -N: 75,910 rejected, same proportions. With `--allow-fusions`
(default `--min-support 2`): polyA hidden transcripts recovered 5,087 ->
5,721, novel chains 27,984 -> 33,162; RiboZero -N 4,129 -> 4,695 and
24,563 -> 28,878. The barrier as written mostly hits annotated overlapping
genes, not bridges between separate genes. Alga data: 4.1% of records.

## 17. Annotation-aware barrier and class codes (2026-11-10, uncommitted)

Maintainer direction: the barrier must follow the annotation. Implemented:

- Meta-genes: guides whose spans overlap on a strand, or share a gene ID, are
  one meta-gene (`merge_guide_metagenes()`); `cross_gene_bridge` needs exon
  overlap with two meta-genes.
- An input with a guide's intron chain is never rejected as a bridge, so it is
  credited to the guide.
- gffcompare class codes on generated transcripts (`merge_class_code()`, gclib
  `getOvlData()` strict, `classcode_rank()`): `cmp_ref`, `class_code`, and
  `ref_gene_id`/`gene_name` for same-strand exon overlaps. `reference_id` and
  `=` only on guides.

```text
                              polyA                RiboZero -N
                              before    after      before    after
cross_gene_bridge records     78,690   11,880      75,914   12,414
  with a guide's own chain    24,888        0      20,878        0
hidden transcripts recovered   5,087    5,705       4,129    4,671
novel chains                  27,984   32,154      24,563   28,010
(--allow-fusions: recovered    5,721 / novel 33,162; 4,695 / 28,878)
```

Alga: 38,353 -> 24,154 rejected records, 70,972 -> 71,284 transcripts, guides
exact, output identical at `-p 1`/`-p 8`, 0 unaccounted.

Class codes against `gffcompare --strict-match -e 0` on the polyA merge
(42,644 generated transcripts): same code 98.2%, same reference 93.9% of
those, same gene 99.9%. Differences: priority among `i`/`x`/`y`, and
single-exon pairs, where gclib `singleExonTMatch()` in this tree ends with
`return 0` and never reports `=`/`~` (also on master; not changed).
849 of 6,729 class `k` transcripts refer to a single-exon guide (for example
a host-gene isoform that contains an embedded small RNA gene): that is the
gffcompare priority (`k` before `j`), kept as asked.

Follow-up (same day): the `singleExonTMatch()` fix was ported from upstream
gclib (commit 94c50ea, with the `cdsMatch` argument for `transcriptMatch()`,
`singleExonTMatch()` and the call in `getOvlData()`; `GFF_VERSION` 130). The
vendored `gff.cpp`/`gff.h` now differ from upstream only by `setTrackName()`
and one warning text. `k` against a single-exon guide now yields, for a
multi-exon transcript, to a splice-sharing relation with more exon overlap
(304 transcripts on the polyA merge). Agreement with
`gffcompare --strict-match -e 0`: same code 97.8% (98.6% without those 304);
what is left is the order of `i`, `x`, `y`.
