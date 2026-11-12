# StringTie codebase handoff (for work on master PRs)

Written 2026-11-10 at the end of a long session that reviewed and reworked
`--merge` on the `merge-dev` branch. It is for a new session that will work on
`master` without `merge-dev` merged. Statements are marked where they are
inferred rather than checked. Line numbers are for `master` at `d1dc38d`.

## 1. Where things stand

```text
Branch      Head     State
----------  -------  ---------------------------------------------------------
master      d1dc38d  StringTie 3.0.3. Also the merge base of merge-dev.
merge-dev   71e27c9  11 commits ahead, pushed to origin, no PR opened. Version
                     3.0.4. --merge behavior, docs and tests differ, plus a
                     gclib fix (see section 5).
```

- `merge-dev` is BOTH a branch and an untracked directory in the repo root
  (`merge-dev/`, datasets and survey output). Git commands that take a
  revision and paths need `--`: `git diff master merge-dev --`. The
  `rev:path` form is unambiguous: `git show merge-dev:AGENTS.md`.
- These exist only on `merge-dev`; read them from `master` with `git show`:
  `AGENTS.md` (code map, C++ conventions, glossary), `MERGE_MODE.md` (merge
  workflow), `.gitignore`, `test_merge_policy.sh`, `compare_merge_outputs.pl`,
  `worklogs/merge-issues.md`, `worklogs/merge-default-vs-3.0.3.md`.
- `master` has no `.gitignore`. Object files, `stringtie`, `htslib/` build
  products, `tests/`, `merge-dev/`, `worklogs/` all show as untracked there.
  Never `git add -A` or `git add .`; add files by name.
- `merge-dev` changes outside merge mode are small but real: `Makefile`
  (`make test` also runs the merge suite), `run_tests.sh` (exits 1 on a
  mismatch), `gclib/GSam.cpp` (4 lines, debug-only pointer reset),
  `gclib/gff.cpp`/`gff.h` (transcript match fix ported from upstream gclib),
  `bundle.h`/`rlink.h` (extra fields). A master PR touching those files will
  conflict with `merge-dev` later; keep such edits minimal.
- Open PRs against master seen on 2026-11-10 (`gh pr list`):
  #510 `-A` duplicate gene rows (1 line, `stringtie.cpp`),
  #503 `prepDE.py3` KeyError, #495 point features (`stringtie.cpp`),
  #482 `noexcept` in `gclib/GList.hh`.

## 2. Code map (master)

```text
File             Lines   Role
---------------  ------  ----------------------------------------------------
stringtie.cpp     1,782  CLI, options, guide loading, read loop, bundle queue,
                         worker threads, final output pass
bundle.{h,cpp}      242  BundleData, CReadAln, CJunction, CPrediction, guides
rlink.{h,cpp}    21,220  everything else: read processing, graphs, flow,
                         transcript inference, printing
tmerge.{h,cpp}      172  multi-input reader; --merge converts GTF to temp BAM
tablemaker.*        317  Ballgown tables (-B/-b) and read-count data
gclib/                   vendored copy of gpertea/gclib (not a submodule)
htslib/                  vendored htslib, built by htslib/build_lib.sh
prepDE.py, .py3          count-matrix helper scripts, independent of the C++
```

`stringtie.cpp`: `main` 301, `processOptions` 997, `processBundle` 1423,
`workerThread` 1554, `loadPtFeatures` 1648, `writeUnbundledGenes` 1692,
`writeUnbundledGuides` 1707.

`rlink.cpp` in call order:

```text
processRead            659   one alignment into BundleData (called by main thread)
infer_transcripts    17161   entry point per bundle (worker thread)
count_good_junctions 16656   junction filtering, coverage (normal mode)
build_graphs         14165   groups -> bundles -> splice graphs -> transcripts
  add_read_to_group   1633
  merge_fwd_groups    1246
  create_bundle       1084
  create_graph        3579
  get_fragment_pattern 4870, update_abundance 4680, process_transfrags 5768
  find_transcripts   13633   dispatch: guides_pushmaxflow 13079, parse_trf 11992,
                             get_trf_long 11070, get_trf_long_mix 10639
  push_max_flow       9075, long_max_flow 8864, bfs 7216
  store_transcript    9688
build_merge          16022   merge-mode counterpart of build_graphs
  merge_transfrags    7534, store_merge_prediction 7065/7142
printResults         20034   normal output; print_predcluster 18402 does the
                             final filtering (included_pred 17502,
                             retainedintron 17459/17472)
printMergeResults    19613   merge output
```

`bundle.h` structs: `GRefData` 38, `Ref_RC_Data` 77, `GPtFeature` 92,
`TAlnInfo` 129 (merge-only per-record info), `CJunction` 141, `CReadAln` 169,
`GReadAlnData` 232, `CBundle` 260, `CPrediction` 270, `BundleData` 314.

## 3. How it runs

- The main thread is the producer: it reads coordinate-sorted alignments,
  fills one `BundleData` per locus (`evalReadAln`, `processRead`), and hands
  the bundle to a worker. A locus ends when the next alignment starts beyond
  the current end plus the bundle gap, or the reference changes.
- A worker runs `processBundle` -> `infer_transcripts` -> print. Printing and
  the global gene counter `GeneNo` are under `printMutex`.
- Normal mode does not write the final GTF directly. Workers print to
  `<tmpdir>/<name>.tmp`; after all bundles, `main` rereads it to fill in
  FPKM/TPM (they need genome-wide totals) and writes the output. `-A` gene
  abundance (`g_out`) is written in that same final pass. Guides that never
  entered a bundle are written last by `writeUnbundledGuides`.
- Options are plain globals in `stringtie.cpp`, declared `extern` at the top
  of `rlink.cpp`. There is no options struct on master. Adding an option
  means: usage text, the `GArgs` format string, `processOptions`, a global,
  and an `extern` wherever it is read.
- Merge mode (`--merge`) reuses the assembly machinery: each input GTF is
  converted to a temporary BAM in which every transcript is one "read"
  (`TInputFiles::convert2BAM`), references are sorted alphabetically, and
  `TInputFiles::next` does the multiway merge. Record metadata travels in aux
  tags and lands in `TAlnInfo`.

## 4. Conventions that are not obvious from reading one function

- Strand is encoded two ways. Reads, junctions and groups use `-1/0/1`
  (`CReadAln::strand`); arrays indexed by strand use `0` = minus, `1` =
  unstranded, `2` = plus; predictions and GTF use `'-'`, `'.'`, `'+'`.
  Conversions are written inline everywhere.
- `CReadAln::nh == 0` means "this read is deleted or collapsed into another",
  not "zero hits". Code skips such reads with `if(!rd.nh) continue`.
- `readlist` is a sorted `GList`. Adding to it or calling `setSorted` can
  reorder it, so an index saved earlier is stale afterwards. I hit this as a
  real bug (guides are appended and the list is re-sorted in `build_merge`).
  Keep pointers, not indices, across such points.
- `CPrediction::t_eq` is the guide a prediction equals; `geneno` groups
  predictions into genes; `mergename == "n"` marks nascent models, which
  `predCmp` sorts last.
- Guide state lives in the guide object: `getGuideStatus`/`setGuideStatus`
  with `GBST_*` values decide whether `writeUnbundledGuides` still has to
  print it. `uptr` of a guide holds its `RC_TData`.
- Splice-graph code identifies edges through `gpos` (hash of node pair to bit
  position) and paths through `GBitVec` patterns over nodes and edges. Node 0
  is the source and the last node the sink.
- Many functions take 10 to 20 positional parameters. When changing a
  signature, grep every call site; the compiler will not always catch a
  swapped pair of same-typed arguments.

## 5. gclib pitfalls (each one cost me a debugging round)

- `GStr::operator=(int)` exists, so `s = some_char` stores the character code
  as decimal text. Use `s = ""; s += c;`.
- `GStr(const char* s, uint n)`: `n` is extra capacity, not a length. There is
  no substring constructor; build substrings with `substr` or a loop.
- `GStr::appendfmt` formats into a buffer of `strlen(fmt) + 1024` bytes with
  `vsprintf`. Never pass it an unbounded `%s`.
- `GVec<T> v(n, T(x))` is the "n copies" constructor, next to
  `GVec(int, T*, bool)`. Write the cast (`GVec<bool> keep(n, bool(true))`,
  `GVec<const char*> f(n, (const char*)NULL)`) so a literal 0 or NULL does not
  select the pointer overload.
- `GHash<T>(false)` does not own its values. `Find` returns `T*` (NULL when
  absent). Iteration (`startIterate`/`NextData`) is in hash order; keep a
  parallel `GVec<GStr>` of names when output order matters (see
  `writeMergeReport` on merge-dev).
- `GPVec::Sort()` and `GList` ordering use the element's `operator<` and
  `operator==`; both must be defined.
- `GffObj::names` (the shared name-to-id table, including reference sequence
  ids) did not keep the names registered by the merge input readers once those
  readers were gone: lookups at output time returned -1. I did not read the
  gclib code that frees it; treat ids from one `GffReader` as invalid after it
  is destroyed.
- `singleExonTMatch()` in `gclib/gff.cpp` computes `=`/`~` and then ends with
  `return 0;`, so single-exon pairs never match through `getOvlData()` (they
  come out as `c`, `k` or `o`). Present on master. Fixed on merge-dev by
  porting upstream gclib commit 94c50ea (`~/src/gclib`); nothing on master
  calls these functions.
- gclib and htslib are vendored. A gclib fix (PR #482) arguably belongs
  upstream in gpertea/gclib too.
- AGENTS.md on merge-dev asks for C++11, C-style code, and gclib containers
  over STL. That is the house style on master as well.

## 6. Build and test

- Use `make release`. Plain `make` is a DEBUG build (`-DDEBUG -D_DEBUG
  -DGDEBUG`), much slower, and its objects are not link-compatible with
  release objects: mixing them gave `free(): invalid pointer` in
  `GSamRecord::clear`. When switching mode, or switching between master and
  merge-dev (their `bundle.h`/`rlink.h` differ), run `make clean release`.
- Header dependencies in the Makefile are incomplete (`stringtie.o` does not
  depend on `rlink.h` or `bundle.h`). After a header edit, `touch` the `.cpp`
  files or clean.
- Run `ulimit -c 0` before experiments. A crash otherwise writes `core` into
  the repo root.
- `make test` runs `run_tests.sh`: it downloads test data (tag `v3.0.2` of
  gpertea/stringtie-testdata, needs network) into `tests/` and runs 9 assembly
  cases against expected GTFs, ignoring `#` lines. On master a mismatch prints
  `Error: test failed` but the script still exits 0 (the `exit 1` is commented
  out). Read the output; do not trust the return code.
- The test data is tiny: 5 to 13 transcripts per case (47 for the `--nasc`
  case), 1.7 MB of BAM in total. Passing proves the common path still runs,
  nothing more. There is no merge test, no `-A`, `-e`, `-B` or `-p` test on
  master.
- `run_tests_valgrind.sh` exists on master. On merge-dev,
  `VALGRIND=1 ./test_merge_policy.sh` runs the merge fixtures under valgrind.
- The code graph MCP (`codebase-memory-mcp`, project
  `home-gpertea-src-stringtie-git`) indexes the working tree. After a branch
  switch wait for it to refresh, or use `rg`, before trusting line numbers.

## 7. How to check a change (what worked)

- Compare transcript structure sets, not files: reduce each GTF to
  `ref strand exon,exon,...` per transcript, sort, and `cmp`. IDs and order
  can differ for reasons unrelated to the change.
- Run every change at `-p 1` and at `-p 8`. Before merge-dev, generated IDs
  and record order in merge output depended on which worker finished first.
  For assembly mode I expect the same from the code (`GeneNo` is assigned at
  print time under `printMutex`) but could not show it: the bundled test data
  is too small to produce a difference.
- Keep a baseline binary. `/home/gpertea/src/stringtie/stringtie` is a 3.0.3
  release build used as the reference in this session.
- For anything that changes output, a real dataset says more than the unit
  tests. Two are set up (see the memory note `brainseq-merge-eval-dataset`):
  the alga set from issue #247 under `merge-dev/merge-dev/example_dataset/`
  (20 GTFs plus a 26,059-transcript guide) and 30 human brain assemblies under
  `merge-dev/brainseq_eval/`, scored against GENCODE 49 with
  `/opt/sw/bin/gffcompare`. Scripts are in
  `merge-dev/brainseq_eval/survey/scripts/`.
- Write the negative assertion too. The merge suite on the reviewed branch
  passed while asserting almost nothing, because its helper for "transcript
  absent" could not fail.
- Have an independent agent review a non-trivial diff. One found five real
  defects in my own changes, including a regression.

## 7a. What the assemblies themselves look like

`worklogs/assembly-observations-polyA-ribozero.md` holds measurements of
single-sample assemblies (no `-G`) from 10 human samples sequenced as polyA
and as RiboZero, the latter with and without `-N`: per-sample precision, what
`-N` adds and removes, how replication across samples predicts correctness,
how far 5' and 3' ends are from the annotation, fragments and intron overhang.
Read it before changing the read-alignment assembler; it also describes a
ready dataset and scripts for measuring the effect of such a change.

## 8. What I know about master's behavior (checked with the 3.0.3 binary
   unless marked)

Merge mode:

- Abundance filters default to `-F 1 -T 1 -f 0.01`, so low-expressed input
  transcripts are dropped before merging.
- Guides can be altered: on the alga data 13,958 of 26,059 guide IDs came out
  with a changed structure (mostly moved ends).
- An unstranded multi-exon input is written as one single-exon transcript per
  exon.
- Output order: references alphabetical regardless of input order (#177);
  guides without input support appended at the end; 1,337 out-of-order
  transcript lines on the alga data at `-p 1`.
- Generated `MSTRG.N` numbers and record order vary with `-p`.
- A guide named like a generated ID (`MSTRG.1.1`) collides with it (noted in
  the review as present on master; fixed on merge-dev).
- On 10 human samples the default merge kept 20,360 intron chains matching
  GENCODE at 49.5% precision, against 33,021 at 49.3% for merge-dev with
  `--min-support 2`.

General (inferred from the code, not tested on master):

- Numeric options are read with `GStr::asInt`/`asDouble`, which accept
  trailing garbage and turn non-numbers into 0.
- Alignments on `-x` excluded references, and in merge mode input transcripts
  shorter than `-m`, are skipped without a message.

## 9. Working with the maintainer

Standing instructions are in `~/.claude/CLAUDE.md` and memory; the ones that
shaped this session:

- Concise, ASCII-only output. No AI attribution in commits. Commit and push
  only when asked, to the branch named; the maintainer opens PRs themselves.
- Data before opinion. Defaults were accepted or overturned on measurements,
  and every open policy choice had to be a command-line option.
- `pushover-notify -p stringtie-git "..."` at decision points and when long
  work finishes.
- Bulk compute (assemblies, batch conversions) goes to an Opus subagent.
- `worklogs/planning.txt` must not be modified or committed. The review
  documents `worklogs/review-branch-*.md` and this file are untracked on
  purpose; they contain local data paths.
- `~/work/cbrain/` is the maintainer's data; treat it as read-only.

## 10. If a master PR touches merge code

Merge-dev rewrote how `printMergeResults` decides what to emit and how
`build_merge` collects evidence, and it changed `writeUnbundledGuides` and the
merge output path in `stringtie.cpp`. A master fix in those functions will
conflict and may be moot. Check `git show merge-dev:MERGE_MODE.md` and
`git diff master merge-dev -- <file>` first, and say so in the PR discussion
rather than fixing the same thing twice.
