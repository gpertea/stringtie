# Critical assessment of the merge branch review

Review date: 2026-11-09. Assessed `merge-dev` at `a27c05d` against `master`
at `d1dc38d`, in response to
[review-branch-to-fix-before-merge.md](review-branch-to-fix-before-merge.md).
This assessment does not approve or implement the proposed fixes.

The merge should be held, but the review's proposed fix plan should not be
adopted wholesale. Most reported failures reproduce; several conclusions and
acceptance criteria overreach.

The assessment covered the README, merge architecture, planning notes,
implementation diff, and tests. The codebase MCP index was checked for
freshness and coverage. Fresh release builds of the application code were
made from both revisions, linked against the existing vendored dependency
libraries. Review experiments used scratch files and did not change tracked
project files.

## Confirmed blockers and additional failures

1. **Endpoint controls miss terminal-exon extensions (F1).** Two files ending
   at 4300 plus one ending at 5300 produce 5300 even with
   `--max-extension-3p 100`. The left/right asymmetry and synthesized endpoint
   pair also reproduce, including minus-strand examples. The check uses
   graph-node bounds, while terminal variants can share nodes.
   Source: `rlink.cpp:7592`, `merge_can_join()`.

2. **Opposite-strand SETs are incorrectly discarded (F5).** The reported
   antisense fixture survives master and disappears on this branch. The final
   containment check lacks a strand condition. This is a clear regression.
   Source: `rlink.cpp:20294`, final SET containment in `printMergeResults()`.

3. **Unstranded loss is broader than F6 describes.** The unstranded MET
   disappears, but contrary to the review, it produces two `minimum_support`
   report rows with support 0. Legacy code splits it into exon predictions;
   exact-structure support recovery cannot associate those predictions with
   the original MET. Two overlapping unstranded SETs, 1000-1300 and 1200-1500,
   also disappear entirely; master emits their union. Fixing only unstranded
   METs would miss this second regression. Sources: `rlink.cpp:16765`,
   unstranded bundle output, and `rlink.cpp:20046`,
   `merge_fill_prediction_support()`.

4. **Guide controls have another bypass beyond F3.** Guide exons
   `1000-1199,1300-1499`, input exons `900-1199,1300-1399`: the variant survives
   both `--guide-policy reference` and `--max-extension-5p 0`, without a
   reference relationship. `merge_related_guide()` requires one complete span
   to contain the other, excluding variants that extend one end and shorten
   the other. Source: `rlink.cpp:20068`.

5. **Direct mode can destroy valid observed paths.** With endpoint thresholds
   relaxed to 1, inputs `1000-1100,2000-2100,3000-3300` and
   `2000-2100,3000-3300,4000-4300` form an unobserved four-exon path.
   `--path-mode direct` rejects that path and emits neither original, although
   each original meets minimum support 1. The mode is enforced after
   destructive coassembly. This needs a fallback or an earlier constraint,
   plus an explicit regression test. Sources: `rlink.cpp:7792`, path
   absorption, and `rlink.cpp:20220`, final direct-path filtering.

6. **The tests provide substantially less protection than their success
   message suggests (F7).** All nine `! has_tx` assertions are ineffective
   under `set -e`. Making one deliberately false still lets the suite print
   `OK`. A corrected scratch copy fails that mutation and passes the real
   assertions. Also, `make test` excludes merge tests, and assembly output
   mismatches do not fail `run_tests.sh`. Sources: `test_merge_policy.sh:110`
   and the other eight negated assertions, `run_tests.sh:75`, and
   `Makefile:181`.

The guide-swallowing examples A/D/E and support-loss example Q also reproduce,
including empty reports. Numeric validation is weak beyond the examples
listed: `--min-support 2x`, `--set-overlap nan`, and `--slop-frac-3p nan` all
exit successfully.

## Qualifications needed before adopting the proposed plan

- **F2 describes an unsuitable policy, not a violation of the documented
  exact-structure policy.** The shifted-end fixture fails exactly as reported,
  but README explicitly defines `--min-support` using exact structures.
  Intron-chain support is a reasonable redesign; it needs a deliberate
  contract change. Calling support above 1 universally "unusable" is too
  strong.

- **F4's observation is solid; its causal claim is not established.**
  Replicated shorter models disappear, and support is lost during cleanup.
  But aggregate retention does not prove weak long models are the principal
  cause. Missing chains can also be redundant fragments. Monotone retention
  across support bins is not a correctness invariant.

- **F1's zero-outlier acceptance criterion conflicts with its own caveat.**
  The 12,851 count excludes endpoint evidence from other compatible chains.
  A fix allowed to use that evidence need not drive this same-chain-only
  statistic to zero.

- **Do not change the default to `auto` yet.** One guided-assembly dataset
  cannot establish precision, and adding an unrelated fourth file would
  change filtering everywhere. Keep default 1 while fixing accounting and
  endpoint behavior; evaluate threshold changes separately.

- **Per-chain absorption rows cannot reconcile every input without provenance
  or membership counts.** The proposed report format and the 927,737-input
  accounting target need to agree.

## What the verification establishes

- All nine assembly outputs match expected results.
- The corrected merge suite passes, including under Valgrind.
- Dataset totals reproduce: 192,173 guided outputs; 26,059 unchanged guide
  IDs; zero duplicate IDs on this dataset.
- Independent recounts reproduce the input chain-support distribution,
  branch retention figures, and 12,851/27,208 endpoint statistic.
- Structures agree at 1, 4, and 8 threads; reports are byte-identical.

These are useful results, not broad correctness proof. The merge suite lacks
direct-mode testing and support-value assertions; most structural checks
inspect outer bounds rather than complete exon chains. Zero duplicate IDs is
also dataset-specific: a guide named `MSTRG.1.1` can collide with a generated
ID, although that defect also exists on master.

## Recommended order of work

First repair test failure propagation and add exact-structure fixtures for
the reproduced failures. Then fix endpoint/guide enforcement, strand
handling, unstranded support, and destructive direct-mode filtering. Support
semantics and new defaults should remain separate decisions.

## Reproduction artifacts

The experiments are retained in the local scratch directory
`/tmp/stringtie-critical-review-m07al9tr`. These are temporary artifacts, not
committed test coverage or a durable benchmark package.

```text
Artifact                  Contents
------------------------  -------------------------------------------------
fixtures.log              Reported fixture results against branch and master
counts.txt                Independent input and output chain counts
comparison.tsv            Full dataset comparison using repository script
valgrind.log              Corrected merge suite Valgrind result
strict-tests.sh           Merge suite with effective negative assertions
mutated-tests.sh          Original suite with a deliberately false assertion
strict-mutated.sh         Corrected suite rejecting the same mutation
shiftguide/left-out.gtf   Same-chain variant bypassing guide policy and cap
direct/                   Shared, direct, and replicated mode fixture output
ids/                      Guide/generated ID collision on both revisions
assembly/                 Nine assembly outputs checked against expectations
guided1.gtf              Dataset output at one thread
guided4.gtf              Dataset output at four threads
guided8.gtf              Dataset output at eight threads
guided1.tsv              Report identical to guided4.tsv and guided8.tsv
```

Reproduction scripts also remain at:

```text
/tmp/stringtie-critical-fixtures.py
/tmp/stringtie-critical-validation.py
/tmp/stringtie-critical-dataset.py
/tmp/stringtie-critical-counts.py
/tmp/stringtie-critical-ends.py
```

Source line numbers in this assessment refer to `a27c05d`.
