# PR 513 code review and integration considerations

Date: 2026-11-12

Scope: code-only review of the issue #511 patch in PR #513 and comparison with the
`i511-long-reads` branch. This document also considers accepting
the PR first, improving its performance afterward, and coordinating those
changes with `lr-junction-deletions`.

No builds, benchmarks, dataset modifications, or GitHub actions were performed
for this review. Results quoted from the PR or existing worklogs belong to their
authors; they were not independently reproduced here. Writing this document does
not authorize any merge or implementation work.

## Revisions and sources

```text
Code                         Revision
---------------------------  ----------------------------------------
master / PR base              2c3fd9cfaeb07850f67e92fc3cae5a2b289c377d
PR #513                       b022e913a595ee851a39a27da2e5c834f97872b5
i511-long-reads                fec648b9aa48a1a5f003defaf3927c8cced4708f
lr-junction-deletions          85b71ae6c6af57ed1a511783aab0f31d00050f35
```

- [PR #513: Fix slow read collection at deep same-start pile-ups](https://github.com/gpertea/stringtie/pull/513)
- [PR #513 patch](https://github.com/gpertea/stringtie/commit/b022e913a595ee851a39a27da2e5c834f97872b5)
- [`i511-long-reads` at the reviewed revision](https://github.com/gpertea/stringtie/tree/fec648b9aa48a1a5f003defaf3927c8cced4708f)
- [Deletion-fix branch at the reviewed revision](https://github.com/gpertea/stringtie/tree/85b71ae6c6af57ed1a511783aab0f31d00050f35)

Local context reviewed:

- `astra-i511-long-read-feasibility-testing.md`
- `fable-i511-review-and-test.md`, including its later corrections and follow-ups
- `fable-junction-deletions-issue-draft.md`
- `fable-i511-followup-leads-and-alignment-flag.md`
- Relevant GHash analysis and its assessment in `astra-GHash-Clear-optimizations.md`

The code graph was checked for freshness and used to locate and trace PR helper
functions in the existing PR worktree. Exact branch differences and source were
read through Git. The graph's current master view was not treated as evidence
for code existing only on another branch.

## Overall assessment

`i511-long-reads` is the stronger performance-only fix against master:
it preserves the existing search window and representative choice, avoids
indexing short runs, and distinguishes important nonmatching states in its keys.
The PR's implementation is substantially smaller, uses existing containers, and
retains the original exact matching predicate. It also intentionally expands
collapsing to the whole bundle and leaves some relevant searches quadratic.

Accepting PR #513 first is a reasonable integration route if its intentional
behavior change is accepted on the available evidence. That choice establishes
a different baseline: subsequent performance improvements should preserve its
whole-bundle, newest-first matching semantics, not silently restore master's
trailing-run restriction.

The distinction matters more than the choice of hash function. The decisions are:

1. Which stored reads are eligible to match?
2. Which fields determine whether an eligible read matches?
3. How can candidates be found without scanning large incompatible groups?

PR #513 changes the first and accelerates part of the third. `i511-long-reads`
preserves the first and accelerates the third using the full current predicate.
The deletion branch changes the second and downstream junction accounting.

## Comparison of the two implementations

```text
Property              PR #513                         i511-long-reads
--------------------  -----------------------------  -----------------------------
Eligible records      Entire current bundle          Trailing same-start run
Candidate order       Newest first                   Newest first
Index key             Exon-coordinate fingerprint    Incoming/stored matching fields
Exact verification    Original predicate             Original predicate
Index activation      Every non-merge input          Runs of at least 64 entries
Index ownership       Per BundleData                 One producer-owned static index
Index lifetime        Through bundle processing      Reused between runs/bundles
Output contract       Additional collapsing allowed  Preserve master's choices
GHash cleanup fixes   Not included                   Included in branch tip
Diagnostic checker    Not added by this patch        RUNIDX_CHECK against old scan
```

The PR's changes are in `bundle.h` and `rlink.cpp`, with 42 additions and seven
deletions. `i511-long-reads` also includes the independent GHash changes;
its total diff must not be interpreted as the cost or complexity of the read
index alone.

## 1. The exon-only index retains quadratic cases

The PR's [lookup loop](https://github.com/gpertea/stringtie/blob/b022e913a595ee851a39a27da2e5c834f97872b5/rlink.cpp#L783-L791)
retrieves the newest stored read with a matching exon fingerprint, then walks
the linked list of earlier reads with that fingerprint. Every candidate still
passes through `samereadmatch()`.

That is effective when exon structures differ. It is not sufficient when many
reads have the same exons but each incoming read fails to match its stored
counterparts. Every failed arrival appends another entry to the same chain.
For N such reads, the candidate checks total N*(N-1)/2, without requiring any
hash collision.

### Cases established by the source

- **Deletion immediately after an intron.** The existing parser can produce two
  `juncsdel` entries for one intron. Incoming and stored junction counts differ,
  so `deljuncmatch()` rejects all candidates. The PR acknowledges this case.
- **Deletion at a guide-matched junction.** Stored junction creation suppresses
  deletion offsets when the guide matches. The current predicate compares the
  resulting zero offsets with nonzero incoming offsets. Repeated identical
  alignments can therefore keep accumulating in one exon chain.
- **Unknown incoming strand, inferred stored strand.** A long read without
  strand information can be stored with a strand inferred from polyA/T evidence.
  The next equivalent incoming read still has strand 0, whereas the stored
  copy has strand 1 or -1. The strand check rejects it before exon comparison,
  but the PR still visits every earlier same-exon entry.

The third case can be constructed with multi-exon reads and an unaligned polyA
tail, avoiding the terminal aligned-exon artifact screen. The worklogs separately
report reproducers confirming that such reads remain separate on master.

Repairing the parser's doubled-entry bug alone does not remove the guide-offset
or inferred-strand cases. The first two are addressed by the deletion branch;
the inferred-strand behavior is not resolved by that branch.

The `i511-long-reads` index derives stored keys from the stored state and incoming
keys from the incoming state. It includes strand and junction count/offsets.
Consequently, these nonmatching states normally have different hashes and avoid
exact comparisons against the entire same-exon chain. Like any hash structure,
it still depends on hash distribution; this is not a deterministic worst-case
constant-time guarantee.

For a later optimization of the PR's code, the useful idea to retain is this
distinction between incoming and stored effective keys. The old trailing-run
restriction is not necessary to obtain that benefit.

## 2. Additional collapsing affects more than sort order

### What PR #513 deliberately changes

Master searches backward only through the consecutive entries at the end of
`readlist` that share the incoming start. Terminal-exon trimming can move a
stored read's start forward and interrupt that run, even though input alignments
were coordinate-sorted before trimming.

The PR's bundle-wide lookup finds older matching entries across those interruptions.
This is intentional, documented behavior. The PR reports unchanged transcript
coordinates in its tested long-read runs, with coverage differences in five runs
and a largest reported change of 8%. Those observations are useful evidence,
but do not establish a general guarantee that only coverage can change.

### Sorting is one possible mechanism

Long-read inference calls the ordinary `readlist.Sort()` without the explicit
tie-breaker used by the mixed-mode path. Changing the number of retained
representatives can therefore affect the ordering of tied records. The PR's
explanation of coverage differences through that ordering is plausible.

This review did not trace the reported differences to a particular downstream
calculation. In particular, it does not establish that sort order caused all
the measured differences, or that making the sort stable would remove them.

### Junction support can change directly

There is another concrete interaction in the existing
[deletion-rescue path](https://github.com/gpertea/stringtie/blob/b022e913a595ee851a39a27da2e5c834f97872b5/rlink.cpp#L17639-L17660).
For each stored read being rescued to an existing junction, the code adds the
shared tentative junction's entire `nreads`, `mm`, and `nm` totals to the target.
It does not transfer only that stored read's contribution.

Consider an illustrative case in which two stored representatives share a
tentative junction with total support C and both enter this transfer branch:

```text
Representation                  Existing transfer to the target
------------------------------  --------------------------------
Two stored representatives      C + C
One combined representative     C
```

The input support remains C in both representations. The difference arises
because the target is updated once per stored representative. If broader
collapsing combines those representatives, it can change target support even
with identical input counts and a fixed iteration order. This example assumes
both representatives reach the existing-target transfer branch; it is not a
newly executed reproducer or a claim about every collapse.

Junction support influences splice-site decisions. Therefore a representation
change can potentially affect transcript structures as well as abundance. The
worklogs' separate count-transfer reproducers demonstrate that this existing
accounting defect can affect splice-site selection; they do not by themselves
demonstrate that PR #513 triggered that outcome on the PR's datasets.

### Implications for accepting the PR

This is an existing defect, not a new incorrect addition inside the PR's hash loop.
Broader collapsing changes how often it is exercised and may reduce some
overcounting. It would be equally unjustified to call every resulting difference
a regression or to assume every difference is an improvement.

The `lr-junction-deletions` branch changes the two deletion-rescue transfers to
add `rd.read_count` and a proportional share of `mm` and `nm`. It does not fix
the analogous shared-total transfer in the later unstranded-junction path.
Nor does correcting these sites prove that all remaining inference operations
are invariant to how equivalent reads are grouped.

If PR #513 lands first:

- Treat the changed collapse behavior as part of the accepted baseline.
- Keep its output differences distinct from the additional correctness changes
  introduced by the deletion branch.
- Use stable sorting, if pursued, as a separately evaluated behavior change;
  do not present it as a complete explanation or remedy for all collapse effects.
- Preserve targeted cases that distinguish candidate selection, count transfer,
  and final inference. A single final-GTF comparison cannot localize a mismatch.

## 3. Coordination with lr-junction-deletions

### What the deletion branch changes

The reviewed branch contains four related corrections:

1. Preserve the direction of deletion offsets instead of constructing them with
   `GSeg(start,end)`, which can swap the values.
2. Add one deletion-offset entry per intron, rather than adding one after `D`
   and another after the following aligned bases.
3. Make guide-matched incoming offsets consistent with the unadjusted junction
   that `processRead()` would store.
4. Transfer per-read support in the two deletion-rescue count-update branches.

The [updated predicate](https://github.com/gpertea/stringtie/blob/85b71ae6c6af57ed1a511783aab0f31d00050f35/rlink.cpp#L239-L259)
takes `GReadAlnData&`, not merely `GVec<GSeg>&`. It needs both the incoming
alignment's deletion offsets and `alndata.juncs[i]->guide_match`.

### Integrating with the PR's present exon-only index

The exon fingerprint remains a valid candidate filter: matching still requires
identical exon coordinates. Pass `GReadAlnData` through `samereadmatch()` and
its call site to the corrected `deljuncmatch()`.

There is no need to add guide status to the PR's current exon-only key to make that
specific predicate change work. The parser corrections alter the incoming
offsets and counts used by the predicate, while the count-transfer correction
operates downstream.

The patches touch overlapping matching code, so this is a semantic integration
step even if a textual merge appears straightforward. Confirm that all helper
calls use the new context and that the final predicate matches the intended
stored-junction construction.

### Integrating with a richer key

The present `keyAln()` of `i511-long-reads` hashes the incoming raw deletion offsets.
That is appropriate for the current master predicate. It becomes insufficient
when the guide-aware predicate is adopted.

For example:

```text
Quantity                                       Value
---------------------------------------------  ------------------------------
Incoming raw offset at a guide-matched intron   Nonzero
Stored effective offset for that junction      Zero
Corrected predicate's expected incoming offset Zero
```

If the predicate is updated but the key remains raw, the incoming and stored
hashes differ. A valid candidate is discarded by the hash filter before the
corrected predicate can accept it. Retaining an exact comparison after lookup
does not protect against this kind of false negative.

A richer index must consequently use:

- Incoming strand and read-type/unitig state as consumed by the predicate.
- Incoming exon coordinates and the predicate's junction-count constraint.
- Zero incoming offsets for guide-matched junctions, otherwise the incoming
  deletion offsets, when using the corrected guide-aware predicate.
- Stored strand/state and offsets derived from the stored exon/junction
  coordinates, as the predicate derives them.

Key construction and exact comparison must agree on which pairs can match.
The essential property is that any pair accepted by the predicate has equal
keys; exact verification then handles fingerprint collisions. Stored keys
must be computed after the fields used by them have reached their final
producer-side state, including inferred strand and unitig status.

Parser corrections also alter offset counts and values. Test the integrated
parser, key, and predicate together rather than relying only on separate
branch results.

### What a merge-first route changes about validation

The current `RUNIDX_CHECK` compares against master's trailing-run scan. That
oracle is intentionally wrong for validating a whole-bundle implementation
after the PR's behavior is accepted: it will reject legitimate new matches.

A follow-up optimization needs a slow reference lookup that visits the entire
current bundle newest first and uses the same predicate as the build under
test. Compare match/no-match and the chosen representative index before
applying accumulation once. For builds with the deletion fixes, both lookups
must use the corrected predicate.

Compare against master separately to document deliberate behavioral changes.
Do not demand master equivalence from a patch stack whose accepted purpose
includes changing collapsing or repairing junction accounting.

## Other implementation observations

### Details that look sound in PR #513

- `samereadmatch()` preserves the original matching checks. Exon-fingerprint
  collisions can add work but cannot silently merge unequal structures.
- Candidate links are newest first, matching the intended representative order
  within the enlarged eligibility set.
- The pointer returned by `last_read_with_exons.Find()` remains valid until its
  use: intervening junction insertion affects a different container. The code
  does not keep using that pointer after an insertion into its own map.
- `prev_read_same_exons` receives one entry for each newly appended non-merge
  read. Producer-side append order keeps these indices aligned with `readlist`.
- Workers may later reorder reads, but the producer no longer uses this bundle's
  index. Cleanup clears both index containers before bundle reuse.
- Merge mode bypasses the new lookup and insertion.

No new memory-safety or index-lifetime defect was identified in these paths by
this code review.

### Costs and maintenance tradeoffs

PR #513 hashes every accepted non-merge input and retains an index per bundle until
cleanup, including while workers process completed bundles. The `i511-long-reads`
index activates only for long runs and reuses one producer-owned allocation.
These are source-level differences; their actual time and memory effects belong
to the ongoing benchmark comparison.

The run tracking of `i511-long-reads` appears consistent with current early-return
paths and bundle reuse. New entries are indexed at the end of `processRead`,
after relevant state is finalized. Its custom storage and duplicated key logic
require more maintenance than the PR's small exon-only index, particularly when
the predicate changes.

The reviewed branch includes GHash bitmap traversal, empty-clear handling,
peak-based shrinking, an iterator correction, and the khashl allocation-check
fix. These are independent of read-candidate eligibility and can be ported
without adopting the trailing-run algorithm. Avoid merging the entire branch
merely to obtain those container improvements if whole-bundle semantics have
already been accepted.

Neither reviewed patch adds a committed targeted regression suite for its new
lookup behavior. The existing worklog experiments and `RUNIDX_CHECK` provide
valuable evidence, but selected small reproducers should eventually become
durable tests alongside the chosen implementation.

## Reading the new worklog evidence

The new notes resolve several uncertainties from the original feasibility review:

- The producer scan was instrumented, and a synthetic workload demonstrated its
  quadratic growth. The local issue #476 bundle exercises it, but is much less
  extreme than the reported issue #511 workload.
- The `RUNIDX_CHECK` per-read checker was exercised with multiple thresholds and
  an intentional negative control. The worklog reports identical single-thread
  GTF output on its tested cases.
- Deletion-offset and count-transfer defects were separated into their own
  branch and reproducers, rather than hidden inside the performance change.
- The GHash cleanup improvement was merged into `i511-long-reads`, so branch-tip
  timing differences include that independent change.
- Later notes correct earlier timing claims and identify variability and code
  placement effects. Use those corrections rather than the superseded single-run
  conclusions. The function-alignment flag remains a separate consideration,
  not a prerequisite established by this code review.

These are recorded experimental findings, not new measurements from this review.

## Possible integration route if PR #513 is accepted first

This route is under consideration; no merge decision is recorded here.

1. **Accept and record the baseline.** Treat PR #513's whole-bundle, newest-first
   collapsing as intentional. Preserve the exact revision and the available
   output comparison evidence, including known drift from master.
2. **Bring over independent container improvements.** Apply the GHash and khashl
   commits separately, without replacing the PR's lookup with the trailing-run
   implementation. This keeps their effect independently reviewable.
3. **Integrate deletion correctness changes.** Update the matching helper to
   pass guide context and validate the parser/predicate/count-transfer changes
   as correctness fixes. Keep the unresolved unstranded transfer site explicit.
4. **Improve candidate selection.** Preserve whole-bundle eligibility while
   indexing effective matching fields. This can draw on the separate incoming
   and stored key construction from `i511-long-reads`, normalized for the
   corrected predicate. Keep exact checks and newest-first selection.
5. **Evaluate common-case overhead and index lifetime.** Thresholds, a small-input
   fast path, or releasing producer-only index state at handoff may be useful,
   but are later implementation choices. A fast path must not accidentally
   restore the trailing-run restriction or forget earlier eligible records.

If performance work precedes the deletion fixes instead, derive keys from the
current predicate and change the incoming key and predicate together when the
guide-aware fix lands. The semantic dependency is mandatory; the order of
independent commits is not.

This route preserves the PR's small initial contribution while allowing broader
optimization. It also means the eventual performance patch should be judged
against the accepted PR-plus-correctness baseline, not solely against master.

## Focused validation considerations for later work

The other agent is already benchmarking. These are code-review-derived cases
to retain or add when appropriate, not tests run for this document:

- A read matching an earlier representative across a trimming-induced run break.
- Repeated same-exon reads that always fail because of inferred stored strand.
- Deletions before and after introns, with and without guide matches, including
  several introns so offset-entry displacement is observable.
- Multiple representatives sharing a tentative junction, with a preexisting
  stranded target, to distinguish per-read transfer from repeated shared totals.
- Equal exon structures with different read type, unitig state, or deletion state.
- Forced fingerprint collisions, ensuring exact checks still reject unequal reads.
- Bundle reuse and worker processing at multiple thread counts.
- A whole-bundle diagnostic oracle using the same predicate as each patch stack.

Keep intermediate representative-selection checks, junction-evidence checks,
and final transcript/abundance comparisons distinct. This makes it possible to
attribute changes to lookup semantics, corrected evidence handling, or inference
instead of treating every final-output difference as a hash-index failure.

## Conclusion

There is no code-review basis here to reject PR #513 as an unsafe hash-table
implementation. The substantive concerns are incomplete elimination of quadratic
nonmatching chains and the intentional enlargement of collapse eligibility.

Accepting it first is compatible with a better final implementation. The critical
conditions are to preserve the accepted whole-bundle behavior in later speedups,
keep effective keys synchronized with the deletion-aware predicate, and evaluate
count-transfer corrections separately from lookup acceleration.
