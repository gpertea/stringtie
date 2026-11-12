# GHash Clear optimization investigation and scratch benchmarks

Date: 2026-11-12

Repository: `/home/gpertea/src/stringtie`

Inspected revision: `53d964779c5782e0687ae9907dc0f31f6748d355` (`master`).

Status: source inspection and temporary microbenchmarks only. No repository source,
build settings, git state, or existing experiment files were changed. The user
explicitly requested scratch-space testing and subsequently authorized this worklog.

## Context

The broader issue #511 investigation concerns producer-side CPU bottlenecks before
StringTie hands a bundle to its assembly workers. The same-start duplicate search
in `processRead()` is one bottleneck, but profiling also identified
`GHash<int>::Clear()` as expensive. This side investigation examines the hash
container itself, independently of the same-start candidate-index implementation.

The relevant hash is `hashread`, declared in `stringtie.cpp:769`. It maps an owned
string key containing read name, position, and hit index to a readlist index for
mate pairing. `processRead()` adds and removes entries. At bundle boundaries,
`main()` calls `hashread.Clear()` (`stringtie.cpp:953`). It is not the proposed
same-start structure index.

This is a potentially surprising performance problem: a container that is empty
or nearly empty can remain expensive to clear because it retains the capacity
needed by an earlier large bundle. Repeated small bundles can therefore pay for
an earlier peak allocation repeatedly.

The existing stack-sampling file
`issue_511/work/i511/pmp_br2.txt` contains repeated producer stacks stopped at
`GHashMap.hh:460-461`, in the bucket scan inside `GHash<int>::Clear()`, called from
`main()` at `stringtie.cpp:953`. Worker stacks in the sampled process were waiting
for bundles. These samples identify a concrete hot code path; they do not provide
a statistically established runtime percentage or the hash's live count/capacity
at every clear.

The existing `issue_511/work/i511/perf.log` records a failed hardware-performance
profiling attempt because `perf_event_paranoid` was 4. No permissions or system
configuration were changed in this investigation. The measurements below use
ordinary wall-clock timing inside a standalone program.

Related worklogs:

- `worklogs/astra-i511-long-read-feasibility-testing.md`: input-data survey and
  matching-semantics concerns.
- `worklogs/fable-i511-review-and-test.md`: separate same-start lookup work.

## Conclusions

Two small changes to the non-pointer-value `GHash::Clear()` overload have a strong
case:

1. Reset `lastKey` and immediately return when `count == 0`.
2. For a nonempty hash, scan occupancy bitmap words and enumerate their set bits,
   rather than checking every bucket individually.

In the scratch benchmark, a hash with 1,048,576 retained buckets and only 16 live
keys cleared in about 0.0153 ms with both changes, versus 0.678 ms for the current
code: approximately 44 times faster. Empty clears became approximately timer
overhead. At 400,000 live entries, the improvement was only about 1.13 times.

These are clear-operation microbenchmark results, not whole-StringTie speedups.
They establish that sparse retained capacity is an optimization opportunity.
Dense clears still require a separate free for each owned key, which limits what
bucket traversal improvements can accomplish.

## Current implementation and ownership

`GHash<V>` is defined in `gclib/GHashMap.hh:390`. It derives from
`GHashMap<const char*, V, ...>`, backed by cached-string-key khashl storage.
The occupancy representation is a `uint32_t` array with one bit per bucket.

Relevant source locations:

```text
Location                 Behavior
-----------------------  -------------------------------------------------------
GHashMap.hh:394          GHash constructor and freeItems setting
GHashMap.hh:398          Add: duplicate a newly inserted string key with Gstrdup
GHashMap.hh:431          non-pointer-value Remove: free key and delete bucket
GHashMap.hh:443          pointer-value Clear overload
GHashMap.hh:456          non-pointer-value Clear overload used by GHash<int>
GHashMap.hh:468          Reset: Clear, then release table storage
khashl.hh:42             occupancy-bit helpers
khashl.hh:51             n_buckets: capacity, not live count
khashl.hh:56             underlying clear: bitmap memset and count reset
khashl.hh:140            deletion with backward shifting and occupancy-bit removal
stringtie.cpp:769        hashread declaration
stringtie.cpp:953        hashread.Clear at bundle boundaries
```

The essential current non-pointer-value clear is:

```cpp
khInt_t nb=this->n_buckets();
for (khInt_t i = 0; i != nb; ++i) {
    if (!this->_used(i)) continue;
    GFREE(this->key(i));
}
lastKey=NULL;
this->clear();
```

The inherited lower-case `clear()` does:

```cpp
if (!used) return;
memset(used, 0, __kh_fsize(n_buckets()) * sizeof(uint32_t));
count = 0;
```

Let B denote allocated buckets and N denote live entries. The current wrapper
performs B occupancy tests and N key frees, followed by clearing roughly B/8 bytes
of occupancy storage. Even N=0 still performs the B-iteration loop when capacity
has previously been allocated.

For the benchmark capacity B=1,048,576, the occupancy bitmap is only 131,072 bytes
(128 KiB), or 32,768 bitmap words. The current wrapper nevertheless performs
1,048,576 individual bit tests on each call.

Important ownership distinctions:

- `GHash<int>` owns copies of its string keys. `Add()` duplicates each new key.
- `GHash<int>(false)` does not disable string-key copying or freeing. The flag
  controls deletion of pointer values, when applicable; `int` values need no delete.
- Calling inherited lower-case `clear()` directly on an owning `GHash` bypasses
  freeing the key strings and leaks their allocations.
- `GHashMap<const char*, int>` borrows keys. Substituting it is safe only if another
  owner guarantees stable key storage for the entire time each key is in the table.
  The reused `_id` string in `processRead()` does not meet that requirement.
- Changing the hash function does not address the current clear loop: it does
  not compute hashes or perform key comparisons.

## Optimization 1: empty-hash fast path

The proposed fast path at the start of `GHash::Clear()` is:

```cpp
lastKey=NULL;
if (!this->count) return;
```

Resetting `lastKey` preserves the wrapper's existing postcondition even when no
entries remain. The empty return avoids both the per-bucket scan and the bitmap
`memset`.

This relies on an existing container invariant: a zero live count implies no
occupied bits. khashl deletion uses backward shifting and marks the final vacated
bucket unused; it does not leave tombstones that need cleanup during clear.
Newly allocated occupancy arrays are zeroed, and ordinary clear also zeros them.
The fast path is valid for normal use of these APIs, not for a corrupted table or
external manipulation of internal state.

The scratch `empty_guard` variant prepended these lines to the original method
body. It deliberately retained the original later `lastKey=NULL`, making the
experiment minimal. A production implementation need not keep the redundant
assignment.

An analogous count guard in khashl's lower-case `clear()` could benefit other
wrappers, but would not by itself eliminate `GHash`'s preceding bucket loop. That
broader change was not made or benchmarked here.

## Optimization 2: bitmap-word traversal

The `bitmap` scratch variant replaced the non-pointer-value clear body with:

```cpp
lastKey=NULL;
if (!this->count) return;
khInt_t nb=this->n_buckets();
khInt_t nw=nb<32 ? 1 : nb>>5;
for (khInt_t w=0; w<nw; ++w) {
    uint32_t live=this->used[w];
    while (live) {
        unsigned bit=__builtin_ctz(live);
        khInt_t i=(w<<5)+bit;
        GFREE(this->key(i));
        live &= live-1;
    }
}
this->clear();
```

The outer loop skips empty groups of 32 buckets. For each nonzero word,
`__builtin_ctz` finds the lowest occupied bit, and `live &= live-1` removes that
bit from the local copy. Only occupied keys are visited.

Important details:

- `__builtin_ctz(0)` must never be called; the `while (live)` condition prevents it.
- Tables smaller than 32 buckets still have one bitmap word. Only valid bucket
  bits are set under normal container invariants.
- The index calculation uses `khInt_t` for the word offset, preserving the
  container's bucket-index width.
- The code enumerates words and bits in ascending order, preserving the old
  wrapper's ascending-bucket cleanup order.
- It does not add a list of occupied buckets or modify insertion/deletion behavior.
- The existing lower-level clear still zeros occupancy storage and resets count.
- Only the non-pointer-value overload was changed in the temporary header copies.
  Applying the same traversal to pointer-value cleanup would additionally need to
  preserve conditional value deletion and be tested separately.

Traversal becomes O(B/32 + N), plus N frees and the bitmap reset. The asymptotic
dependence on capacity remains for nonempty hashes, but the costly scalar scan is
much shorter for sparse tables. Empty hashes take the fast path instead.

`__builtin_ctz` worked with the compiler used here. A library patch should follow
the project's compiler-portability requirements, using an existing portable helper
or fallback if required by supported toolchains.

## Scratch environment and benchmark design

Scratch directory:

```text
/tmp/stringtie-ghash-clear-995r_qmx/
    bench.cpp
    check.cpp
    baseline/GHashMap.hh
    baseline/bench
    baseline/check
    empty_guard/GHashMap.hh
    empty_guard/bench
    empty_guard/check
    bitmap/GHashMap.hh
    bitmap/bench
    bitmap/check
```

The baseline header was copied from the repository. The other two copies changed
only the non-pointer-value `GHash::Clear()` overload as described above. Each
temporary include directory preceded the repository gclib include path, so its
`GHashMap.hh` was selected while dependencies remained the same. `GBase.cpp` was
read and compiled from the repository; binaries were written only under `/tmp`.

Environment observed:

```text
Architecture    x86_64
CPU             AMD Ryzen 9 5900X 12-Core Processor, 24 logical CPUs
Compiler        g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
Language        C++11
Benchmark flags -O3, assertions enabled
Timer           std::chrono::steady_clock
```

For each case, a new hash was first filled with 400,000 keys and cleared. This
established 1,048,576 retained buckets. Setup and its initial clear were outside
the timed loop.

The test keys were `read_<integer>_pos_123456_hi_1`, with integer values matching
the key's integer. Each measured round inserted the requested number of keys,
then measured only the call to `Clear()`. The printed `total_s` additionally
included insertion and the rest of the measured-round loop, but excluded initial
capacity establishment.

For the zero-live case, 16 keys were added and individually removed before timing.
Thus the empty table was also exercised after deletion, not only after an earlier
clear. Subsequent timed calls repeatedly cleared the empty retained table.

```text
Live entries  Rounds per process  Retained buckets
------------  ------------------  ----------------
           0                 300           1048576
          16                 300           1048576
       32768                  20           1048576
      400000                   5           1048576
```

There were three repetitions of the complete sequence, each executing the variants
in the order baseline, empty_guard, bitmap. The reported summary is the median
across the three process-level average clear times for each case. It is not the
median of all individual clears.

No CPU pinning, randomized variant ordering, allocator replacement, background-load
control, or frequency control was used. These limitations matter particularly for
small differences and dense cases. The sparse-table improvement was much larger
than the observed run-to-run variation.

## Results

Median milliseconds per `Clear()`:

```text
Live keys   Baseline  Empty guard  Guard + bitmap  Baseline / bitmap
----------  --------  -----------  --------------  -----------------
         0  0.678463     0.000020        0.000020  timer-floor result
        16  0.678386     0.453106        0.015330              44.25x
     32768  1.455458     1.278906        0.681859               2.13x
    400000 14.061292    14.050691       12.467581               1.13x
```

The approximately 20 ns empty results are near clock-measurement overhead, and the
compiler can optimize a known-empty fast path. Do not interpret them as a precise
20 ns cost for an application call or advertise a numerical speedup ratio there.

The empty-guard-only variant also ran faster in some nonempty sparse cases despite
retaining the same bucket loop. That is not an algorithmic reduction for nonempty
hashes. Compiler code generation/layout or other benchmark effects may explain it;
this investigation did not inspect assembly to establish the cause. The bitmap
variant was still about 29.6 times faster than the guard-only variant at 16 entries.

For dense clears, the modest reduction is consistent with per-key deallocation
dominating the remaining work. No separate allocator profile was collected, so the
exact fraction attributable to `free` is not established.

All raw average clear times, in milliseconds:

```text
Repeat  Live keys   Baseline  Empty guard  Guard + bitmap
------  ----------  --------  -----------  --------------
     1           0  0.674887     0.000020        0.000020
     1          16  0.673029     0.451311        0.015131
     1       32768  1.442452     1.278906        0.678116
     1      400000 13.543990    14.182098       12.311341
     2           0  0.678463     0.000020        0.000020
     2          16  0.678386     0.454281        0.015367
     2       32768  1.455458     1.281373        0.684799
     2      400000 14.061292    14.050691       13.009526
     3           0  0.683361     0.000020        0.000021
     3          16  0.686421     0.453106        0.015330
     3       32768  1.484618     1.268603        0.681859
     3      400000 14.563062    13.605347       12.467581
```

Raw measured-loop totals in seconds, including insertion:

```text
Repeat  Live keys  Baseline  Empty guard  Guard + bitmap
------  ---------  --------  -----------  --------------
     1          0    0.2025       0.0000          0.0000
     1         16    0.2022       0.1357          0.0048
     1      32768    0.0802       0.0777          0.0647
     1     400000    0.3346       0.3491          0.3290
     2          0    0.2035       0.0000          0.0000
     2         16    0.2038       0.1365          0.0049
     2      32768    0.0806       0.0770          0.0659
     2     400000    0.3421       0.3363          0.3483
     3          0    0.2050       0.0000          0.0000
     3         16    0.2062       0.1362          0.0049
     3      32768    0.0820       0.0770          0.0658
     3     400000    0.3470       0.3351          0.3418
```

The dense total-loop timings illustrate why a faster isolated clear need not yield
a clear improvement for an entire workload. Insertion and allocator behavior also
contribute.

## Correctness checks performed

The benchmark program included assertions covering repeated insertion, lookup,
removal of alternate keys, clear, lookup after clear, repeated clear, reset, and
reuse. Zero integer values were treated as valid stored values.

A separate program used a custom hash function returning zero for every key to
force collisions. For each population from 0 through 129, it:

1. Inserted distinct string keys and checked successful insertion.
2. Removed every even-numbered key, exercising backward-shift deletion in a cluster.
3. Checked that odd-numbered entries remained accessible with correct values and
   removed entries were absent.
4. Cleared the table and checked zero count, null `lastKey`, and absence of all keys.
5. Cleared again, inserted a key with value zero, and verified that a duplicate add
   did not replace the value.
6. Removed the last key, cleared the empty hash, reset it, and checked zero capacity.

All three variants were compiled and run with:

```text
-std=c++11 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer
```

All completed successfully with no sanitizer diagnostics. They exercise unallocated
tables, capacities below one bitmap word, growth across word boundaries, collisions,
partial deletion, empty-after-deletion, repeated cleanup, and destructor cleanup.

This is not exhaustive validation. Pointer-value ownership, custom value destructors,
all supported compilers, full-range bucket sizes, and whole-StringTie behavior were
not tested. The benchmark did not alter the pointer-value overload.

## Alternatives and tradeoffs

### Release excess capacity at bundle boundaries

An existing scratch source file,
`issue_511/work/i511/src/stringtie_clear.cpp:955`, already contains:

```cpp
if (hashread.n_buckets()>(1<<16)) hashread.Reset();
else hashread.Clear();
```

This was inspected, not edited or rerun. `Reset()` first performs cleanup, then
releases the bucket array and occupancy bitmap. The current large clear still
happens once, but future bundles no longer inherit that capacity.

It is a plausible targeted remedy for an unusually large bundle followed by many
small ones. Consecutive large bundles instead pay repeated allocation and table
growth costs. The 65,536-bucket threshold is an existing experimental choice, not a
threshold established by these microbenchmarks.

Library-level fast paths and bitmap traversal preserve capacity reuse and avoid a
workload-specific threshold. A shrink policy can still be useful for retained
memory, but should be evaluated separately with actual bundle-size distributions.

### Pool owned string keys

If allocation/free overhead dominates dense tables after traversal improvements,
an arena or chunk pool could reduce many individual allocations/frees to bulk
operations. This is a larger design change:

- keys must remain stable while stored;
- `Remove()` semantics and memory reclamation need a deliberate policy;
- frequently removed keys can leave unreclaimed arena storage until clear;
- pointer values, if any, still need their required destruction;
- reset, rehash, and destruction ownership must remain correct.

This was neither implemented nor benchmarked. It should not precede the simple
clear-loop improvements without evidence that allocation remains important.

### Track occupied indices or use generation counters

A separate live-entry list could enable O(N) cleanup, but maintaining indices is
nontrivial because rehash and backward-shift deletion move entries. Generation
counters can avoid some bitmap clearing, but do not eliminate owned-key cleanup
and introduce state/storage costs into lookup and insertion. Neither is needed to
obtain the measured sparse-table gain.

### Preserve ownership rather than bypass cleanup

Replacing `Clear()` with lower-case `clear()`, disabling `freeItems`, or borrowing
temporary key pointers is not a safe shortcut for the current `hashread` usage.
Likewise, changing wyhash to another hash function cannot remove a capacity scan
that does not hash anything.

## Suggested next validation before a repository patch

The smallest justified library patch would add the count guard and bitmap
traversal to the non-pointer-value overload, followed by a separately tested
equivalent cleanup path for pointer values if desired. Avoid unrelated hash changes.

Before claiming application benefit:

- Measure `hashread` capacity, live count, clear-call count, and cumulative clear
  time on the workloads that produced the existing stack samples.
- Compare the baseline, empty guard, bitmap traversal, and optional capacity-reset
  policy using otherwise identical builds and inputs.
- Include one large bundle followed by many small bundles, consecutive large
  bundles, ordinary paired reads, and long-read/mixed-mode inputs.
- Measure total wall time, CPU time, memory, and retained capacity, not only clear
  latency. Check that allocation costs have not just moved elsewhere.
- Confirm identical assembly/count outputs and run the repository's relevant tests.
- Test pointer-value ownership if broadening the change beyond `GHash<int>` and
  other non-pointer values.

Because gclib is vendored, a reusable container improvement should also be considered
for the upstream gclib repository after validation. No upstream changes or publication
were performed here.

## Reproduction

The exact two test programs are included below so the experiment remains reproducible
if `/tmp` is cleared. The temporary header variants can be recreated by copying
the inspected revision's `gclib/GHashMap.hh`, then replacing only its non-pointer
`GHash::Clear()` body as specified above. Leave the separate non-pointer
`GHashMap::Clear()` overload unchanged.

Commands equivalent to the successful compilations and runs, using the existing
scratch directory:

```bash
scratch_dir=/tmp/stringtie-ghash-clear-995r_qmx
repo_dir=/home/gpertea/src/stringtie
for variant in baseline empty_guard bitmap; do
  g++ -std=c++11 -O3 \
    -I"$scratch_dir/$variant" -I"$repo_dir/gclib" \
    "$scratch_dir/bench.cpp" "$repo_dir/gclib/GBase.cpp" \
    -o "$scratch_dir/$variant/bench"
  g++ -std=c++11 -O1 -g -fsanitize=address,undefined \
    -fno-omit-frame-pointer \
    -I"$scratch_dir/$variant" -I"$repo_dir/gclib" \
    "$scratch_dir/check.cpp" "$repo_dir/gclib/GBase.cpp" \
    -o "$scratch_dir/$variant/check"
  "$scratch_dir/$variant/check"
done
for repeat in 1 2 3; do
  for variant in baseline empty_guard bitmap; do
    "$scratch_dir/$variant/bench"
  done
done
```

The original execution printed labels identifying repeat and variant before each
program run. Capture those labels when collecting new results to avoid ambiguous
output. Do not compile against a later changed header and describe it as the same
baseline without recording that revision.

The only pre-existing tracked modifications visible before and after testing were
`.gitignore` and `AGENTS.md`; neither was edited here. `worklogs` resolves to
`../stringtie-git/worklogs`, so this requested document lives through that existing
symlink. No source optimization has been applied.

## Benchmark program

```cpp
#include "GHashMap.hh"
#include <chrono>
#include <cstdio>
#include <cassert>
using Clock=std::chrono::steady_clock;
static void add(GHash<int>& h, int n) {
 char key[48];
 for(int i=0;i<n;i++) { snprintf(key,sizeof(key),"read_%d_pos_123456_hi_1",i); h.Add(key,i); }
}
static void check() {
 GHash<int> h;
 for(int pass=0;pass<20;pass++) {
  int n=(pass+1)*37;add(h,n);
  char key[48];
  for(int i=0;i<n;i++) { snprintf(key,sizeof(key),"read_%d_pos_123456_hi_1",i); assert(h.Find(key) && *h.Find(key)==i); }
  for(int i=0;i<n;i+=2) { snprintf(key,sizeof(key),"read_%d_pos_123456_hi_1",i); assert(h.Remove(key)>=0); }
  h.Clear(); assert(h.Count()==0 && h.getLastKey()==NULL);
  for(int i=0;i<n;i++) { snprintf(key,sizeof(key),"read_%d_pos_123456_hi_1",i); assert(h.Find(key)==NULL); }
  h.Clear();
 }
 h.Reset(); assert(h.n_buckets()==0);add(h,9);h.Clear();
}
int main() {
 check();
 const int live[]={0,16,32768,400000};
 const int rounds[]={300,300,20,5};
 for(int c=0;c<4;c++) {
  GHash<int> h;add(h,400000);h.Clear();
  // exercise empty after removing all entries, retaining capacity
  if(c==0) { add(h,16);char key[48];for(int i=0;i<16;i++){snprintf(key,sizeof(key),"read_%d_pos_123456_hi_1",i);h.Remove(key);} }
  double seconds=0;auto total_start=Clock::now();
  for(int r=0;r<rounds[c];r++) {
   add(h,live[c]);assert(h.Count()==(unsigned)live[c]);
   auto a=Clock::now();h.Clear();auto b=Clock::now();
   seconds+=std::chrono::duration<double>(b-a).count();assert(h.Count()==0);
  }
  printf("live=%d buckets=%llu rounds=%d clear_ms=%.6f total_s=%.4f\n",live[c],(unsigned long long)h.n_buckets(),rounds[c],seconds*1000/rounds[c],std::chrono::duration<double>(Clock::now()-total_start).count());
 }
}
```

## Sanitizer correctness program

```cpp
#include "GHashMap.hh"
#include <cassert>
#include <cstdio>
struct Collide { uint32_t operator()(const char*) const { return 0; } };
int main() {
 for(int n=0;n<=129;n++) {
  GHash<int,Collide> h;
  char key[32];
  for(int i=0;i<n;i++){snprintf(key,sizeof(key),"key_%d",i);assert(h.Add(key,i)>=0);}
  for(int i=0;i<n;i+=2){snprintf(key,sizeof(key),"key_%d",i);assert(h.Remove(key)>=0);}
  for(int i=0;i<n;i++){snprintf(key,sizeof(key),"key_%d",i);int* v=h.Find(key);assert(i%2 ? v && *v==i : !v);}
  h.Clear();assert(h.Count()==0 && h.getLastKey()==NULL);
  for(int i=0;i<n;i++){snprintf(key,sizeof(key),"key_%d",i);assert(!h.Find(key));}
  h.Clear();assert(h.Add("again",0)>=0);assert(h.Add("again",1)==-1);assert(*h.Find("again")==0);
  h.Remove("again");assert(h.Count()==0);h.Clear();h.Reset();assert(h.n_buckets()==0);
 }
 puts("collision, small-capacity, removal, repeated clear/reset checks passed");
}
```
