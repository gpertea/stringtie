# Repository agent instructions

## Code map

```text
File             Role
---------------  -----------------------------------------------------------
stringtie.cpp    CLI, option handling, guide loading, input scan, bundle queue
bundle.{h,cpp}   Per-locus reads, guides, junctions, coverage, and bundle state
rlink.{h,cpp}    Splice graphs, transfrags, transcript inference, result output
tmerge.{h,cpp}   Multi-input merge stream; converts merge-mode GTF to temp BAM
tablemaker.*     Optional Ballgown/raw-count collection and table output
gclib/           Preferred strings, containers, genomics types, I/O, threads
htslib/          Vendored SAM/BAM/CRAM and compression implementation
```

The primary control flow is:

```text
main
  -> processOptions
  -> load optional guide GTF/GFF and point features
  -> TInputFiles/GSamReader coordinate-sorted alignment stream
  -> BundleData::evalReadAln + processRead
  -> close each genomic bundle and queue it to a worker
  -> processBundle
       -> infer_transcripts
            -> normal mode: count_good_junctions -> build_graphs
            -> merge mode: build_merge
            -> group reads -> build splice graphs/transfrags
            -> find transcript paths and estimate abundance
       -> printResults or printMergeResults
       -> optional Ballgown/raw-count updates
  -> finalize unbundled guides, abundance files, and temporary inputs
```

`stringtie.cpp` is the producer and each worker is a bundle consumer. A bundle
is a connected genomic locus expanded by overlapping alignments and guides.
`BundleData` owns its reads, junctions, guide references, coverage, and predicted
transcripts until `processBundle()` prints the result and clears the bundle.
See [MERGE_MODE.md](MERGE_MODE.md) for the detailed `--merge` workflow.

## C++ generation conventions

Use C and C++11. Prefer plain C-style code and the existing gclib types over STL:

```text
Need                 Prefer
-------------------  --------------------------------------------------------
simple text          const char*, char*, and C string functions
mutable/complex text GStr
value sequence       GVec<T> or GArray<T>
pointer sequence     GPVec<T> or GList<T>
lookup               GHash<T>, GHashMap<K,V>, or another existing gclib hash
genomic interval     GSeg
alignment record     GSamRecord, read through GSamReader
transcript model     GffObj
bit set              GBitVec
```

- Search gclib and nearby project code before introducing an STL type or a new
  abstraction. Use STL only when gclib has no practical equivalent or the local
  subsystem already requires it.
- Use indexed `for`/`while` loops, `Count()`, and `operator[]` for `GVec` and
  related containers; these containers do not provide STL iterators. Verify the
  concrete class API rather than assuming full `std::vector` compatibility.
- Reuse interval and genomics operations such as `GSeg::len()`/`overlap()`,
  `GSamRecord::refName()`/strand/exons, and `GffObj` transcript/exon helpers.
- Keep new data structures minimal. Simple structs with public fields and small
  methods are preferred; do not add encapsulation or generic machinery without
  a concrete need.
- Match surrounding ownership, sorting, comparator, and error-handling idioms.
  Do not mechanically mix STL ownership or iterator conventions into gclib
  containers.

## Merge glossary

```text
Term       Meaning
---------  --------------------------------------------------------------
tx         Transcript
SET        Single-exon transcript or transfrag
MET        Multi-exon transcript or transfrag
bridge     Path that connects otherwise distinct tx structures or loci
fusion tx  Transcript joining models assigned to distinct genes
```

## Code graph MCP

Use `codebase-memory-mcp` as the first source-discovery tool. The project name is
`home-gpertea-src-stringtie-git`.

- Confirm the index is current before relying on graph results; refresh the
  repository in `full` mode when it is absent or stale.
- Use `search_graph` for symbols and concepts, `trace_path` for callers,
  callees, impact, and data flow, and `get_architecture` for subsystem boundaries.
- Call `search_graph` before `get_code_snippet` and pass the exact returned
  `qualified_name`.
- Use `search_code` for graph-aware literal or regex searches. Check truncation
  fields and narrow or paginate broad searches.
- Before changing a symbol, inspect inbound and outbound relationships. Refresh
  the index after material source changes before making graph-based claims.
- Use `rg` and direct reads for non-code text, exact verification, unsupported
  files, or when the graph is unavailable. State when falling back.
