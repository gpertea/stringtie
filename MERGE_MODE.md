# StringTie merge mode

`stringtie --merge` combines transcript annotations from multiple GTF/GFF
inputs into one unified, non-redundant isoform set. It operates on previously
assembled transcripts, not raw RNA-seq reads.

## Workflow

```text
Input GTF/GFF files
  -> convert transcripts to temporary alignment records
  -> coordinate-merge all input streams
  -> form genomic bundles
  -> collapse equivalent transcript structures
  -> build strand-specific splice graphs
  -> map transcripts to graph paths
  -> combine compatible paths and remove redundant predictions
  -> filter and write merged GTF
```

### 1. Input preparation

- `TInputFiles::start()` in `tmerge.cpp` accepts individual input files or a
  file containing their names, prepares all inputs, and opens their readers.
- `TInputFiles::convert2BAM()` converts each GTF/GFF transcript into a temporary
  BAM alignment. This lets merge mode reuse the alignment-processing pipeline.
- `TInputFiles::next()` performs the coordinate-ordered multiway merge across
  the input readers.

### 2. Genomic bundles

- `main()` in `stringtie.cpp` scans the merged stream and divides transcripts
  into genomic bundles.
- `BundleData::evalReadAln()` extracts exon, junction, and guide-overlap data.
- `processRead()` stores and consolidates each transcript-shaped alignment.
- `processBundle()` sends each completed bundle through inference and output.

### 3. Core merge algorithm

`infer_transcripts()` selects `build_merge()` instead of the normal
`build_graphs()` assembly path.

`build_merge()` in `rlink.cpp` is the merge-mode driver. It:

1. Adds reference guides as transcript records with `guide_to_read()`.
2. Collapses non-guide transcripts with matching splice junction chains and
   combines their abundance information.
3. Forms overlapping transcript groups with `add_read_to_group()` and
   `merge_read_to_group()`.
4. Assigns groups by strand and creates connected sub-bundles with
   `create_bundle()` and `add_group_to_bundle()`.
5. Stores unstranded, single-exon sub-bundles directly as predictions.
6. Processes stranded sub-bundles through splice graphs and path merging.

The main deeper workers for stranded transcripts are:

- `create_graph()`: Constructs splice-graph nodes and edges from exon boundaries
  and junctions.
- `construct_treepat()`: Builds the pattern lookup used when mapping transcripts
  to graph paths.
- `get_read_to_transfrag()`: Maps each input transcript to a transfrag, the
  internal graph-path representation.
- `get_read_pattern()` and `update_abundance()`: Encode path membership and
  accumulate coverage or abundance from contributing inputs.
- `merge_transfrags()`: Sorts guides first and other paths by abundance, joins
  compatible partial paths, reconciles fuzzy ends, handles containment and
  retained-intron relationships, and removes redundant predictions.
- `store_merge_prediction()`: Converts a selected graph path into a
  `CPrediction`, including its exon structure, coverage, guide association, and
  optional contributing transcript names.

`process_merge_transfrags()` exists in `rlink.cpp`, but its call in
`build_merge()` is commented out. It is not part of the active merge path.

### 4. Final filtering and output

`printMergeResults()` sorts the predictions, removes weaker overlapping
non-guide predictions, joins eligible nearby single-exon predictions, preserves
matching guide identifiers, and writes the merged transcript and exon records.

## Active call path

```text
main
  -> TInputFiles::start / convert2BAM / next
  -> BundleData::evalReadAln
  -> processRead
  -> processBundle
       -> infer_transcripts
            -> build_merge
                 -> guide_to_read
                 -> add_read_to_group / merge_read_to_group
                 -> create_bundle / add_group_to_bundle
                 -> create_graph
                 -> construct_treepat
                 -> get_read_to_transfrag
                      -> get_read_pattern
                      -> update_abundance
                 -> merge_transfrags
                      -> add_transfrag_to_path
                      -> store_merge_prediction
       -> printMergeResults
```
