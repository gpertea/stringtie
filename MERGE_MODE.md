# StringTie merge mode

`stringtie --merge` combines transcript annotations from multiple GTF/GFF
inputs into a structurally non-redundant isoform set. It operates on assembled
transcripts, not raw RNA-seq reads.

## Workflow

```text
Input GTF/GFF files
  -> preserve tx ID, input file, gene ID, CDS, and available abundance fields
  -> convert transcripts to temporary alignment records
  -> coordinate-merge streams and form genomic bundles
  -> reject input bridges across distinct -G genes
  -> count intron chains (and overlapping SETs) by distinct input file
  -> apply optional abundance filters, with replicated rescue
  -> pass unstranded inputs through unchanged
  -> record every observed stranded intron chain with its per-file ends
  -> build strand-specific splice graphs and transfrags
  -> combine intron-compatible paths under the selected evidence policy
  -> make one candidate per observed chain and elect its ends
  -> keep graph paths only as joins of observed chains
  -> enforce guide, CDS-bridge, endpoint, and path-support barriers
  -> remove redundant fragments under the dominance rule, with report rows
  -> deduplicate generated structures and write merged GTF/report
```

### 1. Input preparation and filtering

- `TInputFiles::start()` accepts individual files or a file containing their
  names. `TInputFiles::next()` performs the coordinate-ordered multiway merge.
- `TInputFiles::convert2BAM()` preserves the input exon structure and carries
  transcript ID, input-file index, gene ID, CDS bounds, cov, FPKM, and TPM in
  the temporary record.
- `processRead()` records whether an input passes explicit `-c`, `-F`, or `-T`
  filters. It defers removal until structural support is known.
- `filter_merge_inputs()` counts each structure once per input file. A MET is
  identified by strand and intron chain (`merge_read_key()`), not by its ends;
  SETs support each other by pairwise reciprocal overlap of `--set-overlap`.
  `--min-support` removes insufficiently replicated structures and
  `--rescue-support` can restore a structure rejected only by abundance filters.
- Unstranded inputs are removed from graph construction in `build_merge()` and
  emitted as their own predictions, one per intron chain or per SET.
- `build_merge()` records each observed stranded intron chain in
  `BundleData::mergechains` (`CMergeChain`: introns, files, and the most
  extended start and end seen in each file) before reads are collapsed.
- An explicit `-F` or `-T` requires that field to be present and pass. An
  explicit `-f` is evaluated within input file and input gene, or within the
  current locus when `gene_id` is absent. It prefers TPM, then FPKM, then cov.

The merge defaults are `-F 0`, `-T 0`, and `-f 0`. Cross-sample TPM/FPKM values
are not used as splice-graph weights; each input file contributes one structural
vote.

### 2. Genomic bundles and core merge

`main()` forms bundles and `processBundle()` calls `infer_transcripts()`, which
selects `build_merge()` in merge mode.

`build_merge()` is the main driver. It adds guides with `guide_to_read()`,
collapses exact duplicate input structures, groups overlapping transcripts,
creates strand-specific bundles, and sends graph-backed bundles through these
deeper workers:

- `create_graph()`: builds splice-graph nodes and edges.
- `construct_treepat()`: builds the graph-pattern lookup.
- `get_read_to_transfrag()`: maps each input tx to a graph path.
- `get_read_pattern()` and `update_abundance()`: encode paths and structural
  votes.
- `merge_transfrags()`: joins compatible paths subject to intron and endpoint
  evidence.
- `add_transfrag_to_path()`: extends a candidate graph path.
- `store_merge_prediction()`: creates a `CPrediction`.
- `merge_set_observed_ends()`: restores input terminal coordinates when graph
  nodes were widened by an overlapping SET.
- `set_merge_prediction_support()`: records distinct-file, whole-path,
  intron-piece, and biological 5-prime/3-prime support.

`process_merge_transfrags()` is present but not called by the active path.

### 3. Structural policies

The default shared-path policy requires at least one exact shared intron before
two MET paths can be combined. All overlapping introns must remain compatible.
SETs therefore do not extend METs merely through terminal-exon overlap. SETs
are deduplicated by reciprocal overlap without constructing a transitive union;
a non-guide SET fully contained in a MET exon is removed.

The intron and endpoint checks apply both while `merge_transfrags()` builds a
path and while its legacy local-prediction cleanup considers absorbing one model
into another. This closes a former path around the policy. Prediction ends are
then clamped to observed `CTransfrag::longstart` and `longend` coordinates so a
graph node widened by a neighboring SET cannot silently widen an unrelated MET.

`--min-introns 0` restores terminal-overlap path joining, but it never bypasses
guide-gene, CDS-bridge, hard endpoint, or support barriers. Guide-gene barriers
have precedence over shared-intron evidence.

Biological 5-prime and 3-prime ends are handled independently. A new extension
must meet its end-specific support threshold unless both its base-pair and
terminal-exon-fraction overhang are within the configured slop. Hard extension
caps apply regardless of support.

The path modes are:

```text
Mode        Requirement
----------  --------------------------------------------------------------
shared      Compatible paths share at least --min-introns exact introns
replicated  Every output intron has --min-piece-support input files
direct      The complete ordered intron chain has --min-support input files
```

### 4. Guide and bridge policies

`-G` guide exon structures are loaded without close-exon normalization and are
never rewritten. A changed model never reuses a guide transcript ID.

Barrier order is:

1. Reject a same-strand input or path whose exons overlap exons of guides from
   more than one meta-gene, unless `--allow-fusions` is set.
2. Reject a same-gene path combining unique introns from distinct coding guides
   unless `--allow-cds-bridges` is set.
3. Report opposite-strand guide overlaps without treating them as a barrier.
4. Apply biological 5-prime/3-prime hard extension caps.
5. Apply shared-intron, path, piece, and endpoint support policies.

Meta-genes (`merge_guide_metagenes()`): the guides of a bundle are grouped by
union-find, joining guides whose spans overlap on the same strand and guides
with the same `gene_id`. Overlapping gene families, nested genes, small RNA
genes inside a host gene, and a readthrough gene with its parents therefore
form one meta-gene, and crossing between them is not a fusion. Genes at
non-overlapping loci remain separate barriers. `filter_merge_inputs()` applies
the rule per input record and never rejects a record whose intron chain is a
guide's chain; `merge_generated_barriers()` applies it to synthesized paths.

Guide output policies are:

```text
Policy     Result for a same-chain endpoint variant
---------  ---------------------------------------------------------------
both       Keep the unchanged guide and emit the variant under a new ID
reference  Keep the unchanged guide and suppress the endpoint variant
best       Replace the guide only with replicated, high-core-support evidence
```

Exact guide matches carry `reference_id` and `class_code "="`. Every other
output transcript is classified at print time by `merge_class_code()`, which
runs gclib `getOvlData()` (strict matching) against the guides of the bundle
that its span overlaps and keeps the best `classcode_rank()`; opposite-strand
exon overlaps become `s` (intron match) or `x`. Ties in rank go to the guide
with more matching splice sites, then more overlapping exon bases; and for a
multi-exon transcript a `k` against a single-exon guide yields to a
splice-sharing `j`/`m`/`n` relation with more exon overlap. It writes `cmp_ref` and
`class_code`, and `ref_gene_id`/`gene_name` for codes ranked below
`CLASSCODE_OVL_RANK`. The internal `merge_class`/`merge_ref` fields still
drive `--guide-policy` and are no longer printed.
Identical guide structures keep every reference ID by default. With
`--dedup-guides`, the lexicographically smallest ID is kept and the others are
reported in `reference_aliases`. Generated structures are always deduplicated.

### 5. Output and report

`printMergeResults()` applies final barriers, structural deduplication, guide
identity rules, and SET containment before writing GTF. Merge transcript
records include:

```text
merge_support          distinct files contributing to the prediction
merge_path_support     files containing the complete intron chain
merge_piece_support    minimum distinct-file support over output introns
merge_5p_support       files supporting the biological 5-prime endpoint
merge_3p_support       files supporting the biological 3-prime endpoint
merge_locus_support    distinct files represented in the locus
merge_5p_fraction      5-prime support / merge_support
merge_3p_fraction      3-prime support / merge_support
reference_id           the guide itself (class_code "=")
cmp_ref                best matching guide of a generated transcript
class_code             gffcompare class code of that relation
ref_gene_id, gene_name gene of the guide, for same-strand exon overlaps
reference_aliases      suppressed identical guide IDs, when requested
```

`--report FILE` writes a deterministic TSV after a comment summary. It records
filtering, replicated rescue, fusion and CDS-bridge decisions, opposite-strand
overlaps, containment, duplicate suppression, support, guide/gene IDs, input
transcript ID, and input file number.

## Merge-only controls

```text
Option                    Default
------------------------  ----------
--min-support             2 for 4 or more input files, else 1
--rescue-support          2
--min-introns             1
--path-mode               shared
--min-piece-support       2
--set-overlap             0.80
--min-support-5p          2
--min-support-3p          2
--slop-bp-5p              100
--slop-bp-3p              500
--slop-frac-5p            0.50
--slop-frac-3p            0.50
--max-extension-5p        unlimited
--max-extension-3p        unlimited
--guide-policy            both
--guide-replace-frac      0.90
--guide-slop-bp-5p        same as --slop-bp-5p
--guide-slop-bp-3p        same as --slop-bp-3p
--guide-fragments         absorb
--max-intron-overhang     unlimited
--dedup-guides            off
--allow-fusions           off
--allow-cds-bridges       off
--report                  unset
```

All controls above are rejected outside `--merge` mode.

### Support and assembly options

- `--min-support N` requires an intron chain (or a SET by reciprocal overlap)
  to occur in at least `N` distinct input files. Transcript ends do not have to
  agree. When the option is omitted, `N` is 2 if 4 or more files are merged
  and 1 otherwise; the value used is written to the `# min_support` line of
  the report.
- `--rescue-support N` restores a structure when every observed copy
  failed an explicitly requested abundance filter but the structure occurs in
  at least `N` files. It does not rescue a structural or guide barrier.
- `--min-introns N` requires `N` exact shared introns before compatible MET
  paths can be joined. The default `1` blocks SET-only terminal overlap;
  `0` deliberately restores terminal-overlap joining.
- `--path-mode shared` uses the shared-intron gate. `replicated` also requires
  every output intron to be replicated. `direct` requires replication of the
  complete ordered intron chain.
- `--min-piece-support N` is the per-intron distinct-file threshold used only
  by `--path-mode replicated`.
- `--set-overlap F` sets the reciprocal overlap fraction for choosing one of
  two generated SETs. It does not create a union of their ranges.

### Endpoint options

- `--min-support-5p N` and `--min-support-3p N` set independent distinct-file
  thresholds for an endpoint introduced by path joining.
- `--slop-bp-5p N` and `--slop-bp-3p N` set the maximum unsupported overhang in
  bases. An under-supported join must also pass the corresponding fractional
  slop.
- `--slop-frac-5p F` and `--slop-frac-3p F` set the maximum unsupported
  overhang as a fraction of the candidate terminal exon. Both the base-pair and
  fraction conditions must pass.
- `--max-extension-5p N` and `--max-extension-3p N` are absolute hard caps.
  They apply regardless of support; omission means unlimited.

Endpoint names are biological, not genomic: genomic right is 3-prime on `+`
and 5-prime on `-`. Support and slop govern synthesized joins. A directly
observed same-chain guide variant is controlled by the hard caps and guide
policy because its full endpoint is itself observed.

### Guide and barrier options

- `--guide-policy both` preserves an unchanged guide and emits an accepted
  same-chain endpoint variant under a generated ID.
- `--guide-policy reference` preserves the guide and suppresses its same-chain
  endpoint variant.
- `--guide-policy best` replaces the guide only when at least two files contain
  the variant's complete chain and complete-path support divided by total
  variant support reaches `--guide-replace-frac`.
- `--guide-replace-frac F` sets that `best` core-support fraction.
- `--dedup-guides` emits one canonical record for structurally identical guides.
  The lexicographically smallest transcript ID wins and aliases are retained in
  `reference_aliases`.
- `--allow-fusions` changes a distinct-guide-gene crossing from rejection to an
  accepted, reported fusion. This barrier otherwise precedes intron evidence.
- `--allow-cds-bridges` permits an intra-gene path that combines guide-specific
  introns from two coding guide transcripts.
- `--report FILE` writes the sorted decision events and per-reason counts to a
  TSV file after all worker threads have completed.

## Active call path

```text
main
  -> TInputFiles::start / convert2BAM / next
  -> BundleData::evalReadAln
  -> processRead
  -> processBundle
       -> infer_transcripts
            -> build_merge
                 -> filter_merge_inputs / guide_to_read
                 -> add_read_to_group / merge_read_to_group
                 -> create_bundle / add_group_to_bundle
                 -> create_graph / construct_treepat
                 -> get_read_to_transfrag
                      -> get_read_pattern / update_abundance
                 -> merge_transfrags
                      -> add_transfrag_to_path
                      -> store_merge_prediction
                           -> merge_set_observed_ends
                           -> set_merge_prediction_support
       -> printMergeResults
```

## Regression tests

```sh
make test-merge
VALGRIND=1 ./test_merge_policy.sh
```

`test_merge_policy.sh` builds temporary locus-level GTFs and checks support
rescue, path modes, SET/MET terminal overlap, independent biological endpoints,
flexible guide caps, guide policies, gene and CDS bridge barriers, guide dedup,
and transcript-ID uniqueness. `compare_merge_outputs.pl` compares two merge
GTFs at exact-structure, intron-chain, guide-identity, and source-versus-
synthesized levels.

## Observed chains first

`printMergeResults()` applies the merge policy in this order. The splice-graph
code (`merge_transfrags()`) is unchanged; its output is reinterpreted.

1. Candidates. Each `CMergeChain` becomes a prediction (`merge_origin` `o`).
   `merge_elect_prediction_ends()` sets each end to the k-th most extended end
   over the files that carry the chain (k is `--min-support-5p` or `-3p`, or
   the file count when smaller), extends it to a further observed end only
   within `--slop-bp-*`, `--slop-frac-*` and `--max-extension-*`, and writes a
   `trim` report row for each file end that was not adopted.
2. Graph paths. A graph path whose intron chain was observed is superseded by
   its candidate. A path with an unobserved chain is a join (`merge_origin`
   `j`); its ends are elected from the observed chains that share its terminal
   intron and are part of it. A guide prediction whose ends the graph moved is
   restored to the guide structure, and guide support is recomputed from the
   observed chains.
3. Per-prediction rejections: minimum support, path mode, gene and CDS
   barriers, guide extension caps and `--guide-policy`.
   `merge_related_guide()` relates a prediction to a guide when it has the
   guide's chain, or a contiguous part of it, and extends a guide end beyond
   `--guide-slop-bp-*` with `--min-support-5p/3p` files
   (`merge_supported_extension()`).
4. Redundancy. `merge_remove_redundant()` runs once, after step 3, from the
   longest chains to the shortest. A generated MET is absorbed by a retained
   container when its chain is part of the container's chain, its ends fit
   (`merge_fragment_fits()`), and the container's `merge_piece_support` is at
   least the fragment's direct support. Guides are exempt from that last
   condition unless `--guide-fragments dominance` is set. Absorbed files are
   credited and an `absorb` row is written with record and file counts.
5. `merge_assign_genes()` groups the retained predictions into genes by
   same-strand exon overlap.

Because step 4 runs after step 3 and only considers retained containers, a
rejected join cannot remove the observed chains it was built from.

## Input record reconciliation

With `--report`, every input record gets one final fate and the report header
lists the counts (`# input_records`, `# fate <name> <count>`,
`# fate unaccounted <n>`). `mergeFateAdd()` holds the tally; `unaccounted` is
the record count minus the tallied fates, so a record that leaves the merge
without a decision shows up there. The regression suite requires 0 for every
fixture.

```text
Stage                       Where                        Fates
--------------------------  ---------------------------  ----------------------
input loop, processRead()   before the bundle            excluded_reference,
                                                         min_length
filter_merge_inputs()       per record                   cross_gene_bridge,
                                                         minimum_support,
                                                         abundance_filter
printMergeResults()         per observed chain or        emitted_chain,
  -> merge_tally_fates()    pass-through prediction      contained_fragment,
                                                         retained_intron_fragment,
                                                         rejection reasons
                            per stranded single-exon     emitted_set,
                            record, by position          contained_set,
                                                         set_exon_overlap
```

`printMergeResults()` keeps a `fate` label per prediction, set where the
prediction leaves the output, and `fateof` for a prediction dropped in favor of
an equivalent one (same guide chain, duplicate structure, overlapping
single-exon model), which takes that one's fate. A stranded multi-exon record
is counted through its `CMergeChain` (`nrec`), an unstranded record through its
pass-through prediction (`merge_nrec`). Stranded single-exon records enter the
splice graph and are not traced individually: each is assigned, in this order,
to an emitted single-exon transcript that overlaps it, an emitted multi-exon
transcript whose exon contains it, one whose exon overlaps it, or the fate of
a rejected single-exon prediction that overlaps it.

## Output order and IDs

`printMergeResults()` and `writeUnbundledGuides()` do not write to the output
file. Each transcript record (transcript line plus exon lines) goes to
`<tmpdir>/merge_out.gtf` and `mergeOutputAdd()` keeps an index entry for it:
reference rank, start, end, strand, guide ID, bundle start, gene and transcript
number within the bundle, file offset and length. In the record text the
generated gene ID is the single byte `MERGE_GENE_MARK`.

`writeMergeOutput()` runs after all bundles and unbundled guides. It sorts the
index, assigns gene numbers in sorted order (one per bundle-local gene,
skipping numbers reserved by `reserveMergeGeneNo()` for guide IDs that look
like generated IDs), and copies each record to the output with the mark
replaced by `<label>.<gene>`. The result does not depend on the order in which
worker threads finish bundles.

Reference rank is the order of first appearance: `TInputFiles::convert2BAM()`
registers the references of each input file in file order through
`mergeRefIndex()`; references seen only in the guides are ranked when their
unbundled guides are written.
