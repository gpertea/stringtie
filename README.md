![alt text](https://img.shields.io/badge/License-MIT-blue.svg "MIT License")

## StringTie: efficient transcript assembly and quantitation of RNA-Seq data

Stringtie employs efficient algorithms for transcript structure recovery and abundance estimation from bulk RNA-Seq reads aligned to a reference genome. 
It takes as input spliced alignments in coordinate-sorted SAM/BAM/CRAM format and produces a GTF output which consists of assembled 
transcript structures and their estimated expression levels (FPKM/TPM and base coverage values).

For additional StringTie documentation and the latest official source and binary packages please refer to the official website: <https://ccb.jhu.edu/software/stringtie>

## Obtaining and installing StringTie

Source and binary packages for this software can be directly downloaded from the [Releases](https://github.com/gpertea/stringtie/releases) page for this repository. 
StringTie is compatible with a wide range of Linux and Apple OS systems.
The main program (StringTie) does not have any other library dependencies (besides zlib) and in order to compile it from source it requires
a C++ compiler which supports the C++ 11 standard (GCC 4.8 or newer).

### Building the latest version from the repository 
In order to compile the StringTie source in this GitHub repository the following steps can be taken:
 
```
git clone https://github.com/gpertea/stringtie
cd stringtie
make -j4 release
```

To build with an alternate compiler set the `CC` and `CXX` environment variables, for example:
```
CC=clang CXX=clang++ make -j4 release
CC=icx CXX=icpx make -j4 release
```

During the first run of the above make command a few library dependencies will be downloaded and compiled, but any subsequent stringtie updates (using `git pull`)
should rebuild much faster.

To complete the installation, the resulting `stringtie` binary can then be copied to a programs directory of choice (preferably one that is in the current shell's  PATH).

Building and installing of StringTie this way should take less than a minute on a regular Linux or Apple MacOS 
desktop.

Note that simply running `make` would produce a less optimized executable which is suitable for debugging 
and runtime checking but that is significantly slower than the optimized version which 
is built by using the `make release` command as instructed above.

### Building and testing the offline source package

For HPC environments that do not have online access during the build, please download the latest `.offline.tar.gz` package from <a href="https://github.com/gpertea/stringtie/releases">Releases</a>. Unpack it and cd in the unpacked directory, then:

```
make -j6 release
make test
```
This should build the `stringtie` binary and run the included tests, without having to fetch dependencies and test data.  


### Using pre-compiled (binary) releases
Instead of compiling from source, some users may prefer to download an already compiled binary for Linux 
and Apple MacOS, ready to run. These binary package releases are compiled on older versions of these 
operating systems in order to provide compatibility with a wide range of OS versions not just the most recent distributions. 
These precompiled packages are made available on the <a href="https://github.com/gpertea/stringtie/releases">Releases</a> page for this repository.
Please note that these binary packages do not include the optional [super-reads module](#the-super-reads-module), 
which currently can only be built on Linux machines from the source made available in this repository.

## Running StringTie

The generic command line for the default usage has this format:
```
stringtie [-o <output.gtf>] [other_options] <read_alignments.bam> 
```
The main input of the program (_<read_alignments.bam>_) must be a SAM, BAM or CRAM file with RNA-Seq read 
alignments sorted by their genomic location (for example the `accepted_hits.bam` file produced
by TopHat, or HISAT2 output sorted with `samtools sort` etc.). 

The main output is a GTF file containing the structural definitions of the transcripts assembled by StringTie from the read alignment data. The name of the output file should be specified with the `-o` option. If this `-o` option is not used, the output GTF with the assembled transcripts will be printed to the standard 
output (and can be captured into a file using the `>` output redirect operator).

__Note__: if the `--mix` option is used, StringTie expects two alignment files to be given as positional parameters, in a specific order: the short read alignments must be the first file given while the long read alignments must be the second input file. Both alignment files must be sorted by genomic location.
```
stringtie [-o <output.gtf>] --mix [other_options] <short_read_alns.bam> <long_read_alns.bam> 
```

Note that the command line parser in StringTie allows arbitrary order and mixing of the positional parameters with the other options of the program, so the input alignment files can also precede or be given in between the other options -- the following command line is equivalent to the one above:

```
stringtie <short_read_alns.bam> <long_read_alns.bam> --mix [other_options] [-o <output.gtf>] 
```

### Nascent-aware assembly (new in StringTie3)

Many rRNA-depleted (“Total RNA”) libraries capture a mixture of mature and nascent (incomplete) RNA. StringTie 3 introduces a **nascent-aware** algorithm that accounts for this signal and separates mature from nascent RNA in the assembly process.

| Flag | Behaviour |
| :--- | :--- |
| `-N` | Enables nascent-aware assembly (recommended for any rRNA-depleted or Total RNA library). The primary output GTF contains **only mature** transcript models. |
| `--nasc` | Same algorithm as `-N`, but also retains nascent intermediates in the GTF output. Nascent records carry `nascentRNA` in the **source** column (field 2) and include a `nascent_parent "<mature_ID>"` attribute that links each nascent transcript to its mature parent. |

Full algorithmic details are described in our paper:  
Shinder I, Pertea G, Hu R, Rudnick Z, Pertea M. [*StringTie3 improves total RNA-seq assembly by resolving nascent and mature transcripts.*](https://doi.org/10.1038/s41592-026-03080-3) *Nature Methods* (2026). doi:10.1038/s41592-026-03080-3

### Running StringTie on the provided test data

When building from this source repository, after the program was compiled with `make release` as instructed above, the generated binary can be tested on a small data set with a command like this:
```
make test
```
This will run the included `run_tests.sh` script which downloads a small test data set 
and runs a few simple tests to ensure that the program works and generates the expected output.
It then runs `test_merge_policy.sh`, which checks the `--merge` behavior on small annotations
generated by the script itself (no download needed); `make test-merge` runs only this part.

If a pre-compiled package is used instead of compiling the program from source, the `run_tests.sh` script is included in the binary package as well and it can be run immediately after unpacking the binary package:

```
tar -xvzf stringtie-3.0.0.Linux_x86_64.tar.gz
cd stringtie-3.0.0.Linux_x86_64
./run_tests.sh
```

These small test/demo data set can be downloaded as <a href="https://github.com/gpertea/stringtie/raw/test_data/tests_3.tar.gz">tests_3.tar.gz</a> 
along with the source package and pre-compiled packages on the <a href="https://github.com/gpertea/stringtie/releases">Releases</a> 
page of this repository.

The tests can also be run manually as shown below (after changing to the _test_data_ directory, `cd test_data`):

#### Test 1: Input consists of only alignments of short reads
```
stringtie -o short_reads.out.gtf short_reads.bam
```

#### Test 2: Input consists of alignments of short reads and superreads
```
stringtie -o short_reads_and_superreads.out.gtf short_reads_and_superreads.bam
```
    
#### Test 3: Input consists of alignments of long reads
```
stringtie -L -o long_reads.out.gtf long_reads.bam
```
    
#### Test 4: Input consists of alignments of long reads and reference annotation (guides)
```
stringtie -L -G human-chr19_P.gff -o long_reads_guided.out.gtf long_reads.bam
```
#### Test 5: Input consists of alignments of short reads and alignments of long reads (using `--mix` option)
```
stringtie --mix -o mix_reads.out.gtf mix_short.bam mix_long.bam
```

#### Test 6: Input consists of alignments of short reads, alignments of long reads and a reference annotation (guides)
```
stringtie --mix -G mix_guides.gff -o mix_reads_guided.out.gtf mix_short.bam mix_long.bam
```

For version 3.0.0, three additional tests have been added, please see the `run_tests.sh` scripts for the details.

For large data sets one can expect up to one hour of processing time. A minimum of 8GB of RAM is recommended for running StringTie on regular size RNA-Seq samples, with 16 GB or more being strongly advised for larger data sets.

### StringTie options

The following optional parameters can be specified (use `-h` or `--help` to get the usage message):

```
 --version : print just the version at stdout and exit
 --conservative : conservative transcript assembly, same as -t -c 1.5 -f 0.05
 --mix : both short and long read data alignments are provided
        (long read alignments must be the 2nd BAM/CRAM input file)
 --rf : assume stranded library fr-firststrand
 --fr : assume stranded library fr-secondstrand
 -G reference annotation to use for guiding the assembly process (GTF/GFF)
 --ptf : load point-features from a given 4 column feature file <f_tab>:
         <ref_seq> <position> <strand> <type>, with <type> one of TSS, CPAS,
         JSTART, JEND; TSS and CPAS need a strand (+ or -)
 -o output path/file name for the assembled transcripts GTF (default: stdout)
 -l name prefix for output transcripts (default: STRG)
 -f minimum isoform fraction (default: 0.01)
 -L long reads processing; also enforces -s 1.5 -g 0 (default:false)
 -R if long reads are provided, just clean and collapse the reads but
    do not assemble
 -m minimum assembled transcript length (default: 200)
 -a minimum anchor length for junctions (default: 10)
 -j minimum junction coverage (default: 1)
 -t disable trimming of predicted transcripts based on coverage
    (default: coverage trimming is enabled)
 -c minimum reads per bp coverage to consider for multi-exon transcript
    (default: 1)
 -s minimum reads per bp coverage to consider for single-exon transcript
    (default: 4.75)
 -v verbose (log bundle processing details)
 -g maximum gap allowed between read mappings (default: 50)
 -M fraction of bundle allowed to be covered by multi-hit reads (default:1)
 -p number of threads (CPUs) to use (default: 1)
 -A gene abundance estimation output file
 -N nascent aware assembly: use this option with rRNA-depleted RNAseq
     libraries (e.g. Total RNA, Ribo-Zero)
 --nasc : enables -N but also outputs the assembled nascent transcripts
 -E define window around possibly erroneous splice sites from long reads to
    look out for correct splice sites (default: 25)
 -B enable output of Ballgown table files which will be created in the
    same directory as the output GTF (requires -G, -o recommended)
 -b enable output of Ballgown table files but these files will be 
    created under the directory path given as <dir_path>
 -e only estimate the abundance of given reference transcripts (requires -G)
 --viral only relevant for long reads from viral data where splice sites
    do not follow consensus (default:false)
 -x do not assemble any transcripts on the given reference sequence(s)
 -u no multi-mapping correction (default: correction enabled)
 -h print this usage message and exit
 --ref/--cram-ref reference genome FASTA file for CRAM input
```

## Transcript merge mode

`--merge` reads assembled transcript annotations, not read alignments, and
produces a unified structurally non-redundant GTF. Inputs can be listed directly
or supplied as one text file containing one GTF/GFF path per line.

```sh
stringtie --merge [options] -o merged.gtf samples.list
stringtie --merge [options] -o merged.gtf sample1.gtf sample2.gtf
```

The default merge is structure-driven: `-F`, `-T`, and `-f` are all zero, and
TPM/FPKM values are not pooled across samples as path weights. Each input file
contributes at most one vote to an intron chain, however many of its
transcripts carry that chain.

### Changes in merge mode since version 3.0.3

Version 3.0.4 changes what `--merge` writes with default settings. No set of
options reproduces the 3.0.3 output exactly; the table lists what changed and
the closest way back.

```text
Change                                      Closest to the 3.0.3 behavior
------------------------------------------  ------------------------------------
Input abundance filters are off             -F 1 -T 1 -f 0.01
An intron chain needs 2 input files when    --min-support 1
4 or more files are merged
Every observed intron chain is kept unless  none; 3.0.3 kept only the paths its
a better supported transcript contains it   splice graph assembled
Paths are joined only through a shared      --min-introns 0
intron, not through single-exon overlap
A transcript end needs 2 files; a longer    --min-support-5p 1 --min-support-3p 1
single-file end is trimmed
Reference transcripts given with -G are     none; 3.0.3 could move their ends
written with their exact structure
Fusions of reference genes and bridges      --allow-fusions --allow-cds-bridges
between coding transcripts are rejected
Unstranded multi-exon inputs are passed     none; 3.0.3 wrote each of their exons
through unchanged                           as a separate single-exon transcript
Output is sorted by position; references    none
keep the order of the input files
Generated IDs are assigned in output order  none
and do not depend on -p
```

The merged GTF also carries new attributes (`merge_support` and related fields,
`class_code`, `reference_id`), described under "Merge report and output
evidence".

### General merge options

```text
Option            Default  Meaning
----------------  -------  -------------------------------------------------
-G FILE           unset    Include a GTF/GFF3 reference annotation.
-o FILE           stdout   Write the merged GTF to FILE.
-m N              50       Require input transcript length N or greater.
-c F              0        Require input cov F or greater when filtering.
-F F              0        Require input FPKM F or greater when filtering.
-T F              0        Require input TPM F or greater when filtering.
-f F              0        Require a within-file isoform fraction F.
-g N              250      Bundle transcripts separated by at most N bases.
-i                off      Retain qualifying retained-intron models.
-l LABEL          MSTRG    Prefix generated transcript and gene IDs.
```

An explicitly supplied `-F` or `-T` requires that attribute to be present.
`-f` is calculated within an input file and declared gene, or within the local
bundle if `gene_id` is missing. It uses TPM first, then FPKM, then cov. These
filters mark inputs as failing; a replicated structure may still be kept by
`--rescue-support`.

### How the merge decides

Every structure is identified the same way by every option. A multi-exon
transcript (MET) is its strand plus ordered intron chain; transcript ends are
not part of the identity. Single-exon transcripts (SET) support each other when
they overlap reciprocally by `--set-overlap` (pairwise, never transitively).
Support always counts distinct input files, not records or abundance.

The options act in layers. Each layer only removes or trims what the previous
layer passed:

```text
Layer          Question                        Options
-------------  ------------------------------  -----------------------------------
record         Is this input record usable?    -m -c -F -T -f
structure      Is the structure replicated?    --min-support --rescue-support
                                               --set-overlap
ends           Where does each end go?         --min-support-5p/3p --slop-bp-*
                                               --slop-frac-* --max-extension-*
guide          The guide, a variant, or new?   --guide-policy --guide-slop-bp-*
                                               --guide-replace-frac
joins          May observed chains be joined?  --path-mode --min-introns
                                               --min-piece-support --allow-*
redundancy     A fragment of a kept model?     --max-intron-overhang
                                               --guide-fragments -i
```

Every observed intron chain that passes the structure layer is a candidate
output. Paths that join observed chains are additions; rejecting a join never
removes the chains it was built from. Unstranded inputs are not placed in the
splice graph: they are emitted unchanged, deduplicated and counted by the same
identity rules.

### Structure and path evidence

```text
Option                    Default     Meaning
------------------------  ----------  --------------------------------------
--min-support N           1 or 2      Files required for an intron chain or SET.
                                      Default 2 when 4 or more files are
                                      merged, otherwise 1.
--rescue-support N        2           Files required to rescue a filtered tx.
--min-introns N           1           Exact shared introns required to join.
--path-mode MODE          shared      shared, replicated, or direct evidence.
--min-piece-support N     2           Files per intron in replicated mode.
--set-overlap F           0.80        Reciprocal overlap for SET support and
                                      deduplication.
```

`shared` allows compatible MET paths to be joined when they share
`--min-introns` exact introns. The default therefore does not extend a MET
through terminal SET overlap alone. `--min-introns 0` permits that legacy
behavior, but does not bypass reference gene, CDS-bridge, hard extension, or
support barriers. `replicated` additionally requires every intron of a joined
path to reach `--min-piece-support`. `direct` emits observed intron chains only
and makes no joins.

### Biological endpoint controls

```text
Option                    Default     Meaning
------------------------  ----------  --------------------------------------
--min-support-5p N        2           Files required to support a 5-prime end.
--min-support-3p N        2           Files required to support a 3-prime end.
--slop-bp-5p N            100         Unsupported 5-prime overhang allowance.
--slop-bp-3p N            500         Unsupported 3-prime overhang allowance.
--slop-frac-5p F          0.50        5-prime slop as terminal exon fraction.
--slop-frac-3p F          0.50        3-prime slop as terminal exon fraction.
--max-extension-5p N      unlimited   Absolute 5-prime extension cap.
--max-extension-3p N      unlimited   Absolute 3-prime extension cap.
```

The 5-prime and 3-prime controls follow transcript orientation, including on
the minus strand. One model is emitted per intron chain. Each of its ends is
elected from the files that carry the chain: the end is the most extended
position reached by at least `--min-support-5p` (or `-3p`) files. A chain seen
in fewer files than that keeps its most extended end. A further end seen in
fewer files is adopted only when it lies within both the base-pair slop and the
terminal-exon-fraction slop of the elected end, and never beyond the hard cap.
`--slop-bp-*` accepts `unlimited`. `--min-support-5p 1 --min-support-3p 1`
restores the longest observed ends. Ends that are not adopted appear in the
report as `trim` rows.

### Reference and bridge policies

```text
Option                    Default     Meaning
------------------------  ----------  --------------------------------------
--guide-policy MODE       both        both, reference, or best.
--guide-replace-frac F    0.90        Core support fraction for best mode.
--guide-slop-bp-5p N      slop-bp-5p  Overhang past a guide 5-prime end that
                                      still counts as the guide.
--guide-slop-bp-3p N      slop-bp-3p  Same for the guide 3-prime end.
--guide-fragments MODE    absorb      absorb or dominance (see below).
--max-intron-overhang N   unlimited   Largest overhang into an intron of a
                                      containing transcript for a fragment
                                      to be absorbed; 0 keeps them separate.
--dedup-guides            off         Canonicalize identical guide structures.
--allow-fusions           off         Permit paths crossing guide genes.
--allow-cds-bridges       off         Permit bridges between coding guide tx.
```

Guide structures are immutable and are always emitted under their own ID. An
input related to a guide (the same intron chain, or a contiguous part of it) is
treated as follows:

- If it does not pass a guide end by more than `--guide-slop-bp-*`, or if fewer
  than `--min-support-5p/3p` files reach that far, it is the guide or a
  fragment of it. It is absorbed, its files are credited to the guide
  (`merge_support` for the guide's own chain, `merge_contained_support` for
  fragments), and fragments are listed in the report as `absorb` rows.
- Otherwise it is a guide variant, emitted under a generated ID with `cmp_ref`
  naming the guide and a `class_code` (see below).
  `both` keeps the guide and the variant; `reference` suppresses the variant;
  `best` replaces the guide with a same-chain variant only when the variant has
  at least two complete-path files and its complete-path fraction reaches
  `--guide-replace-frac`. `--max-extension-*` rejects a variant that passes
  the guide end by more than the cap.

With `--dedup-guides`, the lexicographically smallest ID represents identical
guide structures and the other IDs appear in `reference_aliases`.

A generated transcript whose intron chain is part of a longer retained
transcript is redundant and is absorbed, but only when the longer transcript's
weakest intron is supported by at least as many files as the shorter chain
(dominance). A replicated isoform is therefore never removed in favor of a
weaker longer model. Guides absorb their fragments regardless of support unless
`--guide-fragments dominance` is given. A fragment that retains an intron of
the longer transcript follows the same rule, and `-i` keeps all of them.

Same-strand transcripts whose exons overlap exons of distinct guide genes are
rejected (`cross_gene_bridge`) unless `--allow-fusions` is set. The barrier
follows the annotation: guide genes whose spans overlap on the same strand are
one meta-gene and are not a barrier to each other. That covers overlapping gene
families, a gene nested in another, a small RNA gene inside its host gene, and
an annotated readthrough gene together with both of its parents. An input with
exactly the intron chain of a guide always supports that guide, even when its
end runs into a neighboring gene. Within one guide gene, a path carrying
unique introns from two coding guide transcripts is rejected unless
`--allow-cds-bridges` is set. Opposite-strand guide overlap is reported but is
not a barrier.

### Merge report and output evidence

`--report FILE` writes a deterministically sorted TSV of filtering, rescue,
fusion, CDS bridge, containment, duplicate, support, `absorb` and `trim`
decisions, preceded by reason counts. An `absorb` row gives the number of input
records and files of the absorbed chain and the transcript that contains it.
The merged GTF records `merge_support` (files with evidence for the transcript),
`merge_path_support` (files carrying its complete intron chain),
`merge_piece_support` (files supporting its weakest intron),
`merge_5p_support` and `merge_3p_support` (files reaching the emitted end
within the slop), locus support, and optional `merge_contained_support` and
`reference_aliases`.

Relation to the reference annotation. A reference transcript carries
`reference_id` (its own ID) and `class_code "="`. A generated transcript that
overlaps the locus of any guide carries `cmp_ref` (the best matching guide) and
the gffcompare `class_code` of that relation, computed with the same code and
priorities as gffcompare (`~` same intron chain with other ends, `c`, `k`, `m`,
`n`, `j`, `e`, `o`, and for relations without a same-strand exon overlap `s`,
`x`, `i`, `y`). Among guides with equally ranked codes the one sharing more
splice sites, then more exon bases, is taken. There is one departure from the
gffcompare ranking: for a multi-exon transcript, `k` against a single-exon
guide (for example a miRNA gene inside one of its exons) gives way to a guide
that shares a splice site with it (`j`, `m`, `n`) and more exon bases. When the relation is a
same-strand exon overlap the transcript and its exons also get the guide's
`ref_gene_id` and `gene_name`. A transcript with no guide at its locus has none
of these attributes.

The report header also accounts for every input record. `# input_records` is
the number of transcripts read from the input files, and each `# fate` line
counts the records that ended one way; the fates sum to the total, and
`# fate unaccounted` is the remainder (0 when every record is explained).

```text
Fate                      Meaning
------------------------  -------------------------------------------------
emitted_chain             Multi-exon; its intron chain is an output
                          transcript (a reference transcript included).
emitted_set               Single-exon; an output single-exon transcript
                          represents it.
contained_fragment        Chain absorbed by a containing transcript (the
retained_intron_fragment  `records=` values of the `absorb` rows).
contained_set             Single-exon, inside an exon of an output
                          multi-exon transcript.
set_exon_overlap          Single-exon, overlapping such an exon partly.
minimum_support           Below `--min-support`.
abundance_filter          Failed `-c`/`-F`/`-T`/`-f` and was not rescued.
cross_gene_bridge         Rejected by the reference gene or coding
coding_transcript_bridge  transcript barrier (genes at separate loci).
direct_path_support       Rejected by `--path-mode`.
piece_support
guide_extension_cap       Variant of a reference transcript rejected by
guide_policy_reference    `--max-extension-*` or `--guide-policy`.
guide_policy_best
min_length                Shorter than `-m`.
excluded_reference        On a reference excluded with `-x`.
```

Multi-exon records are followed through the decision taken for their intron
chain. Stranded single-exon records are classified by the output transcript
that covers them, because the splice graph joins them before any decision.

Output order and IDs. Transcripts are written sorted by reference, start and
end, each followed by its exons. References appear in the order in which the
input files list them (first file first), followed by references that only the
`-G` annotation has. A reference transcript without input support is written
in place. Generated IDs (`<label>.<gene>.<n>`) are numbered in output order,
so the output is identical for any `-p`. A gene number used by a reference ID
of the form `<label>.<N>` or `<label>.<N>.<M>` is skipped.

All long options in this section are valid only with `--merge`. See
[MERGE_MODE.md](MERGE_MODE.md) for the implementation workflow, precise output
fields, and active call path. More details about other StringTie options are in
the [online manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).

## Input files

StringTie takes as input a SAM, BAM or CRAM file sorted by coordinate (genomic location). 
This file should contain spliced RNA-seq read alignments such as the ones produced by TopHat or HISAT2.
TopHat output is already sorted. Unsorted SAM or BAM files generated by other aligners should be sorted using the `samtools` program:
```
samtools sort -o alns.sorted.bam alns.sam
```
The file resulted from the above command (alns.sorted.bam) can be used as input to StringTie. 

Any SAM record with a spliced alignment (i.e. having a read alignment across at least one junction)
should have the `XS` tag to indicate the transcription strand, i.e. the genomic strand from which the RNA that produced
this read originated. TopHat and HISAT2 alignments already include this tag, but if you use
a different read mapper you should check that this tag is also included for spliced alignment
records. STAR aligner should be run with the option `--outSAMstrandField intronMotif` in order to generate this tag.

The `XS` tags are not necessary in the case of long RNA-seq reads aligned with `minimap2` using the `-ax splice` option. `minimap2` adds the `ts` tags to splice alignments to indicate the transcription strand (albeit in a different manner than the `XS` tag), and StringTie can recognize the `ts` tag as well, if the `XS` tag is missing. 
Thus the long read spliced alignments produced by `minimap2` can be also assembled by StringTie (with the option `-L` or 
as the 2nd input file for the `--mix` option).

As explained above, the alignments must be sorted by coordinate before they can be used as input for StringTie.

When CRAM files are used as input, the original reference genomic sequence can be provided with the `--ref` (`--cram-ref`) option 
as a multi-FASTA file with the same chromosome sequences that were used when aligning the reads. This is optional but recommended because StringTie
can better estimate the quality of some spliced alignments (e.g. noticing mismatches around junctions) and that data can be retrieved 
in the case of some CRAM files only when the reference genome sequence is also provided. 

### Reference transcripts (guides)

A reference annotation file in GTF or GFF3 format can be provided to StringTie 
using the `-G` option which can be used as 'guides' for the assembly process. 

When the `-e` option is used (i.e. expression estimation only), this option is required, 
and in that case StringTie will not attempt to assemble the read alignments but instead it will 
only estimate the expression levels of all the transcripts provided in this file

Note that when a reference transcript is fully covered by reads, the original transcript ID from the reference annotation file will be 
shown in StringTie's output record in the _`reference_id`_ GTF attribute. Output transcripts that lack such `reference_id` attribute 
can be considered "novel" transcript structures with respect to the given reference annotation.

### Point features

Known transcript ends and splice sites can be given with the `--ptf` option, in a text file
with four columns separated by tabs or spaces, one feature per line:

```
<ref_seq> <position> <strand> <type>
chr19   571311  +  TSS
chr19   583492  +  CPAS
chr19   571735  -  TSS
chr19   567222  -  CPAS
```

The position is 1-based, the strand is `+`, `-` or `.` and the type is one of:

 * `TSS`: transcription start site, the first base of a transcript
 * `CPAS`: cleavage and polyadenylation site, the last base of a transcript
 * `JSTART`: the last base of the exon before an intron
 * `JEND`: the first base of the exon after an intron

`TSS` and `CPAS` are the biological ends of a transcript, so they need a strand, and on the `-` strand
the `TSS` is the end with the higher coordinate, as in the example above.

How these features are used:

 * A `TSS` or `CPAS` adds a possible transcript start or end at that position, for transcripts on
   the same strand. It has to fall inside a region covered by reads: it can shorten a transcript, or
   start or end one within an exon, but it never extends a transcript beyond the read coverage.
 * With long reads (`-L`) transcripts starting or ending near such a feature are adjusted to it.
   With short reads, and with `--mix`, the feature is used if the read coverage supports
   a transcript start or end there.
 * `JSTART` and `JEND` are only used when assembling short reads alone (no `-L` or `--mix`):
   in a locus that has such features, only the junctions with both ends given this way
   (or matching a reference transcript given with `-G`) are kept. Their strand is not used.

## The super-reads module

This optional module can be used to de-novo assemble, align and pre-process
RNA-Seq reads, preparing them to be used as "super-reads" by Stringtie.

More usage information is provided in <a href="https://github.com/gpertea/stringtie/blob/master/SuperReads_RNA/README.md">SuperReads_RNA/README.md</a>.
Quick installation instructions for this module from the source available on this repository (assuming main Stringtie installation was already completed as described above):

```
 cd SuperReads_RNA
 ./install.sh
```

### Using super-reads with Stringtie

After running the super-reads module (see the <a href="https://github.com/gpertea/stringtie/blob/master/SuperReads_RNA/README.md">SuperReads_RNA</a> module documentation for usage details), there 
is a BAM file created which contains sorted alignment for both short reads and super-reads, called *`sr_merge.bam`*, 
created in the selected output directory. This file can be directly given as the main input file
to StringTie as described in the [Running StringTie](#running-stringtie) section above.


## License
StringTie is free, open source software released under an <a href="https://opensource.org/licenses/MIT">MIT License</a>.

## Publications
Shinder I, Pertea G, Hu R, Rudnick Z, Pertea M [**StringTie3 improves total RNA-seq assembly by resolving nascent and mature transcripts**](https://www.nature.com/articles/s41592-026-03080-3), _Nature Methods_ (2026), doi:10.1038/s41592-026-03080-3

Shumate A, Wong B, Pertea G, Pertea M [**Improved transcriptome assembly using a hybrid of long and short reads with StringTie**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009730), _PLOS Computational Biology_ 18, 6 (2022), doi.org/10.1371/journal.pcbi.1009730 

Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M  [**Transcriptome assembly from long-read RNA-seq alignments with StringTie2**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1), _Genome Biology_ 20, 278 (2019),  doi:10.1186/s13059-019-1910-1

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL [**Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown**](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html), _Nature Protocols_ 11, 1650-1667 (2016), doi:10.1038/nprot.2016.095

Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT  & Salzberg SL [**StringTie enables improved reconstruction of a transcriptome from RNA-seq reads**](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html), _Nature Biotechnology_ (2015), doi:10.1038/nbt.3122
