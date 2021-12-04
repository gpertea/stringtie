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
make release
```
During the first run of the above make command a few library dependencies will be downloaded and compiled, but any subsequent stringtie updates (using `git pull`) 
should rebuild much faster.

To complete the installation, the resulting `stringtie` binary can then be copied to a programs directory of choice (preferably one that is in the current shell's  PATH).

Building and installing of StringTie this way should take less than a minute on a regular Linux or Apple MacOS 
desktop.

Note that simply running `make` would produce a less optimized executable which is suitable for debugging 
and runtime checking but that is significantly slower than the optimized version which 
is built by using the `make release` command as instructed above.

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

### Running StringTie on the provided test/demo data

When building from this source repository, after the program was compiled with `make release` as instructed above, the generated binary can be tested on a small data set with a command like this:
```
make test
```
This will run the included `run_tests.sh` script which downloads a small test data set 
and runs a few simple tests to ensure that the program works and generates the expected output.

If a pre-compiled package is used instead of compiling the program from source, the `run_tests.sh` script is included in the binary package as well and it can be run immediately after unpacking the binary package:

```
tar -xvzf stringtie-2.2.0.Linux_x86_64.tar.gz
cd stringtie-2.2.0.Linux_x86_64
./run_tests.sh
```

These small test/demo data sets can also be downloaded separately as <a href="https://github.com/gpertea/stringtie/raw/test_data/tests.tar.gz">test_data.tar.gz</a> 
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

These tests should complete in several seconds.

For large data sets one can expect up to one hour of processing time. A minimum of 8GB of RAM is recommended for running StringTie on regular size RNA-Seq samples, with 16 GB or more being strongly advised for larger data sets.


### StringTie options

The following optional parameters can be specified (use `-h` or `--help` to get the usage message):

```
 --mix : both short and long read data alignments are provided
        (long read alignments must be the 2nd BAM/CRAM input file)
 --rf : assume stranded library fr-firststrand
 --fr : assume stranded library fr-secondstrand
 -G reference annotation to use for guiding the assembly process (GTF/GFF)
 --conservative : conservative transcript assembly, same as -t -c 1.5 -f 0.05
 --ptf : load point-features from a given 4 column feature file <f_tab>
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
 -E define window around possibly erroneous splice sites from long reads to
    look out for correct splice sites (default: 25)
 -B enable output of Ballgown table files which will be created in the
    same directory as the output GTF (requires -G, -o recommended)
 -b enable output of Ballgown table files but these files will be 
    created under the directory path given as <dir_path>
 -e only estimate the abundance of given reference transcripts (requires -G)
 --viral : only relevant for long reads from viral data where splice sites
    do not follow consensus (default:false)
 -x do not assemble any transcripts on the given reference sequence(s)
 -u no multi-mapping correction (default: correction enabled)
 --ref/--cram-ref reference genome FASTA file for CRAM input

Transcript merge usage mode: 

  stringtie --merge [Options] { gtf_list | strg1.gtf ...}
With this option StringTie will assemble transcripts from multiple
input files generating a unified non-redundant set of isoforms. In this mode
the following options are available:
  -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
  -o <out_gtf>     output file name for the merged transcripts GTF
                    (default: stdout)
  -m <min_len>     minimum input transcript length to include in the merge
                    (default: 50)
  -c <min_cov>     minimum input transcript coverage to include in the merge
                    (default: 0)
  -F <min_fpkm>    minimum input transcript FPKM to include in the merge
                    (default: 1.0)
  -T <min_tpm>     minimum input transcript TPM to include in the merge
                    (default: 1.0)
  -f <min_iso>     minimum isoform fraction (default: 0.01)
  -g <gap_len>     gap between transcripts to merge together (default: 250)
  -i               keep merged transcripts with retained introns; by default
                   these are not kept unless there is strong evidence for them
  -l <label>       name prefix for output transcripts (default: MSTRG)

```

More details about StringTie options can be found in the [online manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).

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
Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M  [**Transcriptome assembly from long-read RNA-seq alignments with StringTie2**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1), _Genome Biology_ 20, 278 (2019),  doi:10.1186/s13059-019-1910-1

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL [**Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown**](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html), _Nature Protocols_ 11, 1650-1667 (2016), doi:10.1038/nprot.2016.095

Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT  & Salzberg SL [**StringTie enables improved reconstruction of a transcriptome from RNA-seq reads**](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html), _Nature Biotechnology_ 2015, doi:10.1038/nbt.3122
