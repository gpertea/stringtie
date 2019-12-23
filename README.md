## Obtaining and installing StringTie

Source and binary packages for this software, along with a small test data set 
can be directly downloaded from the <a href="https://github.com/gpertea/stringtie/releases">Releases</a> page for this repository. 
StringTie is compatible with a wide range of Linux and Apple OS systems (going as far back as RedHat Enterprise Linux 5.0 and OS X 10.7). 
The main program (StringTie) does not have any other library dependencies and in order to compile it from source it requires only 
a C++ compiler which supports the C++ 0x standard (GCC 4.5 or newer).

### Building the latest version from the repository 
In order to compile the StringTie source in this GitHub repository the following steps can be taken:
 
```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release
```

If the compilation is successful, the resulting `stringtie` binary can then be copied to 
a programs directory of choice.

Installation of StringTie this way should take less than a minute on a regular Linux or Apple MacOS 
desktop.

Note that simply running `make` would produce an executable which is more suitable for debugging 
and runtime checking but which can be significantly slower than the optimized version which 
is obtained by using `make release` as instructed above.

### Using pre-compiled (binary) releases
Instead of compiling from source, some users may prefer to download an already compiled binary for Linux 
and Apple OS X, ready to run. These binary package releases are compiled on older versions of these 
operating systems (RedHat Enterprise Linux 5.0 and OS X 10.7) in order to provide compatibility with 
a wide range of (older) OS versions, not just the most recent versions. These precompiled packages are 
made available on the <a href="https://github.com/gpertea/stringtie/releases">Releases</a> page for this repository.
Please note that these binary packages do not include the optional [super-reads module](#the-super-reads-module), 
which currently can only be built on Linux machines, from the source made available in this repository.


## Running StringTie

Run stringtie from the command line like this:
```
stringtie [options] <aligned_reads.bam>
```
The main input of the program is a SAMTools BAM file with RNA-Seq mappings
sorted by genomic location (for example the accepted_hits.bam file produced
by TopHat).

The main output of the program is a GTF file containing the structural definitions of the transcripts assembled by StringTie from the read alignment data. The name of the output file should be specified by with the `-o` option.

### Running StringTie on the provided test/demo data
When building from this source repository, after the program was compiled with `make release` as instructed above, the generated binary can be tested on a small data set with a command like this:
```
make test
```
This will run the included `run_tests.sh` script which downloads a small test data set 
and runs a few simple tests to ensure that the program works and generates the expected output.

If a pre-compiled package is used instead of compiling the program from source, the `run_tests.sh` script is included in the binary package as well and it can be run immediately after unpacking the binary package:

```
tar -xvzf stringtie-2.0.Linux_x86_64.tar.gz
cd stringtie-2.0.Linux_x86_64
./run_tests.sh
```

These small test/demo data sets can also be downloaded separately as <a href="https://github.com/gpertea/stringtie/raw/test_data/test_data.tar.gz">test_data.tar.gz</a> 
along with the source package and pre-compiled packages on the <a href="https://github.com/gpertea/stringtie/releases">Releases</a> 
page of this repository.

The tests can also be run manually as shown below (after changing to the _test_data_ directory, `cd test_data`):

#### Run 1: Input consists of only alignments of short reads
```
stringtie -o short_reads.out.gtf short_reads.bam
```

#### Run 2: Input consists of alignments of short reads and superreads
```
stringtie -o short_reads_and_superreads.out.gtf short_reads_and_superreads.bam
```
    
#### Run 3: Input consists of alignments of long reads
```
stringtie -L -o long_reads.out.gtf long_reads.bam
```
    
#### Run 4: Input consists of alignments of long reads and reference annotation (guides)
```
stringtie -L -G human-chr19_P.gff -o long_reads_guided.out.gtf long_reads.bam
```

The above runs should take around one second each on a regular Linux or MacOS desktop. 
(see also <a href="https://github.com/gpertea/stringtie/blob/master/test_data/README.md">test_data/README.md</a>).

For very large data sets one can expect up to one hour of processing time. A minimum of 8GB of RAM is recommended for running StringTie on regular size RNA-Seq samples, with 16 GB or more being strongly advised for larger data sets.


### StringTie options

The following optional parameters can be specified (use -h/--help to get the
usage message):
```
 --version : print just the version at stdout and exit
 --conservative : conservative transcriptome assembly, same as -t -c 1.5 -f 0.05
 --rf assume stranded library fr-firststrand
 --fr assume stranded library fr-secondstrand
 -G reference annotation to use for guiding the assembly process (GTF/GFF3)
 -o output path/file name for the assembled transcripts GTF (default: stdout)
 -l name prefix for output transcripts (default: STRG)
 -f minimum isoform fraction (default: 0.01)
 -L use long reads settings (default:false)
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
 -B enable output of Ballgown table files which will be created in the
    same directory as the output GTF (requires -G, -o recommended)
 -b enable output of Ballgown table files but these files will be 
    created under the directory path given as <dir_path>
 -e only estimate the abundance of given reference transcripts (requires -G)
 -x do not assemble any transcripts on the given reference sequence(s)
 -u no multi-mapping correction (default: correction enabled)
 -h print this usage message and exit

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

## Input files

StringTie takes as input a binary SAM (BAM) file sorted by reference position. 
This file contains spliced read alignments such as the ones produced by TopHat or HISAT2.
A text file in SAM format should be converted to BAM and sorted using the 
samtools program:
```
samtools view -Su alns.sam | samtools sort - alns.sorted
```
The file resulted from the above command (alns.sorted.bam) can be used 
directly as input to StringTie. 

Any SAM spliced read alignment (a read alignment across at least one junction)
needs to contain the XS tag to indicate the strand from which the RNA that produced
this read originated. TopHat alignments already include this tag, but if you use
a different read mapper you should check that this tag is also included for spliced alignment
records. For example HISAT2 should be run with the `--dta` option in order to tag spliced 
alignments this way. As explained above, the alignments in SAM format should be sorted and
preferrably converted to BAM.

Optionally, a reference annotation file in GTF/GFF3 format can be provided to StringTie 
using the `-G` option. In this case, StringTie will check to see if the reference transcripts 
are expressed in the RNA-Seq data, and for the ones that are expressed it will compute coverage
and FPKM values.
Note that the reference transcripts should be fully covered by reads in order to be included
in StringTie's output with the original ID of the reference transcript shown in the 
_`reference_id`_ GTF attribute in the output file . Other transcripts assembled from 
the input alignment data by StringTie and not present in the reference file will be 
printed as well ("novel" transcripts).

## The super-reads module

This optional module can be used to de-novo assemble, align and pre-process
RNA-Seq reads, preparing them to be used as "super-reads" by Stringtie.

Mode detailed information is provided in the 
<a href="https://github.com/gpertea/stringtie/blob/master/SuperReads_RNA/README.md">SuperReads_RNA/README.md</a>.
Quick installation instructions for this module from the source available on this repository 
(assuming the above Stringtie installation was completed):

```
 cd SuperReads_RNA
 ./install.sh
```

### Using super-reads with Stringtie

After running the super-reads module (see the <a href="https://github.com/gpertea/stringtie/blob/master/SuperReads_RNA/README.md">SuperReads_RNA</a> module documentation for usage details), there 
is a BAM file which contains sorted alignment for both short reads and super-reads, called *`sr_merge.bam`*, 
created in the selected output directory. This file can be directly given as the main input file
to StringTie as described in the [Running StringTie](#running-stringtie) section above.


## License
StringTie is free, open source software released under an <a href="https://opensource.org/licenses/MIT">MIT License</a>.

## Publications
Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M  [**Transcriptome assembly from long-read RNA-seq alignments with StringTie2**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1), _Genome Biology_ 20, 278 (2019),  doi:10.1186/s13059-019-1910-1

Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL [**Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown**](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html), _Nature Protocols_ 11, 1650-1667 (2016), doi:10.1038/nprot.2016.095

Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT  & Salzberg SL [**StringTie enables improved reconstruction of a transcriptome from RNA-seq reads**](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html), _Nature Biotechnology_ 2015, doi:10.1038/nbt.3122
