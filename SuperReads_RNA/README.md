# StringTie2 Super-Reads
Assembles, aligns, and pre-processes super-reads for use with StringTie2

## Requirements

* python
* perl
* HISAT2
* GMAP
* samtools

## Installation

Run `./install.sh`

## Usage

Requires short-read FASTQ file(s) to assemble and HISAT2 and GMAP indexes of the same reference genome

Example:

`./create_sr_rna.py -1 pair1.fq -2 pair2.fq -H hisat2_ref -G gmap_ref -o out_dir/`

Final super-read + short read alignments can be found in `out_dir/sr_merge.bam`

Use "-U" instead of "-1/-2" if reads unpaired.

If super-reads already assembled by MaSuRCA, can specify location of the "work1" directory using "-w" option.

Full usage:
```
usage: create_rna_sr.py [-h] [-1 SHORT_PAIR1] [-2 SHORT_PAIR2]
                        [-U SHORT_UNPAIRED] [-w WORK_DIR] -H HISAT_INDEX -G
                        GMAP_INDEX [-g GMAP_DIRECTORY] [-p NUM_THREADS]
                        [-o OUT_DIR] [--frag-len FRAG_LEN]
                        [--frag-std FRAG_STD] [--gmap-cmd GMAP_CMD]
                        [--hisat2-cmd HISAT2_CMD]

Quantifies Super-Reads

optional arguments:
  -h, --help            show this help message and exit
  -p NUM_THREADS, --num-threads NUM_THREADS
                        Number of threads each tool can use
  -1 SHORT_PAIR1, --short-pair1 SHORT_PAIR1
                        Paired short read FASTQ. Must also specify '--short-
                        pair2' and NOT '--short-unpaired'
  -2 SHORT_PAIR2, --short-pair2 SHORT_PAIR2
                        Paired short read FASTQ. Must also specify '--short-
                        pair1' and NOT '--short-unpaired'
  -U SHORT_UNPAIRED, --short-unpaired SHORT_UNPAIRED
                        Unpaired short read FASTQ. Must not include '--short-
                        pair1' or '--short-pair2'
  -w WORK_DIR, --work-dir WORK_DIR
                        MaSuRCA super-read assembly 'work1' directory. Will
                        perform super-read assembly if not included.
  -H HISAT_INDEX, --hisat-index HISAT_INDEX
                        HISAT2 index prefix for aligning short reads
  -G GMAP_INDEX, --gmap-index GMAP_INDEX
                        Gmap index for aligning super-reads
  -g GMAP_DIRECTORY, --gmap-directory GMAP_DIRECTORY
                        Gmap directory to find gmap index
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory
  --frag-len FRAG_LEN   Paired fragment length for MaSuRCA to use. Will only
                        be used if '--work-dir' not specify
  --frag-std FRAG_STD   Paired fragment standard deviation for MaSuRCA to use.
                        Will only be used if '--work-dir' not specify
  --gmap-cmd GMAP_CMD   Gmap executable
  --hisat2-cmd HISAT2_CMD
                        HISAT2 executable
```
