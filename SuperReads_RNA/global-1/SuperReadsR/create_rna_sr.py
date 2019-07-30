#!/usr/bin/env python

import argparse
import subprocess
import sys
import re
import os

#SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
BIN_DIR = #__#
WORKING_DIR = os.getcwd()

CONF_TEMPLATE = """DATA
PE = pe %d %d %s
END
PARAMETERS
STOP_AFTER_SUPERREADS=1
NUM_THREADS=%d
GRAPH_KMER_SIZE=auto
END
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantifies Super-Reads")
    parser.add_argument("-1", "--short-pair1", required=False, type=str, default=None, help="Paired short read FASTQ. Must also specify '--short-pair2' and NOT '--short-unpaired'")
    parser.add_argument("-2", "--short-pair2", required=False, type=str, default=None, help="Paired short read FASTQ. Must also specify '--short-pair1' and NOT '--short-unpaired'")
    parser.add_argument("-U", "--short-unpaired", required=False, type=str, default=None, help="Unpaired short read FASTQ. Must not include '--short-pair1' or '--short-pair2'")
    parser.add_argument("-w", "--work-dir", required=False, type=str, default=None, help="MaSuRCA super-read assembly 'work1' directory. Will perform super-read assembly if not included.")
    parser.add_argument("-H", "--hisat-index", required=True, type=str, help="HISAT2 index prefix for aligning short reads")
    parser.add_argument("-G", "--gmap-index", required=True, type=str, help="Gmap index for aligning super-reads")
    parser.add_argument("-g", "--gmap-directory", required=False, type=str, default=None, help="Gmap directory to find gmap index")
    parser.add_argument("-p", "--num-threads", required=False, type=int, default=1, help="Number of threads each tool can use")
    parser.add_argument("-o", "--out-dir", required=False, type=str, default="", help="Output directory")
    parser.add_argument("--frag-len", required=False, type=float, default=100, help="Paired fragment length for MaSuRCA to use. Will only be used if '--work-dir' not specify")
    parser.add_argument("--frag-std", required=False, type=float, default=20, help="Paired fragment standard deviation for MaSuRCA to use. Will only be used if '--work-dir' not specify")
    parser.add_argument("--gmap-cmd", required=False, type=str, default="gmap", help="Gmap executable")
    parser.add_argument("--hisat2-cmd", required=False, type=str, default="hisat2", help="HISAT2 executable")
    args = parser.parse_args()

    if not os.path.isdir(BIN_DIR):
        sys.stderr.write("Error: programs directory not found. Please run './install.sh''\n")
        sys.exit(1)
    
    if args.short_unpaired != None:
        short_reads = os.path.abspath(args.short_unpaired)
        paired = False
    elif args.short_pair1 != None and args.short_pair2 != None:
        short_reads = "%s %s" % (os.path.abspath(args.short_pair1), os.path.abspath(args.short_pair2))
        paired = True
    else:
        sys.stderr.write("Error: must specify short read input")
        sys.exit(1)


    if args.work_dir == None:
        if not os.path.exists(args.out_dir):
            try:
              os.makedirs(args.out_dir)
            except OSError, e:
              sys.stderr.write("\nError creating directory %s (%s)" % (args.out_dir, e));
              sys.exit(1)
        os.chdir(args.out_dir)
        conf_str = CONF_TEMPLATE % (
            args.frag_len, 
            args.frag_std,
            short_reads,
            args.num_threads)
        conf_file = open("super_reads.conf", "w")
        conf_file.write(conf_str)
        conf_file.close()
        subprocess.call([os.path.join(BIN_DIR, "createSuperReads_RNA"), "super_reads.conf"])
        subprocess.call("./assemble.sh")

        os.chdir(WORKING_DIR)
        work1_dir = os.path.join(args.out_dir, "work1")
    else:
        work1_dir = args.work_dir

    gmap_call = [args.gmap_cmd, 
                 "-f", "samse",
                 "-t", str(args.num_threads),
                 "-d", args.gmap_index]

    if args.gmap_directory != None:
        gmap_call += ["--dir", args.gmap_directory]
    gmap_call += ["%s/superReadSequences.fasta" % work1_dir]

    gmap_out_name = os.path.join(args.out_dir, "sr.sam")

    gmap_out = open(gmap_out_name, "w")
    subprocess.call(gmap_call, stdout=gmap_out)
    gmap_out.close()

    cmd = [os.path.join(BIN_DIR, "assign_reads"), 
                     "-p", str(args.num_threads),
                     "-w", work1_dir,
                     "-s", gmap_out_name,
                     "-o", os.path.join(args.out_dir, "sr_quant.sam")]
    subprocess.call(cmd)

    hisat2_call = [args.hisat2_cmd,
                  "-p", str(args.num_threads),
                  "-x", args.hisat_index]

    if paired:
        hisat2_call += ["-1", args.short_pair1, "-2", args.short_pair2]
    else:
        hisat2_call += ["-U", args.short_unpaired]

    hisat2_call += [">", os.path.join(args.out_dir, "short.sam")] 

    subprocess.call(hisat2_call)

    subprocess.call(["samtools", "merge", "-f",
                     "-@", str(args.num_threads),
                     os.path.join(args.out_dir, "sr_merge.bam"), 
                     os.path.join(args.out_dir, "sr_quant.sam"),
                     os.path.join(args.out_dir, "short.sam")])

    subprocess.call(["samtools", "sort", 
                     os.path.join(args.out_dir, "sr_merge.bam"),
                     "-o", os.path.join(args.out_dir, "sr_sort.bam")])
