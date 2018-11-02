## Obtaining and installing StringTie

The current version of StringTie can also be downloaded from
  http://ccb.jhu.edu/software/stringtie
  
In order to build StringTie from this GitHub repository
the following steps can be taken:
 
```
git clone https://github.com/gpertea/stringtie
cd stringtie
make release
```

Note that simply running `make` will produce an executable 
which is more suitable for debugging and runtime checking but which can be
significantly slower than the optimized version which is obtained by using 
`make release`.

## Running StringTie

Run stringtie from the command line like this:
'''
stringtie [options] <aligned_reads.bam>
'''
The main input of the program is a SAMTools BAM file with RNA-Seq mappings
sorted by genomic location (for example the accepted_hits.bam file produced
by TopHat).

A list of parameters and options for the program can be obtained by running
the program with only the -h or --help option.

More detailed usage instructions can be found on the StringTie website: 
http://ccb.jhu.edu/software/stringtie

