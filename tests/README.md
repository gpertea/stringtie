## Test data

The test data can be automatically retrieved by the `run_tests.sh` script included 
with all source or binary distributions of StringTie, or downloaded separately from this url:
https://github.com/gpertea/stringtie/raw/test_data/tests.tar.gz

The `run_tests.sh` script will then run StringTie on these data sets and compare the output with the 
precomputed, expected output for each case. If the output of each test matches the 
expected output, the test is considered successful (and "OK." will be shown on the next line), 
otherwise an error will be reported.

## Running StringTie on the test data separately

The command lines shown below assume that the _stringtie_ program is installed somewhere 
in the shell's PATH. If that is not the case, the `stringtie` comand below should be prepended 
with the directory path to the _stringtie_ executable

### Test 1: Input consists of only alignments of short reads

```
stringtie -o short_reads.out.gtf short_reads.bam
```

### Test 2: Input consists of alignments of short reads and superreads

```
stringtie -o short_reads_and_superreads.out.gtf short_reads_and_superreads.bam
```
    
### Test 3: Input consists of alignments of long reads

```
stringtie -L -o long_reads.out.gtf long_reads.bam
```
    
### Test 4: Input consists of alignments of long reads and reference annotation (guides)

```
stringtie -L -G human-chr19_P.gff -o long_reads_guided.out.gtf long_reads.bam
```

### Test 5: Input consists of short read alignments and long read alignments:

```
stringtie --mix -o mix_reads.out.gtf mix_short.bam mix_long.bam
```

### Test 6: Input consists of short read alignments and long read alignments, with reference annotation (guides):

```
stringtie --mix -G mix_guides.gff -o mix_reads_guided.out.gtf mix_short.bam mix_long.bam
```

