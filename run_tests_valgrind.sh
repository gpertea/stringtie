#!/usr/bin/env bash

function unpack_test_data() {
  t=tests.tar.gz
  if [ ! -f $t ]; then
    echo "Error: file $t not found!"
    exit 1
  fi
  echo "..unpacking test data.."
  echo
  tar -xzf $t
  if [ ! -f tests/human-chr19_P.gff ]; then
     echo "Error: invalid test data archive?"
     exit 1
  fi
  
  #cp tests_exp_out/*.gtf tests/ 2>/dev/null
}

#if [ ! -f tests/human-chr19_P.gff ]; then
  if [ -d ./tests ]; then
    #extract the tarball and rename the directory
    echo "..Using existing ./tests"
    #unpack_test_data
  else
    echo "..Downloading test data.."
    #use curl to fetch the tarball from a specific github release or branch
    curl -ksLO https://github.com/gpertea/stringtie/raw/test_data/tests.tar.gz
    unpack_test_data
  fi
# fi
cd tests
# array element format:
# 
arrins=("short_reads" "short_reads_and_superreads" "long_reads" "long_reads" \
  "mix_short mix_long" "mix_short mix_long")
arrparms=("" "" "-L" "-L -G human-chr19_P.gff" "--mix" "--mix -G mix_guides.gff")
arrout=("short_reads" "short_reads_and_superreads" "long_reads" "long_reads_guided" \
  "mix_reads" "mix_reads_guided")
arrmsg=("Short reads"  "Short reads and super-reads" "Long reads" \
 "Long reads with annotation guides" "Mixed reads" "Mixed reads with annotation guides")
for i in ${!arrmsg[@]}; do
 fout="${arrout[$i]}.out.gtf"
 /bin/rm -f $fout
 fcmp="${arrout[$i]}.out_expected.gtf"
 if [ ! -f $fcmp ]; then
   echo "Error: file $fcmp does not exist! Re-download test data."
   exit 1
 fi
 n=$i
 ((n++))
 echo "Test ${n}: ${arrmsg[$i]}"
 fin=${arrins[$i]}.bam
 if [[ ${arrins[$i]} =~ ^mix ]]; then
   ins=( ${arrins[$i]} )
   fin="${ins[0]}.bam ${ins[1]}.bam"
 fi
 
 valgrind --leak-check=full --show-reachable=yes ../stringtie ${arrparms[$i]} -o $fout $fin
 if [ ! -f $fout ]; then
   echo "Error: file $fout not created! Failed running stringtie on $fin"
   exit 1
 fi
 if diff -q -I '^#' $fout $fcmp &>/dev/null; then
    echo "  OK."
 else
   echo "Error: test failed, output $fout different than expected ($fcmp)!"
   #exit 1
 fi
done
