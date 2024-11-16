#!/usr/bin/env bash
tdir=tests_3
tests=$tdir.tar.gz
function unpack_test_data() {
  if [ ! -f $tests ]; then
    echo "Error: file $tests not found!"
    exit 1
  fi
  echo "..unpacking test data.."
  echo
  tar -xzf $tests
  if [ ! -f $tdir/human-chr19_P.gff ]; then
     echo "Error: invalid test data archive?"
     exit 1
  fi
  /bin/rm -f $tests ## comment this out for testing
}

if [ -f $tests ]; then
    #extract the tarball and rename the directory
    echo "..Using existing $tests"
    unpack_test_data
  else
    echo "..Downloading test data.."
   #use curl to fetch the tarball from a specific github release or branch
    curl -ksLO https://github.com/gpertea/stringtie/raw/test_data/$tests
    unpack_test_data
fi
cd $tdir

arrins=("short_reads" "short_reads_and_superreads" "mix_short" "long_reads" \
   "long_reads" "mix_short mix_long" "mix_short mix_long" "mix_short" "mix_short")
arrparms=("" "" "-G mix_guides.gff" "-L" "-L -G human-chr19_P.gff" "--mix" "--mix -G mix_guides.gff" \
 "-N -G mix_guides.gff" "--nasc -G mix_guides.gff")
arrout=("short_reads" "short_reads_and_superreads" "short_guided" "long_reads" \
       "long_reads_guided" "mix_reads" "mix_reads_guided" "mix_short_N_guided" "mix_short_nasc_guided")
arrmsg=("Short reads"  "Short reads and super-reads" "Short reads with annotation guides" "Long reads" \
 "Long reads with annotation guides" "Mixed reads" "Mixed reads with annotation guides" \
 "Short reads with -N" "Short reads with --nasc")

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
   if [[ ${#ins[@]} -gt 1 ]]; then
     fin="${ins[0]}.bam ${ins[1]}.bam"
   fi
 fi
 echo "Running: ../stringtie ${arrparms[$i]} -o $fout $fin"
 ../stringtie ${arrparms[$i]} -o $fout $fin
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
