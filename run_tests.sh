#!/usr/bin/env bash
#set -euo pipefail

tests_tag=v3.0.2 ## update as needed
ftests=${tests_tag}.tar.gz
tdir=tests

function unpack_test_data() {
  if [[ ! -f $ftests ]]; then
    echo "Error: file $ftests not found!"
    exit 1
  fi
  echo "..unpacking test data.."
  echo
  mkdir -p $tdir
  tar -xzf "$ftests" -C $tdir/ --strip-components=2
  if [ ! -f $tdir/human-chr19_P.gff ]; then
     echo "Error: invalid test data archive?"
     exit 1
  fi
  /bin/rm -f $ftests ## comment this out for testing
}

if [ -f $ftests ]; then
  #extract the tarball and rename the directory
  echo "..Using existing $ftests"
  unpack_test_data
else
  echo "..Downloading test data.."
  #use curl to fetch the tarball from a specific github release or branch
  curl -ksLO https://github.com/gpertea/stringtie-testdata/archive/refs/tags/${tests_tag}.tar.gz
  unpack_test_data
fi

# check a gene abundance (-A) table: a gene ID must be listed only once, and
# all $2 genes of the guide annotation must be there, with or without reads
function check_gene_abund() {
  local ndup=$(tail -n +2 $1 | cut -f1 | sort | uniq -d | wc -l)
  local nref=$(tail -n +2 $1 | cut -f1 | grep -v '^STRG\.' | sort -u | wc -l)
  if [[ $ndup -gt 0 ]]; then
    echo "Error: test failed, $ndup gene IDs are listed more than once in $1!"
    return 1
  fi
  if [[ $nref -ne $2 ]]; then
    echo "Error: test failed, $nref guide genes found in $1, expected $2!"
    return 1
  fi
  return 0
}

# check that a GTF has a transcript with the start, end, strand and exon count given in $2
function check_tx() {
  awk -F'\t' '$3=="transcript" {if (t!="") print t" "n; t=$4" "$5" "$7; n=0}
     $3=="exon" {n++} END {if (t!="") print t" "n}' $1 | grep -q -x "$2"
}

cd $tdir

# point features for the --ptf tests; the - strand has its TSS at the right end
printf "chr16\t789250\t+\tTSS\nchr16\t798000\t+\tCPAS\n" > short_reads.ptf
printf "chr19\t567231\t-\tCPAS\nchr19\t571726\t-\tTSS\n" > long_reads.ptf

arrins=("short_reads" "short_reads_and_superreads" "mix_short" "long_reads" \
   "long_reads" "mix_short mix_long" "mix_short mix_long" "mix_short" "mix_short" \
   "short_reads" "long_reads")
arrparms=("" "" "-G mix_guides.gff" "-L" "-L -G human-chr19_P.gff" "--mix" "--mix -G mix_guides.gff" \
 "-N -G mix_guides.gff" "--nasc -G mix_guides.gff" \
 "--ptf short_reads.ptf" "-L --ptf long_reads.ptf")
# number of guide genes for the tests also checking the gene abundance (-A) output
arrgenes=("" "" 1 "" 14112 "" 1 1 1 "" "")
# for the --ptf tests: start, end, strand and exon count of a transcript that must be found,
# with both ends at the given points and all its introns, instead of an expected output file
arrtx=("" "" "" "" "" "" "" "" "" "789250 798000 + 20" "567231 571726 - 3")
arrout=("short_reads" "short_reads_and_superreads" "short_guided" "long_reads" \
       "long_reads_guided" "mix_reads" "mix_reads_guided" "mix_short_N_guided" "mix_short_nasc_guided" \
       "short_reads_ptf" "long_reads_ptf")
arrmsg=("Short reads"  "Short reads and super-reads" "Short reads with annotation guides" "Long reads" \
 "Long reads with annotation guides" "Mixed reads" "Mixed reads with annotation guides" \
 "Short reads with -N" "Short reads with --nasc" \
 "Short reads with point features" "Long reads with point features")

nfailed=0
for i in ${!arrmsg[@]}; do
 fout="${arrout[$i]}.out.gtf"
 /bin/rm -f $fout
 fcmp="${arrout[$i]}.out_expected.gtf"
 if [[ -n "${arrtx[$i]}" ]]; then
   fcmp="a transcript with start, end, strand, exons: ${arrtx[$i]}"
 elif [ ! -f $fcmp ]; then
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
 parms="${arrparms[$i]}"
 fabund=""
 if [[ -n "${arrgenes[$i]}" ]]; then
   fabund="${arrout[$i]}.out.abund.tab"
   /bin/rm -f $fabund
   parms="$parms -A $fabund"
 fi
 echo "Running: ../stringtie $parms -o $fout $fin"
 ../stringtie $parms -o $fout $fin
 if [ ! -f $fout ]; then
   echo "Error: file $fout not created! Failed running stringtie on $fin"
   exit 1
 fi
 abund_ok=1
 if [[ -n "$fabund" ]] && ! check_gene_abund $fabund ${arrgenes[$i]}; then
   abund_ok=0
   nfailed=$((nfailed+1))
 fi
 if [[ -n "${arrtx[$i]}" ]]; then
   check_tx $fout "${arrtx[$i]}"
 else
   diff -q -I '^#' $fout $fcmp &>/dev/null
 fi
 if [[ $? -eq 0 ]]; then
    [[ $abund_ok -eq 1 ]] && echo "  OK."
 else
   echo "Error: test failed, output $fout different than expected ($fcmp)!"
   nfailed=$((nfailed+1))
 fi
done
if [[ $nfailed -gt 0 ]]; then
  echo "Error: $nfailed test(s) failed."
  exit 1
fi
