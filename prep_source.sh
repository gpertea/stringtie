#!/usr/bin/env bash
# this must be run into an active git repository
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver
echo "preparing $pack.tar.gz"
echo "-------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib
cd samtools-0.1.18
make clean
cd ..
# getting a clean SuperReads_RNA directory
srm=SuperReads_RNA
if [[ -d $srm ]]; then
  mv $srm $srm.prepping
  git checkout -- $srm
  mv $srm $pack/
  mv $srm.prepping $srm
fi
gldir=stringtie-$ver/gclib/
cp Makefile LICENSE README.md run_tests.sh stringtie.cpp {rlink,tablemaker,tmerge}.{h,cpp} $pack/
cp -r samtools-0.1.18 $pack/
/bin/rm -rf $pack/samtools-0.1.18/.svn
cp ./gclib/{GVec,GList,GHash,GIntHash}.hh ./gclib/GBitVec.h $gldir
cp ./gclib/{GArgs,GStr,GBam,GBase,gdna,codons,gff,GFaSeqGet,GFastaIndex,proc_mem,GThreads}.{h,cpp} $gldir
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
