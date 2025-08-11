#!/usr/bin/env bash
# this must be run into an active git repository
#url=https://github.com/gpertea/stringtie.git
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver
echo "preparing $pack.tar.gz"
echo "-------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir -p $pack/gclib
mkdir -p $pack/tests
#cd htslib
#make clean
#/bin/rm -f xlibs config.h lzma bzip2 libdeflate
#cd ..
# getting a clean htslib directory
hts=htslib
/bin/rm -rf $hts.prep
mv $hts $hts.prep
git checkout -- $hts
mv $hts $pack/
mv $hts.prep $hts
# getting a clean SuperReads_RNA directory
srm=SuperReads_RNA
if [[ -d $srm ]]; then
  mv $srm $srm.prep
  git checkout -- $srm
  mv $srm $pack/
  mv $srm.prep $srm
fi
gldir=stringtie-$ver/gclib/
cp Makefile LICENSE README.md run_tests.sh stringtie.cpp prepDE.py prepDE.py3 {rlink,bundle,tablemaker,tmerge}.{h,cpp} $pack/
cp ./gclib/{GVec,GList,khashl,GHashMap}.hh ./gclib/GBitVec.h ./gclib/xxhash.h ./gclib/wyhash.h $gldir
cp ./gclib/{GArgs,GStr,GSam,GBase,gdna,codons,gff,GFaSeqGet,GFastaIndex,proc_mem,GThreads}.{h,cpp} $gldir
cp ./tests/README.md $pack/tests/
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
