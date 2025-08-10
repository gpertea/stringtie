#!/usr/bin/env bash
# Prepare an offline source package including all dependencies and tests
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver.offline

set -e

echo "preparing $pack.tar.gz"
echo "-------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib

# obtain a clean htslib directory
hts=htslib
/bin/rm -rf $hts.prep
mv $hts $hts.prep
git checkout -- $hts
mv $hts $pack/
mv $hts.prep $hts

# include SuperReads_RNA if available
srm=SuperReads_RNA
if [[ -d $srm ]]; then
  mv $srm $srm.prep
  git checkout -- $srm
  mv $srm $pack/
  mv $srm.prep $srm
fi

# fetch htslib dependencies
(
  cd $pack/$hts
  if [[ ! -d libdeflate ]]; then
    git clone https://github.com/ebiggers/libdeflate
    cd libdeflate
    git checkout 7805198
    cd ..
  fi
  if [[ ! -d bzip2 ]]; then
    curl -ksLO https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
    tar -xzf bzip2-1.0.8.tar.gz
    /bin/rm -f bzip2-1.0.8.tar.gz
    mv bzip2-1.0.8 bzip2
  fi
  if [[ ! -d lzma ]]; then
    xzver="5.4.7"
    xz="xz-$xzver"
    xzurl="https://github.com/tukaani-project/xz/releases/download/v$xzver/$xz.tar.gz"
    curl -ksLO "$xzurl"
    tar -xzf $xz.tar.gz
    /bin/rm -f $xz.tar.gz
    mv $xz lzma
  fi

  # patch build_lib.sh to avoid network access
  cat > build_lib.sh <<'EOB'
#!/usr/bin/env bash
##
if [[ "$1" == "clean" ]]; then
  make clean
  /bin/rm -f config.h
  /bin/rm -rf xlibs
  exit
fi

pwd=$(pwd -P)
prefix=$pwd/xlibs
incdir=$prefix/include
libdir=$prefix/lib
mkdir -p $incdir $libdir
cc=${CC:-gcc}
cxx=${CXX:-g++}

if [[ ! -d libdeflate ]]; then
  echo "Error: libdeflate source not found!"
  exit 1
fi
if [[ ! -f $libdir/libdeflate.a ]]; then
  cd libdeflate
  MINGW=''
  libdeflate=libdeflate.a
  if [[ $($cc -dumpmachine 2>/dev/null) == *mingw* ]]; then
   MINGW=1
   libdeflate=libdeflatestatic.lib
  fi
  make -f ../Makefile.libdeflate -j 4 $libdeflate CC="$cc" CXX="$cxx" || exit 1
  cp $libdeflate $libdir/libdeflate.a
  cp libdeflate.h $incdir/
  cd ..
fi

if [[ ! -d bzip2 ]]; then
  echo "Error: bzip2 source not found!"
  exit 1
fi
if [[ ! -f $libdir/libbz2.a ]]; then
  cd bzip2
  make -j 4 libbz2.a CC="$cc"
  cp bzlib.h $incdir/
  cp libbz2.a $libdir/
  cd ..
fi

if [[ ! -d lzma ]]; then
  echo "Error: lzma source not found!"
  exit 1
fi
if [[ ! -f $libdir/liblzma.a ]]; then
  cd lzma
  CC="$cc" CXX="$cxx" ./configure --disable-shared -disable-xz -disable-xzdec --disable-lzmadec \
   --disable-lzmainfo --disable-nls --prefix=$prefix
  make -j 4 CC="$cc" CXX="$cxx"
  make install CC="$cc" CXX="$cxx"
  cd ..
fi

make -j 4 CC="$cc" CXX="$cxx" lib-static
EOB
  chmod +x build_lib.sh
)

# prepare test data
ftpack=tests_v3
ftests=$ftpack.tar.gz
tdir=tests
curl -ksLO https://github.com/gpertea/stringtie/raw/test_data/$ftests
mkdir -p $pack/$tdir
tar -xzf "$ftests" -C $pack/$tdir/ --strip-components=1
/bin/rm -f $ftests

# patch run_tests.sh to use pre-packed tests
cat > $pack/run_tests.sh <<'EOT'
#!/usr/bin/env bash
tdir=tests
if [ ! -d "$tdir" ]; then
  echo "Error: directory $tdir not found!"
  exit 1
fi
cd "$tdir"

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
   echo "Error: file $fcmp does not exist! Re-package test data."
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
EOT
chmod +x $pack/run_tests.sh

# copy general source files
cp Makefile LICENSE README.md stringtie.cpp prepDE.py prepDE.py3 $pack/ >/dev/null 2>&1 || true
cp {rlink,bundle,tablemaker,tmerge}.{h,cpp} $pack/

# copy gclib sources
cp ./gclib/{GVec,GList,khashl,GHashMap}.hh ./gclib/GBitVec.h ./gclib/xxhash.h ./gclib/wyhash.h $pack/gclib
cp ./gclib/{GArgs,GStr,GSam,GBase,gdna,codons,gff,GFaSeqGet,GFastaIndex,proc_mem,GThreads}.{h,cpp} $pack/gclib

# create archive
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
