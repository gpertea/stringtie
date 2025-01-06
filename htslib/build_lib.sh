#!/usr/bin/env bash
## 

if [[ "$1" == "clean" ]]; then
  make clean
  /bin/rm -f config.h
  /bin/rm -rf xlibs
  /bin/rm -rf lzma bzip2 libdeflate
  /bin/rm -f *.tar.gz
  exit
fi

pwd=$(pwd -P)
prefix=$pwd/xlibs
incdir=$prefix/include
libdir=$prefix/lib
mkdir -p $incdir
mkdir -p $libdir

# -- prepare libdeflate
if [[ ! -d libdeflate ]]; then
 git clone https://github.com/ebiggers/libdeflate
 cd libdeflate
 #git checkout '9b565afd996d8b798fc7b94cddcc7cfa49293050'
 git checkout 7805198 # release v1.23
 cd ..
fi
if [[ ! -f $libdir/libdeflate.a ]]; then
  cd libdeflate
  MINGW=''
  libdeflate=libdeflate.a
  if [[ $(gcc -dumpmachine) == *mingw* ]]; then
   MINGW=1
   libdeflate=libdeflatestatic.lib
  fi
  make -f ../Makefile.libdeflate -j 4 $libdeflate || exit 1
  cp $libdeflate $libdir/libdeflate.a
  cp libdeflate.h $incdir/
  cd ..
fi

bzip="bzip2-1.0.8"
if [[ ! -d bzip2 ]]; then
  curl -ksLO https://sourceware.org/pub/bzip2/$bzip.tar.gz || \
    exec echo "Error: failed to fetch $bzip.tar.gz!"
  tar -xzf $bzip.tar.gz || exec echo "Error: failed to unpack $bzip.tar.gz!"
  /bin/rm -f $bzip.tar.gz
  mv $bzip bzip2
fi
if [[ ! -f $libdir/libbz2.a ]]; then 
  cd bzip2
  make -j 4 libbz2.a
  cp bzlib.h $incdir/
  cp libbz2.a $libdir/
  cd ..
fi

# -- prepare liblzma
xzver="5.4.7"
xz="xz-$xzver"
xzurl="https://github.com/tukaani-project/xz/releases/download/v$xzver/$xz.tar.gz"
if [[ ! -d lzma ]]; then
  curl -ksLO "$xzurl" || exec echo "Error: failed to fetch: $xzurl!"
  tar -xzf $xz.tar.gz || exec echo "Error: failed to unpack $xz.tar.gz !"
  unlink $xz.tar.gz
  mv $xz lzma
fi
if [[ ! -f $libdir/liblzma.a ]]; then
  cd lzma
  ./configure --disable-shared -disable-xz -disable-xzdec --disable-lzmadec \
   --disable-lzmainfo --disable-nls --prefix=$prefix
  make -j 4
  make install
  cd ..
fi

make -j 4 lib-static
