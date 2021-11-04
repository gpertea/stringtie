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
 git checkout '9b565afd996d8b798fc7b94cddcc7cfa49293050'
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
  make -j 2 $libdeflate || exit 1
  cp $libdeflate $libdir/libdeflate.a
  cp libdeflate.h $incdir/
  cd ..
fi


bzip="bzip2-1.0.8"
if [[ ! -d bzip2 ]]; then
  curl -sLO https://sourceware.org/pub/bzip2/$bzip.tar.gz
  tar -xzf $bzip.tar.gz
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
xz="xz-5.2.5"

if [[ ! -d lzma ]]; then
  curl -sLO https://tukaani.org/xz/$xz.tar.gz
  tar -xzf $xz.tar.gz
  /bin/rm -f $xz.tar.gz
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
