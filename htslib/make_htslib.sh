#!/usr/bin/env bash
## 
pwd=$(pwd -P)
prefix=$pwd/xlibs
incdir=$prefix/include
libdir=$prefix/lib
mkdir -p $incdir
mkdir -p $libdir

# -- prepare libdeflate
git clone https://github.com/ebiggers/libdeflate
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

# -- prepare liblzma
xz=xz-5.2.5
curl -sLO https://tukaani.org/xz/$xz.tar.gz
tar -xzf $xz.tar.gz
mv $xz lzma
cd lzma
./configure --disable-shared -disable-xz -disable-xzdec --disable-lzmadec \
 --disable-lzmainfo --disable-nls --prefix=$prefix
make -j 4
make install

cd ..

## now build static library for linking with stringtie, htsqc etc.
make -j 4 lib-static
