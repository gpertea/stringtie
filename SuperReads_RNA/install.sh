#!/bin/sh
set -e
ROOT=`pwd -P`
[ -z "$DEST" ] && DEST="$ROOT"

###################
# Check for gmake #
###################
mkdir -p dist-bin
PATH=$PATH:$ROOT/dist-bin
ln -sf $(which make) $ROOT/dist-bin/gmake
ln -sf $ROOT/PkgConfig.pm $ROOT/dist-bin/pkg-config

export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
BINDIR=$DEST/bin
LIBDIR=$DEST/lib
export PKG_CONFIG_PATH=$LIBDIR/pkgconfig:$PKG_CONFIG_PATH
cd global-1
autoreconf -fis
./configure --prefix=$DEST --bindir=$BINDIR --libdir=$LIBDIR && make -j $NUM_THREADS install-special
cd ..
perl -pe 's{^BIN_DIR = #__#}{BIN_DIR = "'$BINDIR'"}' global-1/SuperReadsR/create_rna_sr.py > $DEST/create_rna_sr.py

chmod 755 $DEST/create_rna_sr.py
echo "creating sr_config_example.txt with correct PATHs"
$BINDIR/createSuperReads_RNA -g sr_config_example.txt
echo "You can now copy the script create_rna_sr.py in any directory of PATH for convenience"
echo "or run the script with its full path: $DEST/create_rna_sr.py"
