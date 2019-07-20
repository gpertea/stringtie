#!/bin/sh
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver
macpack=$pack.OSX_x86_64
echo "preparing $macpack.tar.gz"
echo "-------------------"
/bin/rm -rf $macpack
/bin/rm -f $macpack.tar.gz
mkdir $macpack
make clean
make release
cp -p LICENSE README.md run_tests.sh stringtie $macpack/
tar cvfz $macpack.tar.gz $macpack
ls -l $macpack.tar.gz
#echo "If you're on igmN machines you can also update the web files:"
echo "scp $macpack.tar.gz  salz:~/html/software/stringtie/dl/"
#echo "perl -i -pe 's/stringtie-[0123]\.\d+\./stringtie-$ver./g' ~/html/software/stringtie/index.shtml"
