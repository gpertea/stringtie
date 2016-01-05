#!/bin/sh
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver
linpack=$pack.Linux_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
make cleanall
make release
cp -p LICENSE README stringtie $linpack/
tar cvfz $linpack.tar.gz $linpack
ls -l $linpack.tar.gz
echo "If you're on igmN machines you can also update the web files:"
echo "cp $linpack.tar.gz $pack.tar.gz  ~/html/software/stringtie/dl/"
echo "perl -i -pe 's/stringtie\-\d\.\d+\.\d+\./stringtie-$ver./g' ~/html/software/stringtie/home.shtml"
