#!/bin/bash -e
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
srcpack=stringtie-$ver
source prep_source.sh
linpack=$srcpack.Linux_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
cd $srcpack
make clean
make static-cpp
cp LICENSE README.md run_tests.sh stringtie prepDE.py prepDE.py3 ../$linpack/
#cp -r tests_exp_out ../$linpack/
cd ..
tar cvfz $linpack.tar.gz $linpack
ls -l $linpack.tar.gz
echo "scp $linpack.tar.gz $srcpack.tar.gz  salz:~/html/software/stringtie/dl/"
echo "If you're on CCB servers you can also update the release# links:"
echo "perl -i -pe 's/stringtie\-\d\.\d+\.\d+\w?\./stringtie-$ver./g' ~/html/software/stringtie/home.shtml"
