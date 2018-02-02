#!/bin/bash -e
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
git checkout master
git tag -a "v$ver" -m "release $ver"
git push --tags
