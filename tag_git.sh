#!/bin/bash -e
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
if [[ "$1" == "delete" || "$1" == "del" ]]; then
  echo "Deleting tag v$ver .."
  git tag -d v$ver
  git push origin :refs/tags/v$ver
  exit
fi
echo "Tagging with v$ver"
git checkout master
git tag -a "v$ver" -m "release $ver"
git push --tags
