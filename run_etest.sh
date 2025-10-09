#!/usr/bin/env bash
set -euo pipefail

d="e_test"
f="$d/missed_regions.bam"
if [[ ! -f "$f" ]]; then
  command -v curl >/dev/null 2>&1 || { echo "curl req'd" >&2; exit 1; }
  mkdir -p "$d"
  u="https://github.com/gpertea/stringtie-testdata/raw/refs/heads/master/$d"
  for x in missed_regions.bam sel_refs.gtf sel_refs.txlst; do
    echo "dl $x..."
    curl -LfsS "$u/$x" -o "$d/$x"
  done
fi

cd "$d"
../stringtie -e -o out_sel_strg302.gtf -G sel_refs.gtf missed_regions.bam
o=$(comm -23 sel_refs.txlst <(gtfcount -l out_sel_strg302.gtf | sort))
if [[ -n "$o" ]]; then
  echo ">>>>>> ERROR, missing transcripts:"
  printf "%s\n" "$o"
fi
