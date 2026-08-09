#!/bin/sh
set -eu

bin=${STRINGTIE_BIN:-./stringtie}
tmpdir=$(mktemp -d)
cleanup() {
  if test "${KEEP_TMP:-0}" = 1; then
    printf '%s\n' "merge policy files: $tmpdir"
  else
    rm -rf "$tmpdir"
  fi
}
trap cleanup EXIT HUP INT TERM

run_stringtie() {
  if test "${VALGRIND:-0}" = 1; then
    valgrind --quiet --error-exitcode=99 "$bin" "$@"
  else
    "$bin" "$@"
  fi
}

has_tx() {
  awk -v ref="$2" -v start="$3" -v end="$4" -v strand="$5" \
    '$1==ref && $3=="transcript" && $4==start && $5==end && $7==strand {found=1} END {exit !found}' "$1"
}

printf '%s\n' \
'chr1\tref\ttranscript\t100\t399\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";' \
'chr1\tref\texon\t100\t199\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";' \
'chr1\tref\texon\t300\t399\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";' \
'chr1\tref\tCDS\t120\t199\t.\t+\t0\tgene_id "geneA"; transcript_id "txA";' \
'chr1\tref\tCDS\t300\t350\t.\t+\t0\tgene_id "geneA"; transcript_id "txA";' \
'chr1\tref\ttranscript\t100\t399\t.\t+\t.\tgene_id "geneA"; transcript_id "txA2";' \
'chr1\tref\texon\t100\t199\t.\t+\t.\tgene_id "geneA"; transcript_id "txA2";' \
'chr1\tref\texon\t300\t399\t.\t+\t.\tgene_id "geneA"; transcript_id "txA2";' \
'chr1\tref\ttranscript\t600\t899\t.\t+\t.\tgene_id "geneB"; transcript_id "txB";' \
'chr1\tref\texon\t600\t699\t.\t+\t.\tgene_id "geneB"; transcript_id "txB";' \
'chr1\tref\texon\t800\t899\t.\t+\t.\tgene_id "geneB"; transcript_id "txB";' \
'chr1\tref\tCDS\t620\t699\t.\t+\t0\tgene_id "geneB"; transcript_id "txB";' \
'chr1\tref\tCDS\t800\t850\t.\t+\t0\tgene_id "geneB"; transcript_id "txB";' \
| sed 's/\\t/\t/g' > "$tmpdir/guide.gtf"

printf '%s\n' \
'chr1\tin\ttranscript\t80\t450\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext1"; FPKM "0.1";' \
'chr1\tin\texon\t80\t199\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext1";' \
'chr1\tin\texon\t300\t450\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext1";' \
'chr1\tin\ttranscript\t100\t899\t.\t+\t.\tgene_id "bridge"; transcript_id "bridge1"; FPKM "2";' \
'chr1\tin\texon\t100\t199\t.\t+\t.\tgene_id "bridge"; transcript_id "bridge1";' \
'chr1\tin\texon\t300\t399\t.\t+\t.\tgene_id "bridge"; transcript_id "bridge1";' \
'chr1\tin\texon\t600\t699\t.\t+\t.\tgene_id "bridge"; transcript_id "bridge1";' \
'chr1\tin\texon\t800\t899\t.\t+\t.\tgene_id "bridge"; transcript_id "bridge1";' \
| sed 's/\\t/\t/g' > "$tmpdir/in1.gtf"

printf '%s\n' \
'chr1\tin\ttranscript\t80\t450\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext2"; FPKM "0.1";' \
'chr1\tin\texon\t80\t199\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext2";' \
'chr1\tin\texon\t300\t450\t.\t+\t.\tgene_id "sampleA"; transcript_id "ext2";' \
| sed 's/\\t/\t/g' > "$tmpdir/in2.gtf"

run_stringtie --merge -G "$tmpdir/guide.gtf" -F 1 --report "$tmpdir/report.tsv" \
  -o "$tmpdir/out.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"

test "$(grep -c 'transcript_id "txA"' "$tmpdir/out.gtf")" -eq 3
grep -Eq 'transcript[[:space:]]+80[[:space:]]+450.*reference_id "txA"' "$tmpdir/out.gtf"
grep -Eq 'reject[[:space:]]+cross_gene_bridge.*bridge1' "$tmpdir/report.tsv"
grep -Eq 'rescue[[:space:]]+replicated_rescue[[:space:]]+ext1' "$tmpdir/report.tsv"

ids=$(awk '$3=="transcript" {for(i=9;i<=NF;i++) if($i=="transcript_id") {gsub(/[";]/,"",$(i+1)); print $(i+1)}}' "$tmpdir/out.gtf")
test "$(printf '%s\n' "$ids" | sort | uniq -d | wc -l)" -eq 0

run_stringtie --merge --dedup-guides -G "$tmpdir/guide.gtf" -o "$tmpdir/dedup.gtf" \
  "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
test "$(grep -c 'transcript_id "txA2"' "$tmpdir/dedup.gtf")" -eq 0
grep -q 'reference_aliases "txA2"' "$tmpdir/dedup.gtf"

run_stringtie --merge --guide-policy reference -G "$tmpdir/guide.gtf" -F 1 \
  -o "$tmpdir/reference.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
test "$(awk '$3=="transcript" && $4==80 && $5==450 {n++} END{print n+0}' "$tmpdir/reference.gtf")" -eq 0

run_stringtie --merge --guide-policy best --dedup-guides -G "$tmpdir/guide.gtf" -F 1 \
  -o "$tmpdir/best.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
grep -Eq 'transcript[[:space:]]+80[[:space:]]+450.*reference_id "txA"' "$tmpdir/best.gtf"
test "$(awk '$3=="transcript" && /transcript_id "txA"/ {n++} END{print n+0}' "$tmpdir/best.gtf")" -eq 0

run_stringtie --merge --path-mode replicated --min-piece-support 2 -G "$tmpdir/guide.gtf" \
  -o "$tmpdir/replicated.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
grep -Eq 'transcript[[:space:]]+80[[:space:]]+450' "$tmpdir/replicated.gtf"

run_stringtie --merge --allow-fusions -G "$tmpdir/guide.gtf" --report "$tmpdir/fusion.tsv" \
  -o "$tmpdir/fusion.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
grep -q 'fusion_allowed' "$tmpdir/fusion.tsv"

run_stringtie --merge --min-introns 0 -G "$tmpdir/guide.gtf" --report "$tmpdir/fusion-min0.tsv" \
  -o "$tmpdir/fusion-min0.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
grep -Eq 'reject[[:space:]]+cross_gene_bridge.*bridge1' "$tmpdir/fusion-min0.tsv"

## verify that terminal SET overlap cannot bypass the default intron gate
printf '%s\n' \
'chrS\tin\ttranscript\t5000\t5299\t.\t+\t.\tgene_id "geneS"; transcript_id "metS";' \
'chrS\tin\texon\t5000\t5099\t.\t+\t.\tgene_id "geneS"; transcript_id "metS";' \
'chrS\tin\texon\t5200\t5299\t.\t+\t.\tgene_id "geneS"; transcript_id "metS";' \
'chrS\tin\ttranscript\t5250\t5450\t.\t+\t.\tgene_id "geneS"; transcript_id "setS";' \
'chrS\tin\texon\t5250\t5450\t.\t+\t.\tgene_id "geneS"; transcript_id "setS";' \
| sed 's/\\t/\t/g' > "$tmpdir/terminal.gtf"

run_stringtie --merge -o "$tmpdir/terminal-default.gtf" "$tmpdir/terminal.gtf"
has_tx "$tmpdir/terminal-default.gtf" chrS 5000 5299 +
has_tx "$tmpdir/terminal-default.gtf" chrS 5250 5450 +
! has_tx "$tmpdir/terminal-default.gtf" chrS 5000 5450 +

run_stringtie --merge --min-introns 0 -o "$tmpdir/terminal-relaxed.gtf" "$tmpdir/terminal.gtf"
has_tx "$tmpdir/terminal-relaxed.gtf" chrS 5000 5450 +

## verify independent 3-prime support and hard-cap controls on a plus-strand join
printf '%s\n' \
'chrE\tin\ttranscript\t6000\t6499\t.\t+\t.\tgene_id "geneE"; transcript_id "left1";' \
'chrE\tin\texon\t6000\t6099\t.\t+\t.\tgene_id "geneE"; transcript_id "left1";' \
'chrE\tin\texon\t6200\t6299\t.\t+\t.\tgene_id "geneE"; transcript_id "left1";' \
'chrE\tin\texon\t6400\t6499\t.\t+\t.\tgene_id "geneE"; transcript_id "left1";' \
'chrE\tin\ttranscript\t6200\t6699\t.\t+\t.\tgene_id "geneE"; transcript_id "right1";' \
'chrE\tin\texon\t6200\t6299\t.\t+\t.\tgene_id "geneE"; transcript_id "right1";' \
'chrE\tin\texon\t6400\t6499\t.\t+\t.\tgene_id "geneE"; transcript_id "right1";' \
'chrE\tin\texon\t6600\t6699\t.\t+\t.\tgene_id "geneE"; transcript_id "right1";' \
| sed 's/\\t/\t/g' > "$tmpdir/end1.gtf"

printf '%s\n' \
'chrE\tin\ttranscript\t6000\t6499\t.\t+\t.\tgene_id "geneE"; transcript_id "left2";' \
'chrE\tin\texon\t6000\t6099\t.\t+\t.\tgene_id "geneE"; transcript_id "left2";' \
'chrE\tin\texon\t6200\t6299\t.\t+\t.\tgene_id "geneE"; transcript_id "left2";' \
'chrE\tin\texon\t6400\t6499\t.\t+\t.\tgene_id "geneE"; transcript_id "left2";' \
| sed 's/\\t/\t/g' > "$tmpdir/end2.gtf"

run_stringtie --merge -o "$tmpdir/end-default.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
! has_tx "$tmpdir/end-default.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-3p 1 -o "$tmpdir/end-3p.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
has_tx "$tmpdir/end-3p.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-5p 1 -o "$tmpdir/end-5p.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
! has_tx "$tmpdir/end-5p.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-3p 1 --max-extension-3p 100 \
  -o "$tmpdir/end-cap.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
! has_tx "$tmpdir/end-cap.gtf" chrE 6000 6699 +

## verify strand-aware guide extensions and flexible reference-end policies
printf '%s\n' \
'chrP\tref\ttranscript\t1000\t1299\t.\t+\t.\tgene_id "geneP"; transcript_id "txP";' \
'chrP\tref\texon\t1000\t1099\t.\t+\t.\tgene_id "geneP"; transcript_id "txP";' \
'chrP\tref\texon\t1200\t1299\t.\t+\t.\tgene_id "geneP"; transcript_id "txP";' \
'chrM\tref\ttranscript\t2000\t2299\t.\t-\t.\tgene_id "geneM"; transcript_id "txM";' \
'chrM\tref\texon\t2000\t2099\t.\t-\t.\tgene_id "geneM"; transcript_id "txM";' \
'chrM\tref\texon\t2200\t2299\t.\t-\t.\tgene_id "geneM"; transcript_id "txM";' \
| sed 's/\\t/\t/g' > "$tmpdir/end-guide.gtf"

printf '%s\n' \
'chrP\tin\ttranscript\t1000\t1499\t.\t+\t.\tgene_id "sampleP"; transcript_id "extP";' \
'chrP\tin\texon\t1000\t1099\t.\t+\t.\tgene_id "sampleP"; transcript_id "extP";' \
'chrP\tin\texon\t1200\t1499\t.\t+\t.\tgene_id "sampleP"; transcript_id "extP";' \
'chrM\tin\ttranscript\t2000\t2499\t.\t-\t.\tgene_id "sampleM"; transcript_id "extM";' \
'chrM\tin\texon\t2000\t2099\t.\t-\t.\tgene_id "sampleM"; transcript_id "extM";' \
'chrM\tin\texon\t2200\t2499\t.\t-\t.\tgene_id "sampleM"; transcript_id "extM";' \
| sed 's/\\t/\t/g' > "$tmpdir/end-input.gtf"

run_stringtie --merge -G "$tmpdir/end-guide.gtf" --max-extension-3p 50 \
  -o "$tmpdir/guide-cap-3p.gtf" "$tmpdir/end-input.gtf"
! has_tx "$tmpdir/guide-cap-3p.gtf" chrP 1000 1499 +
has_tx "$tmpdir/guide-cap-3p.gtf" chrM 2000 2499 -
run_stringtie --merge -G "$tmpdir/end-guide.gtf" --max-extension-5p 50 \
  -o "$tmpdir/guide-cap-5p.gtf" "$tmpdir/end-input.gtf"
has_tx "$tmpdir/guide-cap-5p.gtf" chrP 1000 1499 +
! has_tx "$tmpdir/guide-cap-5p.gtf" chrM 2000 2499 -
run_stringtie --merge -G "$tmpdir/end-guide.gtf" --max-extension-3p 250 \
  -o "$tmpdir/guide-relaxed.gtf" "$tmpdir/end-input.gtf"
has_tx "$tmpdir/guide-relaxed.gtf" chrP 1000 1499 +
run_stringtie --merge -G "$tmpdir/end-guide.gtf" --guide-policy reference \
  -o "$tmpdir/guide-reference.gtf" "$tmpdir/end-input.gtf"
! has_tx "$tmpdir/guide-reference.gtf" chrP 1000 1499 +
! has_tx "$tmpdir/guide-reference.gtf" chrM 2000 2499 -

## verify same-gene coding bridges require explicit authorization
printf '%s\n' \
'chrC\tref\ttranscript\t100\t399\t.\t+\t.\tgene_id "geneC"; transcript_id "txC1";' \
'chrC\tref\texon\t100\t199\t.\t+\t.\tgene_id "geneC"; transcript_id "txC1";' \
'chrC\tref\texon\t300\t399\t.\t+\t.\tgene_id "geneC"; transcript_id "txC1";' \
'chrC\tref\tCDS\t120\t199\t.\t+\t0\tgene_id "geneC"; transcript_id "txC1";' \
'chrC\tref\tCDS\t300\t350\t.\t+\t0\tgene_id "geneC"; transcript_id "txC1";' \
'chrC\tref\ttranscript\t300\t599\t.\t+\t.\tgene_id "geneC"; transcript_id "txC2";' \
'chrC\tref\texon\t300\t399\t.\t+\t.\tgene_id "geneC"; transcript_id "txC2";' \
'chrC\tref\texon\t500\t599\t.\t+\t.\tgene_id "geneC"; transcript_id "txC2";' \
'chrC\tref\tCDS\t320\t399\t.\t+\t0\tgene_id "geneC"; transcript_id "txC2";' \
'chrC\tref\tCDS\t500\t550\t.\t+\t0\tgene_id "geneC"; transcript_id "txC2";' \
| sed 's/\\t/\t/g' > "$tmpdir/cds-guide.gtf"

printf '%s\n' \
'chrC\tin\ttranscript\t100\t599\t.\t+\t.\tgene_id "sampleC"; transcript_id "bridgeC";' \
'chrC\tin\texon\t100\t199\t.\t+\t.\tgene_id "sampleC"; transcript_id "bridgeC";' \
'chrC\tin\texon\t300\t399\t.\t+\t.\tgene_id "sampleC"; transcript_id "bridgeC";' \
'chrC\tin\texon\t500\t599\t.\t+\t.\tgene_id "sampleC"; transcript_id "bridgeC";' \
| sed 's/\\t/\t/g' > "$tmpdir/cds-input.gtf"

run_stringtie --merge -G "$tmpdir/cds-guide.gtf" --report "$tmpdir/cds-default.tsv" \
  -o "$tmpdir/cds-default.gtf" "$tmpdir/cds-input.gtf"
! has_tx "$tmpdir/cds-default.gtf" chrC 100 599 +
grep -q 'coding_transcript_bridge' "$tmpdir/cds-default.tsv"
run_stringtie --merge -G "$tmpdir/cds-guide.gtf" --allow-cds-bridges --report "$tmpdir/cds-allow.tsv" \
  -o "$tmpdir/cds-allow.gtf" "$tmpdir/cds-input.gtf"
has_tx "$tmpdir/cds-allow.gtf" chrC 100 599 +
grep -q 'cds_bridge_allowed' "$tmpdir/cds-allow.tsv"

printf '%s\n' 'merge policy test: OK'
