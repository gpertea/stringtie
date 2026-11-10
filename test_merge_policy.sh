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

## negative assertion; "! has_tx ..." would be ignored by set -e
no_tx() {
  if has_tx "$@"; then
    printf '%s\n' "unexpected transcript: $*" >&2
    exit 1
  fi
}

## print "exon chain<TAB>attributes" for every transcript on ref/strand
tx_chains() {
  awk -F'\t' -v ref="$2" -v strand="$3" '
    $1==ref && $7==strand && $3=="transcript" {
      if (id!="") print chain "\t" attrs
      id=$9; sub(/.*transcript_id "/,"",id); sub(/".*/,"",id); attrs=$9; chain=""
    }
    $1==ref && $7==strand && $3=="exon" && id!="" {
      chain=chain (chain=="" ? "" : ",") $4 "-" $5
    }
    END { if (id!="") print chain "\t" attrs }' "$1"
}

## has_chain FILE REF STRAND EXONS [ATTR VALUE]: complete exon chain, optional attribute value
has_chain() {
  tx_chains "$1" "$2" "$3" | awk -F'\t' -v chain="$4" -v attr="${5:-}" -v val="${6:-}" '
    $1==chain && (attr=="" || index($2, attr " \"" val "\";")) {found=1}
    END {exit !found}'
}

need_chain() {
  if ! has_chain "$@"; then
    printf '%s\n' "missing transcript: $*" >&2
    exit 1
  fi
}

no_chain() {
  if has_chain "$@"; then
    printf '%s\n' "unexpected transcript: $*" >&2
    exit 1
  fi
}

## tx_count FILE REF: number of transcripts on REF
tx_count() {
  awk -F'\t' -v ref="$2" '$1==ref && $3=="transcript" {n++} END {print n+0}' "$1"
}

## known_fail NAME CMD...: a documented defect. Passes while CMD fails; when a fix
## makes CMD succeed the suite stops so the check is promoted to an assertion.
known_fail() {
  kf_name=$1; shift
  if ( "$@" ) >/dev/null 2>&1; then
    printf '%s\n' "known failure now passes, promote it to an assertion: $kf_name" >&2
    exit 1
  fi
  known_failures="${known_failures:-} $kf_name"
}

## mkgtf FILE SRC LINE...: LINE is "ref strand tx_id s-e,s-e,..."
mkgtf() {
  mk_out=$1; mk_src=$2; shift 2
  for mk_line in "$@"; do
    printf '%s\n' "$mk_line"
  done | awk -v src="$mk_src" '{
    n=split($4,e,","); split(e[1],a,"-"); split(e[n],b,"-")
    at="gene_id \"" $3 "_g\"; transcript_id \"" $3 "\";"
    printf "%s\t%s\ttranscript\t%s\t%s\t.\t%s\t.\t%s\n",$1,src,a[1],b[2],$2,at
    for(i=1;i<=n;i++){ split(e[i],x,"-"); printf "%s\t%s\texon\t%s\t%s\t.\t%s\t.\t%s\n",$1,src,x[1],x[2],$2,at }
  }' > "$mk_out"
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

## an input with a guide's intron chain is a guide variant only beyond the guide slop
GS0="--guide-slop-bp-5p 0 --guide-slop-bp-3p 0"
run_stringtie --merge -G "$tmpdir/guide.gtf" -F 1 --report "$tmpdir/slop.tsv" \
  -o "$tmpdir/slop.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
no_tx "$tmpdir/slop.gtf" chr1 80 450 +
grep -Eq 'transcript_id "txA"; .*merge_support "2"' "$tmpdir/slop.gtf"

run_stringtie --merge $GS0 -G "$tmpdir/guide.gtf" -F 1 --report "$tmpdir/report.tsv" \
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

run_stringtie --merge $GS0 --guide-policy best --dedup-guides -G "$tmpdir/guide.gtf" -F 1 \
  -o "$tmpdir/best.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf"
grep -Eq 'transcript[[:space:]]+80[[:space:]]+450.*reference_id "txA"' "$tmpdir/best.gtf"
test "$(awk '$3=="transcript" && /transcript_id "txA"/ {n++} END{print n+0}' "$tmpdir/best.gtf")" -eq 0

run_stringtie --merge $GS0 --path-mode replicated --min-piece-support 2 -G "$tmpdir/guide.gtf" \
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
no_tx "$tmpdir/terminal-default.gtf" chrS 5000 5450 +

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
no_tx "$tmpdir/end-default.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-3p 1 -o "$tmpdir/end-3p.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
has_tx "$tmpdir/end-3p.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-5p 1 -o "$tmpdir/end-5p.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
no_tx "$tmpdir/end-5p.gtf" chrE 6000 6699 +
run_stringtie --merge --min-support-3p 1 --max-extension-3p 100 \
  -o "$tmpdir/end-cap.gtf" "$tmpdir/end1.gtf" "$tmpdir/end2.gtf"
no_tx "$tmpdir/end-cap.gtf" chrE 6000 6699 +

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

ES1="--min-support-5p 1 --min-support-3p 1"
## a single-file extension of a guide end is an unsupported overhang by default
run_stringtie --merge $GS0 -G "$tmpdir/end-guide.gtf" -o "$tmpdir/guide-unsup.gtf" "$tmpdir/end-input.gtf"
no_tx "$tmpdir/guide-unsup.gtf" chrP 1000 1499 +
no_tx "$tmpdir/guide-unsup.gtf" chrM 2000 2499 -
run_stringtie --merge $GS0 $ES1 -G "$tmpdir/end-guide.gtf" --max-extension-3p 50 \
  -o "$tmpdir/guide-cap-3p.gtf" "$tmpdir/end-input.gtf"
no_tx "$tmpdir/guide-cap-3p.gtf" chrP 1000 1499 +
has_tx "$tmpdir/guide-cap-3p.gtf" chrM 2000 2499 -
run_stringtie --merge $GS0 $ES1 -G "$tmpdir/end-guide.gtf" --max-extension-5p 50 \
  -o "$tmpdir/guide-cap-5p.gtf" "$tmpdir/end-input.gtf"
has_tx "$tmpdir/guide-cap-5p.gtf" chrP 1000 1499 +
no_tx "$tmpdir/guide-cap-5p.gtf" chrM 2000 2499 -
run_stringtie --merge $GS0 $ES1 -G "$tmpdir/end-guide.gtf" --max-extension-3p 250 \
  -o "$tmpdir/guide-relaxed.gtf" "$tmpdir/end-input.gtf"
has_tx "$tmpdir/guide-relaxed.gtf" chrP 1000 1499 +
run_stringtie --merge -G "$tmpdir/end-guide.gtf" --guide-policy reference \
  -o "$tmpdir/guide-reference.gtf" "$tmpdir/end-input.gtf"
no_tx "$tmpdir/guide-reference.gtf" chrP 1000 1499 +
no_tx "$tmpdir/guide-reference.gtf" chrM 2000 2499 -

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
no_tx "$tmpdir/cds-default.gtf" chrC 100 599 +
grep -q 'coding_transcript_bridge' "$tmpdir/cds-default.tsv"
run_stringtie --merge -G "$tmpdir/cds-guide.gtf" --allow-cds-bridges --report "$tmpdir/cds-allow.tsv" \
  -o "$tmpdir/cds-allow.gtf" "$tmpdir/cds-input.gtf"
has_tx "$tmpdir/cds-allow.gtf" chrC 100 599 +
grep -q 'cds_bridge_allowed' "$tmpdir/cds-allow.tsv"


## ---------------------------------------------------------------------------
## Fixtures for the observed-chains-first merge policy.
## known_fail marks a documented defect that is not fixed yet.
t=$tmpdir
FULL=1000-1100,2000-2100,3000-3100,4000-4300
FULLX=1000-1100,2000-2100,3000-3100,4000-5300
## every fixture also checks the reconciliation: each input record has a fate
fx() {
  fx_out=$1; shift
  run_stringtie --merge -o "$t/$fx_out.gtf" --report "$t/$fx_out.tsv" "$@"
  if ! grep -q "^# fate	unaccounted	0$" "$t/$fx_out.tsv"; then
    printf '%s\n' "input records without a fate in $fx_out:" >&2
    grep '^# \(input_records\|fate\)' "$t/$fx_out.tsv" >&2
    exit 1
  fi
}
## fate_is FILE NAME COUNT
fate_is() {
  grep -q "^# fate	$2	$3$" "$1" || { printf '%s\n' "fate $2 is not $3 in $1" >&2; grep '^# fate' "$1" >&2; exit 1; }
}

## G/H/H2/I: a single-file terminal-exon extension must not replace a 2-file end
mkgtf $t/g1.gtf in "fxG + g1 $FULL"
mkgtf $t/g2.gtf in "fxG + g2 $FULL"
mkgtf $t/g3.gtf in "fxG + g3 2050-2100,3000-3100,4000-5300"
fx G $t/g1.gtf $t/g2.gtf $t/g3.gtf
need_chain $t/G.gtf fxG + $FULL
test "$(tx_count $t/G.gtf fxG)" -eq 1
grep -Eq '^fxG.*absorb[[:space:]]+contained_fragment' $t/G.tsv
fate_is $t/G.tsv emitted_chain 2
fate_is $t/G.tsv contained_fragment 1
## the same partial chain in 3 files is a replicated longer end: kept as its own transcript
mkgtf $t/g4.gtf in "fxG + g4 2050-2100,3000-3100,4000-5300"
mkgtf $t/g5.gtf in "fxG + g5 2050-2100,3000-3100,4000-5300"
fx G3 $t/g1.gtf $t/g2.gtf $t/g3.gtf $t/g4.gtf $t/g5.gtf
need_chain $t/G3.gtf fxG + $FULL
need_chain $t/G3.gtf fxG + 2050-2100,3000-3100,4000-5300 merge_support 3

mkgtf $t/h3.gtf in "fxG + h3 $FULLX"
fx H $t/g1.gtf $t/g2.gtf $t/h3.gtf
need_chain $t/H.gtf fxG + $FULL
test "$(tx_count $t/H.gtf fxG)" -eq 1
grep -Eq '^fxG.*trim[[:space:]]+unsupported_end_extension' $t/H.tsv
fx H2 --max-extension-3p 100 $t/g1.gtf $t/g2.gtf $t/h3.gtf
need_chain $t/H2.gtf fxG + $FULL
## --min-support-3p 1 restores the most extended end; the hard cap still wins
fx H3 --min-support-3p 1 $t/g1.gtf $t/g2.gtf $t/h3.gtf
need_chain $t/H3.gtf fxG + $FULLX
## slop: a 200 bp single-file overhang is adopted only if both slop limits allow it,
## and never beyond the hard cap
FULL45=1000-1100,2000-2100,3000-3100,4000-4500
mkgtf $t/h5.gtf in "fxG + h5 $FULL45"
fx H5 $t/g1.gtf $t/g2.gtf $t/h5.gtf
need_chain $t/H5.gtf fxG + $FULL
fx H6 --slop-frac-3p 1 $t/g1.gtf $t/g2.gtf $t/h5.gtf
need_chain $t/H6.gtf fxG + $FULL45
fx H7 --slop-frac-3p 1 --slop-bp-3p 150 $t/g1.gtf $t/g2.gtf $t/h5.gtf
need_chain $t/H7.gtf fxG + $FULL
fx H8 --slop-frac-3p 1 --max-extension-3p 100 $t/g1.gtf $t/g2.gtf $t/h5.gtf
need_chain $t/H8.gtf fxG + $FULL

## minus-strand mirror of H: the extended end is the 5-prime end
mkgtf $t/hm1.gtf in "fxHm - hm1 $FULL"
mkgtf $t/hm2.gtf in "fxHm - hm2 $FULL"
mkgtf $t/hm3.gtf in "fxHm - hm3 $FULLX"
fx Hm $t/hm1.gtf $t/hm2.gtf $t/hm3.gtf
need_chain $t/Hm.gtf fxHm - $FULL
fx Hm5 --min-support-5p 1 $t/hm1.gtf $t/hm2.gtf $t/hm3.gtf
need_chain $t/Hm5.gtf fxHm - $FULLX

mkgtf $t/i3.gtf in "fxG + i3 300-1100,2000-2100,3000-3100,4000-4300"
fx I $t/g1.gtf $t/g2.gtf $t/i3.gtf
need_chain $t/I.gtf fxG + $FULL
test "$(tx_count $t/I.gtf fxG)" -eq 1

## S: --min-support counts an intron chain, not base-identical ends
mkgtf $t/s1.gtf in "fxS + s1 990-1100,2000-2100,3000-3290"
mkgtf $t/s2.gtf in "fxS + s2 1000-1100,2000-2100,3000-3300"
mkgtf $t/s3.gtf in "fxS + s3 1010-1100,2000-2100,3000-3310"
fx S2 --min-support 2 $t/s1.gtf $t/s2.gtf $t/s3.gtf
test "$(tx_count $t/S2.gtf fxS)" -eq 1
fx S3 --min-support 3 $t/s1.gtf $t/s2.gtf $t/s3.gtf
test "$(tx_count $t/S3.gtf fxS)" -eq 1
fx S4 --min-support 4 $t/s1.gtf $t/s2.gtf $t/s3.gtf
test "$(tx_count $t/S4.gtf fxS)" -eq 0

## Sd: without --min-support the threshold is 2 for 4 or more input files, 1 otherwise
ONE=1000-1100,2000-2100,3000-3300
TWO=11000-11100,12000-12100,13000-13300
mkgtf $t/sd1.gtf in "fxSd + sd1 $ONE" "fxSd + sd5 $TWO"
mkgtf $t/sd2.gtf in "fxSd + sd2 $TWO"
mkgtf $t/sd3.gtf in "fxSd + sd3 21000-21100,22000-22300"
mkgtf $t/sd4.gtf in "fxSd + sd4 31000-31100,32000-32300"
fx Sd3 $t/sd1.gtf $t/sd2.gtf $t/sd3.gtf
need_chain $t/Sd3.gtf fxSd + $ONE
need_chain $t/Sd3.gtf fxSd + $TWO
grep -q "^# min_support	1$" $t/Sd3.tsv
fx Sd4 $t/sd1.gtf $t/sd2.gtf $t/sd3.gtf $t/sd4.gtf
no_chain $t/Sd4.gtf fxSd + $ONE
need_chain $t/Sd4.gtf fxSd + $TWO merge_support 2
test "$(tx_count $t/Sd4.gtf fxSd)" -eq 1
grep -q "^# min_support	2$" $t/Sd4.tsv
grep -Eq '^fxSd.*minimum_support' $t/Sd4.tsv
grep -q "^# input_records	5$" $t/Sd4.tsv
fate_is $t/Sd4.tsv emitted_chain 2
fate_is $t/Sd4.tsv minimum_support 3
fx Sd41 --min-support 1 $t/sd1.gtf $t/sd2.gtf $t/sd3.gtf $t/sd4.gtf
need_chain $t/Sd41.gtf fxSd + $ONE
test "$(tx_count $t/Sd41.gtf fxSd)" -eq 4

## SET support is pairwise reciprocal overlap, not exact ends and not transitive
mkgtf $t/ss1.gtf in "fxSS + ss1 1000-2000" "fxSS + ss3 5000-5300"
mkgtf $t/ss2.gtf in "fxSS + ss2 1050-2050" "fxSS + ss4 5200-5500"
fx SS2 --min-support 2 $t/ss1.gtf $t/ss2.gtf
fate_is $t/SS2.tsv emitted_set 2
fate_is $t/SS2.tsv minimum_support 2
test "$(tx_count $t/SS2.gtf fxSS)" -eq 1
grep -Eq '^fxSS.*(5000[[:space:]]+5300|5200[[:space:]]+5500).*minimum_support' $t/SS2.tsv

## Q: a replicated isoform is not erased by a single-file longer model
mkgtf $t/q1.gtf in "fxQ + q1 1000-1100,2000-2100,3000-3300" "fxQ + q1b 1950-2100,3000-3300"
mkgtf $t/q2.gtf in "fxQ + q2 1950-2100,3000-3300"
mkgtf $t/q3.gtf in "fxQ + q3 1950-2100,3000-3300"
fx Q $t/q1.gtf $t/q2.gtf $t/q3.gtf
need_chain $t/Q.gtf fxQ + 1000-1100,2000-2100,3000-3300
need_chain $t/Q.gtf fxQ + 1950-2100,3000-3300 merge_support 3
need_chain $t/Q.gtf fxQ + 1000-1100,2000-2100,3000-3300 merge_support 1
## Q2: with equal support the fragment is absorbed, unless intron overhang is disallowed
mkgtf $t/q4.gtf in "fxQ + q4 1000-1100,2000-2100,3000-3300"
mkgtf $t/q5.gtf in "fxQ + q5 1000-1100,2000-2100,3000-3300"
fx Q2 $t/q1.gtf $t/q2.gtf $t/q3.gtf $t/q4.gtf $t/q5.gtf
need_chain $t/Q2.gtf fxQ + 1000-1100,2000-2100,3000-3300 merge_support 5
need_chain $t/Q2.gtf fxQ + 1000-1100,2000-2100,3000-3300 merge_path_support 3
test "$(tx_count $t/Q2.gtf fxQ)" -eq 1
grep -Eq '^fxQ.*absorb[[:space:]]+contained_fragment.*records=3 files=3' $t/Q2.tsv
fx Q2s --max-intron-overhang 0 $t/q1.gtf $t/q2.gtf $t/q3.gtf $t/q4.gtf $t/q5.gtf
need_chain $t/Q2s.gtf fxQ + 1950-2100,3000-3300
test "$(tx_count $t/Q2s.gtf fxQ)" -eq 2

## L: a SET inside an exon of an opposite-strand MET is kept
mkgtf $t/l1.gtf in "fxL + l1 1000-1100,2000-2600,3000-3100" "fxL - l2 2100-2500"
fx L $t/l1.gtf
need_chain $t/L.gtf fxL + 1000-1100,2000-2600,3000-3100
need_chain $t/L.gtf fxL - 2100-2500
fate_is $t/L.tsv emitted_chain 1
fate_is $t/L.tsv emitted_set 1

## N/U: unstranded inputs pass through unchanged
mkgtf $t/n1.gtf in "fxN . n1 8000-8100,9000-9100"
fx N $t/n1.gtf
need_chain $t/N.gtf fxN . 8000-8100,9000-9100
mkgtf $t/u1.gtf in "fxU . u1 1000-1300" "fxU . u2 1200-1500"
fx U $t/u1.gtf
need_chain $t/U.gtf fxU . 1000-1300
need_chain $t/U.gtf fxU . 1200-1500
test "$(tx_count $t/U.gtf fxU)" -eq 2
## unstranded chain in 2 files: one model, support 2; rejected at --min-support 3
mkgtf $t/n2.gtf in "fxN . n2 8010-8100,9000-9090"
fx N2 --min-support 2 $t/n1.gtf $t/n2.gtf
need_chain $t/N2.gtf fxN . 8000-8100,9000-9100 merge_support 2
test "$(tx_count $t/N2.gtf fxN)" -eq 1
fx N3 --min-support 3 $t/n1.gtf $t/n2.gtf
test "$(tx_count $t/N3.gtf fxN)" -eq 0

## A/E/V: guide relation for partial chains and shifted variants
mkgtf $t/ag.gtf ref "fxA + guideA $FULL"
mkgtf $t/a1.gtf in "fxA + a1 2050-2100,3000-3100,4000-5300"
mkgtf $t/a2.gtf in "fxA + a2 2050-2100,3000-3100,4000-5300"
mkgtf $t/a3.gtf in "fxA + a3 2050-2100,3000-3100,4000-5300"
AV=2050-2100,3000-3100,4000-5300
## A: a single-file partial chain past a guide end is absorbed, credited and reported;
## with endpoint support it is a guide variant, never silently swallowed
fx A -G $t/ag.gtf $t/a1.gtf
need_chain $t/A.gtf fxA + $FULL merge_contained_support 1
test "$(tx_count $t/A.gtf fxA)" -eq 1
grep -Eq '^fxA.*absorb[[:space:]]+contained_fragment.*guideA' $t/A.tsv
fx D --min-support-3p 1 -G $t/ag.gtf $t/a1.gtf
need_chain $t/D.gtf fxA + $FULL
need_chain $t/D.gtf fxA + $AV reference_id guideA
test "$(tx_count $t/D.gtf fxA)" -eq 2
fx Ar --min-support-3p 1 --guide-policy reference -G $t/ag.gtf $t/a1.gtf
test "$(tx_count $t/Ar.gtf fxA)" -eq 1
grep -Eq '^fxA.*reject[[:space:]]+guide_policy_reference' $t/Ar.tsv
fx Ac --min-support-3p 1 --max-extension-3p 500 -G $t/ag.gtf $t/a1.gtf
test "$(tx_count $t/Ac.gtf fxA)" -eq 1
grep -Eq '^fxA.*reject[[:space:]]+guide_extension_cap' $t/Ac.tsv
fx E -G $t/ag.gtf $t/a1.gtf $t/a2.gtf $t/a3.gtf
need_chain $t/E.gtf fxA + $FULL
need_chain $t/E.gtf fxA + $AV reference_id guideA
need_chain $t/E.gtf fxA + $AV merge_support 3
test "$(tx_count $t/E.gtf fxA)" -eq 2
## a fragment inside the guide is absorbed, credited and reported
mkgtf $t/af.gtf in "fxA + af 2050-2100,3000-3100,4000-4200"
fx Af -G $t/ag.gtf $t/af.gtf
test "$(tx_count $t/Af.gtf fxA)" -eq 1
need_chain $t/Af.gtf fxA + $FULL merge_contained_support 1
grep -Eq '^fxA.*absorb[[:space:]]+contained_fragment.*guideA' $t/Af.tsv

## V: same chain as the guide, longer at one end and shorter at the other
mkgtf $t/vg.gtf ref "fxV + guideV 1000-1199,1300-1499"
mkgtf $t/v1.gtf in "fxV + v1 900-1199,1300-1399"
VV=900-1199,1300-1399
fx V -G $t/vg.gtf $t/v1.gtf
need_chain $t/V.gtf fxV + 1000-1199,1300-1499 merge_support 1
test "$(tx_count $t/V.gtf fxV)" -eq 1
fx Vs --min-support-5p 1 --guide-slop-bp-5p 50 -G $t/vg.gtf $t/v1.gtf
need_chain $t/Vs.gtf fxV + $VV reference_id guideV
fx Vr --min-support-5p 1 --guide-slop-bp-5p 50 --guide-policy reference -G $t/vg.gtf $t/v1.gtf
no_chain $t/Vr.gtf fxV + $VV
fx Vc --min-support-5p 1 --guide-slop-bp-5p 50 --max-extension-5p 0 -G $t/vg.gtf $t/v1.gtf
no_chain $t/Vc.gtf fxV + $VV
need_chain $t/Vc.gtf fxV + 1000-1199,1300-1499

## W: --path-mode direct must not destroy the observed inputs of a rejected join
mkgtf $t/w1.gtf in "fxW + w1 1000-1100,2000-2100,3000-3300"
mkgtf $t/w2.gtf in "fxW + w2 2000-2100,3000-3300,4000-4300"
fx Wd --path-mode direct --min-support-5p 1 --min-support-3p 1 $t/w1.gtf $t/w2.gtf
no_chain $t/Wd.gtf fxW + 1000-1100,2000-2100,3000-3300,4000-4300
need_chain $t/Wd.gtf fxW + 1000-1100,2000-2100,3000-3300
need_chain $t/Wd.gtf fxW + 2000-2100,3000-3300,4000-4300
fx Ws --min-support-5p 1 --min-support-3p 1 $t/w1.gtf $t/w2.gtf
need_chain $t/Ws.gtf fxW + 1000-1100,2000-2100,3000-3300,4000-4300
test "$(tx_count $t/Ws.gtf fxW)" -eq 1

## X: malformed numeric option values are rejected
bad_opt() {
  if run_stringtie --merge "$@" -o $t/bad.gtf $t/w1.gtf >/dev/null 2>&1; then return 1; fi
}
bad_opt --min-support 2x
bad_opt --set-overlap nan
bad_opt --slop-frac-3p x
bad_opt --min-introns abc

## thread determinism: same structures and same report at -p 1 and -p 4
fx P1 -p 1 -G "$tmpdir/guide.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf" $t/g1.gtf $t/q1.gtf $t/q2.gtf
fx P4 -p 4 -G "$tmpdir/guide.gtf" "$tmpdir/in1.gtf" "$tmpdir/in2.gtf" $t/g1.gtf $t/q1.gtf $t/q2.gtf
for r in chr1 fxG fxQ; do
  test "$(tx_chains $t/P1.gtf $r + | cut -f1 | sort)" = "$(tx_chains $t/P4.gtf $r + | cut -f1 | sort)"
done
cmp -s $t/P1.tsv $t/P4.tsv
## the output itself, IDs included, does not depend on the number of threads
grep -v '^#' $t/P1.gtf > $t/P1.body
grep -v '^#' $t/P4.gtf > $t/P4.body
cmp -s $t/P1.body $t/P4.body

## regression fixtures from the second code review
## unstranded pass-through with a guide in the bundle (the read list is re-sorted)
mkgtf $t/r1a.gtf in "fxR1 . r1u 1000-1500" "fxR1 + r1p 2000-2200,2400-2600"
mkgtf $t/r1b.gtf in "fxR1 + r1q 2000-2200,2400-2600"
mkgtf $t/r1g.gtf ref "fxR1 + guideR1 100-300,2400-2600"
fx R1 -G $t/r1g.gtf $t/r1a.gtf $t/r1b.gtf
need_chain $t/R1.gtf fxR1 . 1000-1500
need_chain $t/R1.gtf fxR1 + 100-300,2400-2600
test "$(tx_chains $t/R1.gtf fxR1 . | wc -l)" -eq 1
## --min-introns 0 does not bypass end election
mkgtf $t/r5a.gtf in "fxR5 + r5a 500-1200,1500-2000"
mkgtf $t/r5b.gtf in "fxR5 + r5b 1100-1200,1500-2000"
fx R5 --min-introns 0 $t/r5a.gtf $t/r5b.gtf
need_chain $t/R5.gtf fxR5 + 1100-1200,1500-2000
## a guide whose full chain is unobserved still has piece support from its fragments
mkgtf $t/r8g.gtf ref "fxR8 + guideR8 1000-1200,1500-1700,2000-2300"
mkgtf $t/r8a.gtf in "fxR8 + r8a 1000-1200,1500-1700" "fxR8 + r8b 1500-1700,2000-2300"
mkgtf $t/r8b.gtf in "fxR8 + r8c 1000-1200,1500-1700" "fxR8 + r8d 1500-1700,2000-2300"
fx R8 --guide-fragments dominance -G $t/r8g.gtf $t/r8a.gtf $t/r8b.gtf
need_chain $t/R8.gtf fxR8 + 1000-1200,1500-1700,2000-2300 merge_piece_support 2
need_chain $t/R8.gtf fxR8 + 1000-1200,1500-1700,2000-2300 merge_contained_support 2
test "$(tx_count $t/R8.gtf fxR8)" -eq 1
## SET support counts identical records that were collapsed
mkgtf $t/r2a.gtf in "fxR2 + r2a 1000-2000"
mkgtf $t/r2b.gtf in "fxR2 + r2b 1000-2000"
mkgtf $t/r2c.gtf in "fxR2 + r2c 1050-2100"
mkgtf $t/r2d.gtf in "fxR2 + r2d 1050-2100"
fx R2 $t/r2a.gtf $t/r2b.gtf $t/r2c.gtf $t/r2d.gtf
test "$(tx_count $t/R2.gtf fxR2)" -eq 1
grep -q 'merge_support "4"' $t/R2.gtf
## a SET within the guide slop of a single-exon guide is the guide
mkgtf $t/r3g.gtf ref "fxR3 + guideR3 1000-2000"
mkgtf $t/r3a.gtf in "fxR3 + r3a 900-2100"
mkgtf $t/r3b.gtf in "fxR3 + r3b 900-2100"
fx R3 -G $t/r3g.gtf $t/r3a.gtf $t/r3b.gtf
test "$(tx_count $t/R3.gtf fxR3)" -eq 1
need_chain $t/R3.gtf fxR3 + 1000-2000 merge_support 2
fate_is $t/R3.tsv emitted_set 2
## records dropped before the bundle are counted too: -m, and a SET inside a MET exon
mkgtf $t/f1.gtf in "fxF + f1 1000-1100,2000-2300" "fxF + f2 2050-2250" "fxF + f3 5000-5030"
fx F $t/f1.gtf
grep -q "^# input_records	3$" $t/F.tsv
fate_is $t/F.tsv emitted_chain 1
fate_is $t/F.tsv contained_set 1
fate_is $t/F.tsv min_length 1

## O: output order. References keep the order of the input files (not alphabetical), then
## guide-only references; transcripts are sorted by start, a guide without input support is
## written in place and not appended at the end.
mkgtf $t/o1.gtf in "chr2 + o1 5000-5100,6000-6300" "chr10 + o2 1000-1100,2000-2300" "chr10 - o5 1500-1900"
mkgtf $t/o2.gtf in "chr10 + o3 8000-8100,9000-9300" "chr1 + o4 1000-1100,2000-2300"
mkgtf $t/og.gtf ref "chr10 + guideO1 4000-4100,5000-5300" "chrZ + guideO2 100-200,300-400" "chr2 + guideO3 100-200,300-400"
fx O -G $t/og.gtf $t/o1.gtf $t/o2.gtf
test "$(awk -F'\t' '$3=="transcript" {printf "%s:%s ",$1,$4}' $t/O.gtf)" = \
  "chr2:100 chr2:5000 chr10:1000 chr10:1500 chr10:4000 chr10:8000 chr1:1000 chrZ:100 "
fx O4 -p 4 -G $t/og.gtf $t/o1.gtf $t/o2.gtf
test "$(grep -v '^#' $t/O.gtf)" = "$(grep -v '^#' $t/O4.gtf)"

## C: a guide named like a generated ID does not collide with one. Guide transcript MSTRG.1.1
## and guide gene MSTRG.2 reserve gene numbers 1 and 2.
mkgtf $t/cg.gtf ref "fxC + MSTRG.1.1 1000-1100,2000-2300"
awk -F'\t' -v OFS='\t' '{sub(/gene_id "[^"]*"/,"gene_id \"MSTRG.2\"",$9); print}' $t/cg.gtf > $t/cg2.gtf
mkgtf $t/c1.gtf in "fxC + c1 11000-11100,12000-12300" "fxC + c2 21000-21100,22000-22300"
fx C -G $t/cg2.gtf $t/c1.gtf
ids=$(awk -F'\t' '$3=="transcript" {id=$9; sub(/.*transcript_id "/,"",id); sub(/".*/,"",id); print id}' $t/C.gtf)
test "$(printf '%s\n' "$ids" | sort | tr '\n' ' ')" = "MSTRG.1.1 MSTRG.3.1 MSTRG.4.1 "
test "$(awk -F'\t' '$3=="transcript" && $9 ~ /gene_id "MSTRG.2"/' $t/C.gtf | wc -l)" -eq 1

if test -n "${known_failures:-}"; then
  printf '%s\n' "merge policy test: known failures:$known_failures"
fi
printf '%s\n' 'merge policy test: OK'
