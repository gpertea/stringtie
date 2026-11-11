# What single-sample StringTie assemblies look like: polyA vs RiboZero, with and without -N

Written 2026-11-11. These are measurements made while evaluating `--merge`,
kept because they describe the per-sample assemblies themselves and may help
when tuning the read-alignment assembler. Numbers are facts from the runs
named below; the last section is interpretation and is labeled as such.

## 1. Data

- 10 human brain RNA samples, each sequenced twice: a polyA library
  (unstranded) and a RiboZero total-RNA library (fr-firststrand).
- Assembled per sample with StringTie 3.0.3, default parameters, no `-G`:

```text
Set            Flags      What it is
-------------  ---------  -------------------------------------------
polyA          (none)     unstranded polyA
riboZ_default  --rf       stranded total RNA
riboZ_N        --rf -N    stranded total RNA, nascent-aware assembly
```

- Truth proxy: GENCODE 49, 475,887 distinct intron chains. "Matches GENCODE"
  means the same strand and intron chain as an annotated transcript; ends are
  ignored. GENCODE is incomplete and most of it is not expressed in one
  tissue, so the precision figures are lower bounds and the sensitivity
  figures are small by construction. Compare rows, not absolute values.
- "Chain" is a distinct stranded intron chain. "Record" is one assembled
  transcript in one sample.
- Files: assemblies in `merge-dev/brainseq_eval/{polyA,riboZ_default,riboZ_N}/`,
  results and scripts in `merge-dev/brainseq_eval/survey/` (`per_sample/`,
  `*.supk.txt`, `*.endsurvey.txt`, `*.overhang.txt`, `scripts/`).

## 2. Per-sample assemblies (mean of 10 samples; `per_sample/summary.txt`)

```text
Set            Transcripts  Single-exon  Multi-exon  Mean   Chains in  Chain
                                         chains      exons  GENCODE    precision
-------------  -----------  -----------  ----------  -----  ---------  ---------
polyA               57,206   2,434 (4%)      54,771    8.7     24,128      44.2%
riboZ_default       63,610   3,346 (5%)      60,264    7.8     20,504      34.2%
riboZ_N             56,447   2,748 (5%)      53,699    8.9     20,401      38.1%
```

- polyA finds about 18% more annotated chains per sample than RiboZero from
  the same RNA, at higher precision. Range over samples: polyA 40 to 48%,
  riboZ_default 31 to 42%, riboZ_N 35 to 44%.
- In the unstranded polyA assemblies every multi-exon transcript has a strand
  (from its splice sites) and every single-exon transcript is unstranded
  (1,600 to 4,400 per sample). In the stranded RiboZero assemblies almost none
  are unstranded (0 to 15 per sample).

## 3. What -N changes (same RiboZero sample assembled both ways)

```text
Multi-exon chains per sample        Count    In GENCODE
----------------------------------  -------  ----------
in both assemblies                   42,940       43.6%
only without -N                      17,324       10.6%
only with -N                         10,760       16.0%
```

- `-N` is not a filter on the default assembly: it drops about 17,000 chains
  and produces about 10,800 different ones per sample.
- The number of annotated chains is unchanged (20,504 without, 20,401 with).
- The chains it drops are the least reliable group (10.6% annotated). The ones
  it adds are better (16.0%) but still far below the shared ones.
- Transcripts get longer: mean exons per multi-exon transcript 7.8 -> 8.9.
- Over 10 samples the gain persists but shrinks: distinct chains seen once,
  211,099 with `-N` against 261,145 without; precision of all chains 16.0%
  against 13.7%.

## 4. Replication across samples is the strongest predictor of correctness

Distinct chains over the 10 polyA samples, by the number of samples that
assembled them (`polyA.supk.txt`; RiboZero is within 2 to 3 points):

```text
Samples  Chains in GENCODE  Novel chains  Share in GENCODE
-------  -----------------  ------------  ----------------
1                   17,863       174,736              9.3%
2                    7,247        23,362             23.7%
3                    4,434         8,138             35.3%
5                    2,494         2,146             53.8%
8                    2,127           777             73.2%
10                  10,148           677             93.7%
```

- 71% of all distinct chains (192,599 of 271,077) come from a single sample,
  and 91% of those are not annotated.
- A chain assembled in all 10 samples is annotated 94% of the time.
- 10 samples together contain 54,316 GENCODE chains; one sample has about
  24,000. Of 94,819 GENCODE chains hidden in a hold-out, 10,704 were present
  in at least one polyA assembly, 3,600 of them in only one sample.

## 5. Transcript ends

Records whose intron chain is a complete GENCODE chain, end compared with the
annotated end of that transcript (`*.endsurvey.txt`, rows `fullchain_*`):

```text
                     polyA 5'   polyA 3'   riboZ_N 5'   riboZ_N 3'
-------------------  ---------  ---------  -----------  -----------
identical                 0.9%       5.5%         0.7%         2.7%
longer than annot.       47.5%      32.4%        36.0%        38.8%
shorter than annot.      51.6%      62.1%        63.3%        58.5%
within 100 bp              70%        52%          66%          49%
within 500 bp              92%        68%          92%          64%
> 1 kb longer             2.5%      13.0%         2.0%        15.1%
> 1 kb shorter            0.9%      10.4%         1.3%        11.8%
```

- 5' ends are within 500 bp of the annotation for 92% of records in both
  libraries. 3' ends are the weak point: a third are off by more than 500 bp.
- 3' ends are more often too short than too long, but when they are too long
  they are far too long: 13 to 15% of records run more than 1 kb past the
  annotated end.
- The same chain assembled in two samples almost never has identical ends:
  among chains seen in 2 or more samples, one exact structure recurs for 5.2%
  of annotated and 1.2% of novel chains.
- Between samples, for novel chains seen in 3 or more samples, the most
  extended end is within 500 bp of the second most extended one for:

```text
Set            5'     3'
-------------  -----  -----
polyA          90%    85%
riboZ_N        88%    78%
riboZ_default  86%    75%
```

  So 3' ends vary more than 5' ends, more in total RNA than in polyA, and `-N`
  recovers part of the difference.

## 6. Fragments and intron overhang

- A quarter of multi-exon records are a contiguous part of an annotated chain
  rather than the whole chain: 24.9% (polyA), 26.6% (riboZ_N), 28.8%
  (riboZ_default).
- Where such a fragment ends next to an internal exon of the annotated
  transcript, the end stays inside that exon in under half of the cases
  (`*.overhang.txt`):

```text
Set            End  Inside the exon  Overhangs into the intron  Spans the intron
-------------  ---  ---------------  -------------------------  ----------------
polyA          5'             46.4%                      39.6%             14.1%
polyA          3'             43.5%                      37.0%             19.6%
riboZ_N        5'             42.8%                      46.9%             10.3%
riboZ_N        3'             41.7%                      43.1%             15.2%
riboZ_default  5'             38.2%                      48.3%             13.6%
riboZ_default  3'             33.1%                      49.4%             17.5%
```

- The overhang exceeds 1 kb in 5.5% (5') and 10.6% (3') of polyA fragment ends,
  and in 10.7 to 17.5% for RiboZero.

## 7. Transcripts that touch more than one annotated gene

In the polyA assemblies, 13.8% of records have exons overlapping exons of two
or more GENCODE genes on the same strand (tested against a GENCODE subset
holding 80% of its transcripts). Most of that is the annotation, not
the assembler: 86% involve genes whose spans overlap each other (nested genes,
small RNA genes inside a host, readthrough genes). What remains:

- 10,745 records (1.9% of all) connect genes at separate loci.
- 8,955 of those (1.6% of all) have an intron chain that is not annotated.

RiboZero -N gives the same proportions.

## 8. Reading these numbers (interpretation, not measurement)

- The per-sample precision problem is mostly sample-specific chains. Seven in
  ten distinct chains appear in one sample only and nine in ten of those are
  unannotated, while replicated chains are far better. Any assembler change
  that removes chains the other samples do not reproduce should move
  precision more than changes aimed at well-supported loci. The replication
  count is available as a free, annotation-independent quality label for
  evaluating such changes on this dataset.
- `-N` helps total RNA mainly by not producing a class of unreliable chains
  and by yielding longer transcripts, without costing annotated chains. It
  does not make RiboZero as clean as polyA: about 6 points of precision and
  3,700 annotated chains per sample remain between them.
- 3' end placement is the least reliable part of an assembled transcript, in
  both directions, and worse in total RNA. Long 3' overruns (more than 1 kb in
  13 to 15% of records) are the kind that fuse neighboring genes downstream.
- Fragment ends more often sit in an intron than in the exon, also in polyA.
  That fits incompletely processed RNA and coverage that decays into introns,
  and suggests the rule that ends a transcript inside intronic coverage as a
  place to look.
- Single-exon transcripts from an unstranded library carry no strand, so every
  downstream step has to treat them as compatible with both strands or leave
  them alone.

These are hypotheses from one tissue and one annotation. Each is checkable on
the same files with the scripts in `merge-dev/brainseq_eval/survey/scripts/`.
