# Safer and configurable `--merge` meta-assembly

## Approved defaults

- use `-F 0`, `-T 0`, and `-f 0` in merge mode only;
- count structural support by distinct input file;
- retain singleton structures with `--min-support 1`;
- require one compatible exact shared intron for MET coassembly;
- evaluate `-G` gene and CDS barriers before shared-intron requirements;
- reject same-strand cross-gene and same-gene CDS bridges by default;
- keep guides unchanged and emit allowed extensions under new IDs by default;
- swallow contained SETs and require replicated evidence for non-minor SET
  extensions;
- provide independent 5-prime and 3-prime support, slop, and hard-cap knobs;
- optionally write a consolidated merge event and rejection report.

## Interface

Add merge-only options without a redundant `--merge-` prefix:

```text
--min-support N                 --rescue-support N
--min-introns N                 --path-mode shared|replicated|direct
--min-piece-support N           --set-overlap F
--min-support-5p N              --min-support-3p N
--slop-bp-5p N                  --slop-bp-3p N
--slop-frac-5p F                --slop-frac-3p F
--max-extension-5p N            --max-extension-3p N
--guide-policy both|reference|best
--guide-replace-frac F          --dedup-guides
--allow-fusions                 --allow-cds-bridges
--report FILE
```

Emit transparent `merge_*` support attributes, reference relationship tags,
and fusion metadata in GTF output.

## Implementation

1. Preserve source file, gene, CDS, and quantitative metadata through GTF input
   conversion, and keep structural evidence separate from normalized abundance.
2. Aggregate input support by canonical SET or MET structure and distinct file.
   Quantitative failures can be rescued by replicated structural support.
3. Enforce guide gene/CDS barriers first. Within allowed regions, apply SET/MET,
   shared-intron, path-mode, endpoint-support, and hard-cap policies.
4. Preserve reference IDs only for unchanged guides. Deduplicate generated
   structures and apply the selected guide representative policy.
5. Produce a deterministic optional TSV report with summary comments and one
   consolidated row per accepted or rejected event.
6. Update CLI documentation, `MERGE_MODE.md`, and the `AGENTS.md` glossary.

## Verification

- add focused fixtures for option validation, support counting, path modes,
  intron compatibility, SET behavior, per-end behavior, guide policies, gene
  and CDS barriers, fusion overrides, missing abundance, and report output;
- verify deterministic results across thread counts;
- run the issue #247 dataset and require unchanged guide structures, unique
  transcript IDs, no merge-created duplicate structures, and reconciled report
  totals;
- verify ordinary assembly retains `-f 0.01`.

Do not modify or commit `worklogs/planning.txt`.
