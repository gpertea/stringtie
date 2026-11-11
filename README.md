# StringTie development notes

This branch holds working notes from the development of StringTie: measurements,
analyses and plans that explain why the code is the way it is. It has no history
in common with `master` and is never merged into it.

```text
File                                       About
-----------------------------------------  -------------------------------------------
merge-issues.md                            problems found in --merge before the 3.0.4 rework
merge-default-vs-3.0.3.md                  what the new --merge defaults change, measured
assembly-observations-polyA-ribozero.md    what single-sample assemblies look like
                                           (polyA vs RiboZero, effect of -N)
plan_26-08-09_02-30_merge-dev.md           plan for the merge rework
work_26-08-09_09-57_merge-dev.md           work log of the merge rework
```

To have the notes next to a checkout of the code:

```sh
git worktree add ../stringtie-worklogs worklogs
```
