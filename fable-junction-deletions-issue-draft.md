Long reads: deletions next to introns give wrong junction evidence (4 defects in the juncsdel / "DEL AWARE" code)

Found while working on #511. All four are in the code that handles a deletion placed
by the aligner right before or right after an intron (`-L` and `--mix`). They are
independent of the #511 slowdown, but the first three inflate the number of
"distinct" alignments, and all four can change the assembly. Checked on master
53d9647 (3.0.4); the code is the same in 3.0.3.

## How this code is meant to work (as I read it)

`GSamRecord::setupCoordinates()` records in `juncsdel`, for each intron, the deleted
bases before it (`start`) and after it (`end`). For a junction with such a deletion,
`processRead()` stores a tentative unstranded junction widened by the deleted bases
(`rlink.cpp:850`). In `count_good_junctions()` (`rlink.cpp:17609`), if a stranded
junction with the widened coordinates exists (other reads were aligned without the
deletion), the read is moved to it and its exon ends are adjusted; otherwise the
read's own junction is restored. `deljuncmatch()` keeps alignments with different
deletion offsets from being collapsed into one stored read.

## Defects

1. **A deletion before an intron is recorded as a deletion after it.**
   `gclib/GSam.cpp:228` builds the entry with `GSeg deljunc(prevdel,0)`. The `GSeg`
   constructor swaps its arguments when the first is larger (`GBase.h:381`), so the
   entry becomes `(0,prevdel)`. For CIGAR `..100M 2D 300N 100M..` the junction end is
   moved 2 bases into the next exon instead of the junction start being moved back
   by 2. The widened junction never matches the clean one, so the rescue does not
   happen for this case. `GSeg deljunc(prevdel,del)` at line 239 is swapped the same
   way when `prevdel>del`.

2. **A deletion after an intron adds two `juncsdel` entries for one intron.**
   For `..M 498N 2D M..` the `D` case adds `(prevdel,del)` and leaves `intron` set,
   then the `M` case adds `(prevdel,0)` (`gclib/GSam.cpp:227-241`). Consequences:
   `processRead()` reads the offsets of every later junction of that read from the
   wrong entry (`brec.juncsdel[i-1]`), and `deljuncmatch()` fails on the count
   comparison, so identical alignments of this kind are never collapsed.

3. **With `-G`, identical alignments with a deletion at a guide junction are never
   collapsed.** A junction matching a guide intron is stored without the deletion
   offsets (`rlink.cpp:851`), but `deljuncmatch()` compares the offsets derived from
   the stored junction (0) with those of the incoming alignment (non-zero).

4. **The rescue adds the shared junction's total count once per read.**
   At `rlink.cpp:17617` and `17629`, `nj->nreads+=rd.juncs[i-1]->nreads` (and `mm`,
   `nm`) adds the full counts of the tentative junction for every read that has it.
   With k such reads the target junction gains k times their total. This changes
   which of two nearby splice sites is treated as the error.
   (`rlink.cpp:17788` has the same pattern for unstranded junctions that get a
   strand; I did not change or test that one.)

## Reproducers

All reads start at chrT:1000 on the + strand. `mkreads.py`:

```python
import sys,random,re
random.seed(7)
HDR="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:100000\n"
def seq(n): return "".join(random.choice("ACGT") for _ in range(n))
def qlen(cig): return sum(int(n) for n,op in re.findall(r'(\d+)([MIDNS=X])',cig) if op in "MIS=X")
out=sys.argv[1]; groups=[a.split(":") for a in sys.argv[2:]]
with open(out,"w") as f:
    f.write(HDR)
    for pfx,n,cig in groups:
        for i in range(int(n)):
            f.write("%s%d\t0\tchrT\t1000\t60\t%s\t*\t0\t0\t%s\t*\tXS:A:+\n"%(pfx,i+1,cig,seq(qlen(cig))))
```

```sh
P=100M500N100M302N100M            # clean: exons 1000-1099, 1600-1699, 2002-2101
C=100M500N100M2D300N100M          # same molecule, 2D+300N instead of 302N
T=100M498N2D100M2D300N100M        # 498N+2D instead of 500N, and 2D+300N
python3 mkreads.py r_nd_identical.sam nd:3:50M100N2D50M
python3 mkreads.py r_dn_identical.sam dn:3:50M2D100N50M
python3 mkreads.py r_dn_rescue.sam       p:6:$P c:4:$C
python3 mkreads.py r_dn_rescue_major.sam p:4:$P c:6:$C
python3 mkreads.py r_shift.sam           p:6:$P t:4:$T
python3 mkreads.py r_shift_major.sam     p:4:$P t:6:$T
# 8 reads use the acceptor at 1600, 6 the one at 1606
X=100M500N100M302N100M; Y=100M506N94M302N
python3 mkreads.py r_count_alt_ctl.sam x:8:$X y:3:${Y}100M y1:1:${Y}90M y2:1:${Y}80M y3:1:${Y}70M
K=100M504N2D94M302N
python3 mkreads.py r_count_alt.sam     x:8:$X y:3:${Y}100M k1:1:${K}90M k2:1:${K}80M k3:1:${K}70M
```

`guide.gtf` (one transcript, exons 1000-1051 and 1152-1201 on chrT, + strand) is
used with `r_dn_identical.sam`.

Run each with `stringtie -L -v -o out.gtf FILE` (add `-G guide.gtf` for the guide
case). "distinct" is from the `-v` bundle line; transcripts are given as
`longcov/cov: exons`.

```text
Case (defect)                              master 3.0.4                          expected
-----------------------------------------  ------------------------------------  ---------------------------
r_nd_identical (2)                         3 distinct                            1 distinct
r_dn_identical -G guide.gtf (3)            3 distinct                            1 distinct
r_dn_rescue, 4 of 10 with deletion (1)     6.0/10.0: 1000-1099 1600-1699 2002-   10.0/10.0: 1000-1099
                                                                                 1600-1699 2002-2101
r_dn_rescue_major, 6 of 10 (1)             6.0/11.3: 1000-1099 1600-1701 2002-   same as above
                                           (exon end off by 2, last exon cov 14)
r_shift, 4 of 10 (2)                       5 distinct; 6.0/10.0, exons correct   2 distinct; 10.0/10.0
r_shift_major, 6 of 10 (2)                 7 distinct; 6.0/11.3: .. 1600-1701    2 distinct; 10.0/10.0,
                                                                                 1600-1699
r_count_alt_ctl (control, no deletions)    8.0/15.9: 1000-1099 1600-1699 2002-   (already correct)
r_count_alt, 3 of the 6 with deletion (4)  two transcripts: 8.0/7.9 (acceptor    one transcript, as in the
                                           1600) and 6.0/5.9 (acceptor 1606)     control
```

In `r_count_alt` the junction to 1606 ends up with a count of 12 (3 + 3x3) instead
of 6, so it is no longer treated as an error of the 8-read junction, while the
same 14 reads without deletions give one transcript.

## How common these alignments are

Share of spliced alignments with a deletion next to an intron, first 400,000
records of two long-read BAMs:

```text
Data                                         D before N   D after N
-------------------------------------------  ----------   ---------
NA12878 direct RNA, chr1, minimap2 2.17        11.5%        12.4%
a second long-read sample                       0.05%        0.5%
```

## Proposed fixes

One commit per defect (branch `lr-junction-deletions`, not pushed yet):

1. `GSam.cpp`: assign `deljunc.start`/`deljunc.end` directly instead of using the
   `GSeg(s,e)` constructor.
2. `GSam.cpp`: add the `juncsdel` entry only with the first aligned bases after the
   intron, as `(prevdel,del)`; the `D` case no longer adds one.
3. `rlink.cpp`: `deljuncmatch()` takes the `GReadAlnData` and expects offsets 0 for
   a junction with `guide_match`, i.e. it compares with what `processRead()` would
   store for the incoming alignment.
4. `rlink.cpp`: in the two rescue branches add `rd.read_count` to `nreads` and the
   read's share (`rd.read_count/nreads` of the tentative junction) of `mm` and `nm`.

With all four, the cases above give the expected results.

## Effect on real data

Assemblies with each fix added in turn, compared with master (same options), scored
with gffcompare against the RefSeq annotation of the chromosome. "lines" is the
number of transcript structures that differ from master (added plus removed).

```text
Data, options                    Build      Distinct   Tx     Lines  IC sens  IC prec
-------------------------------  ---------  ---------  -----  -----  -------  -------
NA12878 dRNA chr1, -L            master      456,123   4,161      0    20.3     40.5
                                 fix 1       456,471   4,153     54    20.3     40.6
                                 fixes 1-2   445,880   4,169     76    20.3     40.4
                                 fixes 1-4   445,880   4,078    185    20.4     41.5
NA12878 dRNA chr1, -L -G         master      469,461   3,749      0    36.8     87.6
                                 fixes 1-2   469,345   3,743     62    36.8     87.8
                                 fixes 1-3   426,883   3,743     62    36.8     87.8
                                 fixes 1-4   426,883   3,743     62    36.8     87.8
NA12878 dRNA chr21, -L           master       54,594     526      0    18.4     34.7
                                 fixes 1-4    53,960     521     15    18.1     34.5
chr21 short + long, --mix        master    1,043,071     812      0    22.9     30.7
                                 fixes 1-4 1,042,425     809     61    22.9     30.8
```

The effect on accuracy is small in both directions. Fix 3 removes 9% of the stored
reads with `-G` on the chr1 data without changing any transcript structure (coverage
values change). Fix 4 changed nothing in the `-G` runs: the GTFs with fixes 1-3 and
1-4 are identical there.

`make test`: tests 6 and 7 (`--mix`) differ from the expected GTFs in coverage
values only (10 and 13 transcripts, same structures, largest cov change 0.23% and
3.5%), starting with fix 1. The other tests pass.

## Questions

- Is the description of the intended rescue above right, in particular that a
  deletion before an intron should be treated like one after it?
- Should fix 4 also be applied at `rlink.cpp:17788`?
- `gclib/GSam.cpp` is vendored; the first two fixes probably belong in gclib too.
- Whichever lands second, this or the #511 run index, needs a small adjustment:
  the index key of an incoming alignment must use the same offsets as
  `deljuncmatch()`.
