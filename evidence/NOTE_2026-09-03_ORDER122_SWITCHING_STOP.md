# Order-122 conference switching search: exact stopping record

Date: 2026-09-03

Status: no nonregularizable symmetric conference matrix was found. These
computations close only the named switching ansatzes; they do not classify
all order-122 conference matrices and do not close residual (ii).

## Exact completed searches

`pn_wqh_lines.cpp` enumerates every balanced six-direction subset of the 12
directions of `AG(2,11)` and every unordered pair of affine lines. It checks
the Wang--Qiu--Hu hypotheses exactly and validates every prospective output
against the `srg(121,60,29,30)` identity. The replay returned

```text
NONE tested=7988904 wqh_nontrivial=0
```

`pn_wqh3_search.cpp` enumerates type `(3,3,115)` WQH pairs. The four
`PGL(2,11)` representatives of balanced direction sets returned

```text
mask=63  NONE candidates=1520760
mask=95  NONE candidates=1774740
mask=111 NONE candidates=1529760
mask=119 NONE candidates=1709400
```

The total is 6,534,660 exactly checked candidate pairs. The commands were

```text
g++ -O3 -std=c++20 evidence/order122_switching/pn_wqh3_search.cpp -o /tmp/pn_wqh3
for m in 63 95 111 119; do /tmp/pn_wqh3 $m; done
g++ -O3 -std=c++20 evidence/order122_switching/pn_wqh_lines.cpp -o /tmp/pn_wqh_lines
/tmp/pn_wqh_lines
```

`pn_wqh4_linesubsets.cpp` exhausts pairs of four-point subsets, each lying
on an affine line. Translation fixes one selected vertex at zero; the four
balanced direction-set orbit representatives then returned

```text
mask=63  exact_candidates=8  NONE
mask=95  exact_candidates=0  NONE
mask=111 exact_candidates=0  NONE
mask=119 exact_candidates=1  NONE
```

Every prefilter survivor was checked against the full WQH conditions,
regularity, and the exact `srg(121,60,29,30)` identity. Thus this entire
line-supported WQH `(4,4)` ansatz is empty.

An independent replay, including ASan/UBSan, returned the same four rows.
An independent `PGL(2,11)` enumeration has exactly four orbits on the 924
balanced direction masks, of sizes `330,264,220,110` and with minima
`63,95,111,119`. The 43,560 collinear four-sets are unique; translation puts
zero in one selected set. A valid pair has identical neighbor-count
signatures away from its eight vertices, so at least two of the ten index
blocks are clean and the 45 block-pair indices must retrieve it. The eight
mask-63 and one mask-119 `exact_candidates` are prefilter candidates, not
valid WQH pairs: six overlap, and all three disjoint/internal-valid pairs
fail the full outside-vertex condition.

The three disjoint rejection certificates are concrete.  For mask `63`,
`{0,36,84,120}` versus `{32,47,73,88}` fails at vertex `8` with
counts `(2,1)`, while `{0,62,70,113}` versus `{35,54,97,105}` fails at
vertex `2` with counts `(2,3)`.  For mask `119`, `{0,3,5,8}` versus
`{11,14,17,19}` fails at vertex `1` with counts `(0,2)`.

The independent CP-SAT model `wqh121_general.py` removes the line-support
restriction. For each of the Paley and Peisert graphs it fixes one vertex by
Cayley transitivity and exhausts every internal degree `0..3` and every
cross degree `0..4`. All 40 cells returned `INFEASIBLE`, generally in under
three seconds. Consequently neither starting graph has any nontrivial WQH
partition of type `(4,4,113)` satisfying the standard sufficient theorem.
The Boolean eigenshell verifier now fixes only the global sign coordinate;
it deliberately makes no second-coordinate symmetry assumption after a
switch.

## Bounded GM scan

`gm_oa121.py` searches ordinary Godsil--McKay sets in the representative
six-direction OA graph. The split mesh run covered requested even sizes
`4..60` with a two-second CP-SAT budget per `(size,internal degree)` cell.
It returned 25 nontrivial regular switched candidates. Exact Boolean
eigenshell models returned `OPTIMAL` for both signs on every candidate, so
none meets Proposition 15.762's counterexample criterion. Because individual
CP-SAT cells can time out, this is a bounded candidate scan, not an
exhaustive nonexistence certificate.

## Source hashes

```text
4a2d33a5cdff484a9a54d5faa25b54961f67b57200854fda06494e791ba4fa24  gm_oa121.py
a9001052233f6b9d0aef67de1e58263e660dda234608bb126b342c2861872df8  gm121_scan.py
2ae21c03bbe1421fe019ae714fc3e895c2b8e4544a50358044927a1f5d2e639c  peisert121_exact.py
999fb0b45d18487adfc5712134d4a45bf541c5e25b770beada399a41524d4d64  pn_wqh3_search.cpp
747fcc05e6e92499d585eef33dd433ab94cdf654d1db1f9b22f986ec54af2b51  pn_wqh4_linesubsets.cpp
6585819a7c9642b4da0989784cb7e39bcfc7d8dafa4c310c64091d0bb5af3e5c  pn_wqh_lines.cpp
1df403a9be9d088439cf632da03f03cccc4f66cdd71e1bf8ff3169b30b105a19  wqh121_general.py
```

## Exact boundary

The Paley/Peisert/linear-OA switching ansatzes above remain regularizable.
The live conference route is an explicit square-order symmetric conference
class outside the audited families with no Boolean `+p` or `-p` eigenvector.
Proposition 15.762 then supplies the full cube-norm gap without a near-shell
optimizer. Residual (ii), E1, `L=1/2`, and the original MO limit remain OPEN.
