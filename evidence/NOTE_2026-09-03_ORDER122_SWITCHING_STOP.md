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

## Bounded GM scan

`gm_oa121.py` searches ordinary Godsil--McKay sets in the representative
six-direction OA graph. The split mesh run covered requested even sizes
`4..60` with a two-second CP-SAT budget per `(size,internal degree)` cell.
It returned 25 nontrivial regular switched candidates. Exact Boolean
eigenshell models returned `OPTIMAL` for both signs on every candidate, so
none meets Proposition 15.762's counterexample criterion. Because individual
CP-SAT cells can time out, this is a bounded candidate scan, not an
exhaustive nonexistence certificate.

## Explicit unfinished branch

`pn_wqh4_linesubsets.cpp` compiles and encodes pairs of four-point subsets
lying on affine lines. Its full enumeration was not completed in this
checkpoint, and no mathematical result is claimed from that source.

## Source hashes

```text
4a2d33a5cdff484a9a54d5faa25b54961f67b57200854fda06494e791ba4fa24  gm_oa121.py
999fb0b45d18487adfc5712134d4a45bf541c5e25b770beada399a41524d4d64  pn_wqh3_search.cpp
7ef82ea5fa6924d98e5456b7cfabd62c324ca195930836986822bec363d5720d  pn_wqh4_linesubsets.cpp
6585819a7c9642b4da0989784cb7e39bcfc7d8dafa4c310c64091d0bb5af3e5c  pn_wqh_lines.cpp
```

## Exact boundary

The Paley/Peisert/linear-OA switching ansatzes above remain regularizable.
The live conference route is an explicit square-order symmetric conference
class outside the audited families with no Boolean `+p` or `-p` eigenvector.
Proposition 15.762 then supplies the full cube-norm gap without a near-shell
optimizer. Residual (ii), E1, `L=1/2`, and the original MO limit remain OPEN.
