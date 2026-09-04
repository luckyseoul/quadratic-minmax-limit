# A strict incremental cut bound for p31 transverse atom rows

Date: 2026-09-04

**Status:** proved a GPU-friendly necessary cost for each actual transverse
row.  It strictly strengthens the scalar `l1` budget, updates exactly in
`O(p)` work after one row-cell mutation, and is evaluated through the public
479-edge graph API.  Zero cost is not sufficient for an integral atom
decomposition, so this does not close branch C or residual (ii).

## 1. Exact cut interval

Let

\[
 K(a,b;c)=e_{ab}-e_{ac}-e_{bc}
\]

be a compact atom, with `a,b,c` distinct, and let

\[
 A(a,b,c)=e_{ab}+e_{ac}+e_{bc}
\]

be an all-positive distinct-label triangle.  For every label subset `U`, a
direct three-vertex check gives

\[
 K(a,b;c)(\delta U)\in\{-2,0\},\qquad
 A(a,b,c)(\delta U)\in\{0,2\}.                           \tag{1}
\]

The compact value is `-2` exactly when the distinguished label `c` is on one
side and both positive labels are on the other.  The positive triangle value
is `2` exactly when `U` splits the triangle.

Therefore every row made from `h` positive triangles and `b` compact atoms
satisfies

\[
                   -2b\le x(\delta U)\le2h              \tag{2}
\]

for every `U`.  The live top rows have

\[
 (h,b)=
 \begin{cases}
 (0,P_D-3),&D\text{ hard},\\
 (6,Q_D-9),&D\text{ opposite}.
 \end{cases}                                             \tag{3}
\]

For a hard row, `x` in (2) is the normalized transverse row after adding
back its literal star `S_j`; for an opposite row it is the full transverse
row.

## 2. The 496-cut bank

Retain all 31 singleton cuts and all 465 two-label cuts.  If

\[
 d_u=\sum_{v\ne u}x_{uv},
\]

then

\[
 x(\delta\{u\})=d_u,
 \qquad
 x(\delta\{u,v\})=d_u+d_v-2x_{uv}.                       \tag{4}
\]

Thus (2) gives 496 integer interval checks per row.  In a compact-only row,
putting `n_u=-d_u/2` turns the two-label inequalities into the transparent
edge box

\[
             -(n_u+n_v)\le x_{uv}\le b-n_u-n_v.          \tag{5}
\]

These inequalities are stronger than signed degrees and scalar occurrence
budgets.

For example, at the live hard count `b=11`, start with

\[
 +02+13-01-24-45-35
\]

and add nine copies of `K(6,7;8)`.  The resulting vector has exactly
positive mass 11, negative mass 22, `l1=33`, edge sum `-11`, and nonpositive
even degrees summing to `-22`.  It passes the complete scalar budget and the
signed-degree projection.  But

\[
                       x(\delta\{0,1\})=2>0,             \tag{6}
\]

so it is not compact-decomposable.  At least two unit coefficient edits are
needed to repair (6).  This is a concrete strict advance over the previous
`l1` diagnostic.

## 3. Rigorous edit lower bound

Besides (2), the implementation retains the exact occurrence identities

\[
 \sum_e x_e=3h-b,
 \quad \sum_e x_e^+\le3h+b,
 \quad \sum_e x_e^-\le2b,
 \quad \lVert x\rVert_1\le3(h+b),                        \tag{7}
\]

and the fact that every signed degree is even.  A unit cell edit changes
each scalar in (7) by at most one, toggles two degree parities, and changes a
fixed cut by at most one.  Consequently the maximum of

- the four scalar defects in (7),
- half the number of odd signed degrees, and
- the largest singleton/two-label interval violation

is a rigorous lower bound on coefficient `l1` distance to every
atom-decomposable row.  Taking a maximum avoids double counting.  The code
also exposes the sum of all interval hinges as a smoother search energy; it
is an exact deterministic cost but is not labeled as an edit lower bound.

For opposite rows the degree projection is sharpened at negligible cost.
Write `d_v/2=A_v-C_v`, where `A_v` counts occurrences among the six positive
triangles and `C_v` counts distinguished compact occurrences.  Such degree
data are realizable at this projection exactly when the degrees are even and
one can choose

\[
 \max(0,d_v/2)\le A_v\le\min(6,b+d_v/2),
 \qquad\sum_vA_v=18.                                     \tag{8}
\]

The interval-sum test for (8) is included in the incremental search energy.
It is only a degree projection, not the full edge-level atom transport.

## 4. Exact incremental update

For one mutation

\[
                         x_{uv}\leftarrow x_{uv}+\Delta,
\]

only `d_u,d_v` change.  The two-label cut of `{u,v}` itself is unchanged:

\[
 (d_u+\Delta)+(d_v+\Delta)-2(x_{uv}+\Delta)
 =d_u+d_v-2x_{uv}.                                       \tag{9}
\]

For every `w` outside `{u,v}`, the two cuts `{u,w}` and `{v,w}` each change
by `Delta`; all remaining two-label cuts stay fixed.  Hence one update
touches exactly

```text
2 singleton cuts + 2(p-2) pair cuts = 2 + 58 at p=31.
```

The implementation updates the cell, degrees, four occurrence totals,
degree parity count, all affected cut hinges, and the complete additive
search cost in `O(p)`.  Focused tests compare a sequence of incremental
updates against fresh full reconstruction after every mutation.

## 5. Actual 479-edge graph replay and scope

The public `centered_physical_graph()` API with graph hash

```text
c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d
```

produces 16 hard and 16 opposite normalized rows.  All 32 have two-label cut
violations.  Their rigorous coefficient-edit lower bounds range from 122 to
194.  The old `l1` defect already dominates those large bounds on this bad
seed; the strict witness (6) proves that the new cut bank becomes genuinely
stronger as center search approaches the scalar budget surface.

A GPU center search should keep, per row, 465 cell values, 31 degrees, 465
two-label cuts, the positive/negative masses, odd-degree count, and additive
hinge sum.  A half-center replacement is a sparse batch of the cell updates
above; in its target hard row the literal-star correction is included in the
same batch.  The parallel profile is unchanged by nonzero center scaling.

This cost is a presolver and search objective.  A zero-cost candidate still
requires the exact integer compact/positive-triangle transport checker and
the simultaneous global conditions.

Replay with

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_row_atom_incremental_bound.py
```

Branch C, residual (ii), E1, and `L=1/2` remain **OPEN**.
