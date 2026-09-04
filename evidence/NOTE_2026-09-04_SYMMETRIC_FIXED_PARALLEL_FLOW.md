# Fixed/parallel network flow and transverse cycle steering

Date: 2026-09-04

**Status:** proved an exact all-prime network-flow normal form for the fixed
cells and parallel quotas of the fixed-edge-eliminated symmetric Boolean
fibre.  The projected constraint matrix is totally unimodular after an
arbitrary structured Mobius deletion.  Its fibre moves are alternating
cycles, and clean radial four-cycles span all mixed degree-six/eight global
moment corrections over `F_p` in the large-prime branch.

This does **not** prove that a conformal cycle sequence exists inside one
binary point, does not fix the pure radial marginals or all nonfixed
transverse cells, constructs no common simple graph, and does not close
residual (ii).

## 1. The divided symmetric box

Let `p=2h+1` be odd, let `d=p+1`, and retain the notation of
`NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md`.  After the uniquely
forced antipodal edge vector has been subtracted, the live variables are

\[
 b_O\in\{0,1\},\qquad O=([a],[\delta])\notin U,             \tag{1}
\]

where `U` is the actual support of the chosen antisymmetric Mobius lift.
The divided columns are

\[
 \widehat B_O=\left((B_O)_{\rm fix}/2,(B_O)_{\rm pair}\right).
                                                                    \tag{2}
\]

The fixed word `a(T_U)` is already determined by `U`; the variables in (1)
cannot change it.  They must still reproduce the divided fixed-cell values,
the exact parallel quotas, and every nonfixed paired-cell coordinate.

For a nonzero vector `v`, write `L_v` for its unique projective annihilator.
Put

\[
 A=L_a,\qquad D=L_\delta,
\]

and let `tau_delta` be the source sign, constant on the inversion orbit.  In
the Paley application it depends only on the projective difference direction
because multiplication of `delta` by a nonzero scalar changes its anisotropic
norm by a square.

## 2. Every fixed/parallel column is a network arc

**Column lemma.** In ordinary row coordinates, the fixed part of
`widehat B_O` is

\[
 \begin{cases}
 \tau_\delta e_{P_D},&A=D,\\[1mm]
 \tau_\delta\left(e_{P_D}+e_{K_A(0,A(\delta)^2)}\right),
     &A\ne D.
 \end{cases}                                                \tag{3}
\]

There are no other fixed entries.

**Proof.** In row `D`, the difference projects to zero, so both edges in the
central source pair are parallel.  Their common coefficient is
`2 tau_delta`, which becomes `tau_delta` after division.  If `A!=D`, then in
row `A` the midpoint projects to zero and the endpoints project to
`+-A(delta)`.  Both inverted edges therefore hit the same fixed antipodal
cell, again with divided coefficient `tau_delta`.  In every other row both
the midpoint and the difference are nonzero, so the image belongs to a
nonfixed target-cell pair.  When `A=D`, the second displayed fixed cell is
the already-counted parallel cell.  This proves (3).  `square`

Negate each parallel row.  Equation (3) becomes

\[
 \boxed{\tau_\delta
 \left(e_{K_A(0,A(\delta)^2)}-e_{P_D}\right)}               \tag{4}
\]

when `A!=D`.  This is the incidence column of an oriented edge between the
parallel node `P_D` and the fixed-cell node `K_A(0,A(delta)^2)`; changing
`tau_delta` only reverses its orientation.  When `A=D`, adjoin one root node
and regard the one-entry column as an edge from `P_D` to the root, whose row
is deleted.

**Network-flow theorem.** The projection of (1)--(2) onto all divided fixed
cells and all parallel rows is a capacitated directed-incidence system.
Consequently:

1. its matrix is totally unimodular;
2. deleting every column in an arbitrary used set `U` preserves total
   unimodularity;
3. adding the bounds `0<=b_O<=1` leaves an integral polytope for every
   integral target; and
4. projected Boolean feasibility is exactly a capacitated max-flow/min-cut
   problem.

The total Hamming equation is redundant here because summing all exact
parallel quotas already counts every chosen double orbit once.

For completeness, total unimodularity needs no external classification.
Every column of (4) has one `+1` and one `-1`, or only one nonzero after the
root row is deleted.  Expanding a square subdeterminant along any column with
at most one retained nonzero reduces its order.  If every column has two
retained nonzeros, the row sum is zero and the determinant vanishes.  By
induction every square subdeterminant is `0,+1,-1`.

This theorem is stronger than mod-two surjectivity for this projection: it
gives an actual zero-one solution whenever the bounded real flow is
feasible.  It says nothing yet about the omitted nonfixed transverse rows.

## 3. The graph object and its alternating cycles

The flow network gives the promised graph interpretation.  A Boolean
completion is a set of capacity-one arcs.  Two completions with the same
fixed and parallel data differ by an integral circulation in the augmented
network, hence by alternating cycles.  Toggling one conformal alternating
cycle preserves

* the already-forced fixed-edge word;
* the number of selected double orbits;
* every divided fixed-cell value; and
* every projective parallel quota.

It can change only the nonfixed transverse cells.

There is a particularly explicit family.  Fix projective midpoint and
difference directions `(A,D)`, choose representatives `a_0` and `delta_0`,
and write

\[
 a=x a_0,\quad \delta=y\delta_0,\qquad
 \xi=x^2,\quad\eta=y^2\in Q=(\mathbf F_p^*)^2.             \tag{5}
\]

For the rectangular formulas below, assume that the source sign is constant
on a projective difference direction.  This is exactly true for the branch-C
Paley sign; it is not asserted for an arbitrary inversion-invariant signing.
The `h^2` variables in this block form an `h` by `h` radial table.  For two
midpoint radii and two difference radii, the move

\[
 E_{\xi_1,\eta_1}+E_{\xi_2,\eta_2}
 -E_{\xi_1,\eta_2}-E_{\xi_2,\eta_1}                       \tag{6}
\]

has zero fixed/parallel image.  The source sign is constant on the
projective difference direction `D`, so it does not spoil the cancellation.
The basepoint moves

\[
 (e_\xi-e_{\xi_0})(e_\eta-e_{\eta_0})^{\mathsf T}
\]

form an integral basis of all radial tables with zero row and column sums.

An even smaller legal operation is a parallel-arc exchange: for fixed
`[delta]`, replace a selected `([a_1],[delta])` by an unselected
`([a_2],[delta])` with `a_1` and `a_2` in the same projective midpoint
direction.  The two arcs have identical network endpoints.  This swap
removes one central pair of physical graph edges and adds another while
preserving all data listed above.

If used orbits puncture the radial table, the legal moves are the alternating
even cycles of the remaining bipartite graph.  Four-cycles generate the
unpunctured table lattice, but they need not connect a punctured binary fibre:
a chordless longer cycle may be indispensable.  Thus the theorem does not
silently assume a complete Markov basis after deletion.

## 4. Exact action on transverse cells and even moments

For a row `L` distinct from `A,D`, put

\[
 \alpha=L(a_0)^2\xi,\qquad
 \beta=L(\delta_0)^2\eta.                                \tag{7}
\]

Both leading factors are nonzero squares, so (7) independently permutes the
rows and columns of the radial table.  The move (6) therefore appears in
every transverse projection as the same alternating cell rectangle after
two multiplicative relabelings.  It is not a full edge--Radon kernel move;
it is an exact way to steer the still-unmatched transverse target while the
solved projection remains frozen.

For projected endpoints `s=x-y,t=x+y`, the moment channel from the compact
ray notes is

\[
 Q_{2n,k}(s,t)
 =(s-t)^2(st)^k(s+t)^{2n-2-2k}
 =2^{2n-2k}\beta\alpha^{n-1-k}(\alpha-\beta)^k.           \tag{8}
\]

Thus the change of a mixed monomial under (6) factors exactly as

\[
 (\alpha_1^r-\alpha_2^r)(\beta_1^s-\beta_2^s).           \tag{9}
\]

Pure-alpha or pure-beta statistics vanish in (9); mixed moments are exactly
what the rectangle can change.

There is no hidden radial rank loss.  The functions

\[
 1,\xi,\ldots,\xi^{h-1}
\]

are a Vandermonde basis on the `h` nonzero squares.  Therefore
`xi^r`, `1<=r<=h-1`, form the dual basis of the zero-sum radial functions.
Tensoring the two statements shows that

\[
 \xi^r\eta^s,\qquad1\le r,s\le h-1,                      \tag{10}
\]

are a dual basis of the four-cycle space modulo `p`.  Hence clean
four-cycles span every mixed radial correction.

For fixed `n`, the first `n-1` polynomials in (8), with
`0<=k<=n-2`, are triangular in

\[
 \alpha^{n-1}\beta,\alpha^{n-2}\beta^2,
 \ldots,\alpha\beta^{n-1}.                               \tag{11}
\]

Their diagonal coefficients are nonzero.  The final channel `k=n-1` adds
the pure `beta^n` marginal.  Consequently four-cycles control the complete
mixed part of every even degree through `p-1`; they deliberately do not
claim control of the pure marginal.

## 5. Clean blocks after arbitrary Mobius deletion

Now specialize to the all-active large-prime branch, where

\[
 |U|\le {p^2-1\over2}=dh.                                \tag{12}
\]

There are `d^2` projective blocks `(A,D)`, and each used orbit contaminates
only one.  Thus at least

\[
 C=d^2-dh={(p+1)(p+3)\over2}                             \tag{13}
\]

blocks are completely clean.

Fix positive exponents `a,b<d`.  A linear functional on binary forms of
degree `a+b`, evaluated on `A^aD^b`, is a bihomogeneous polynomial
`F(A,D)` of bidegree `(a,b)` on
`P^1(F_p) x P^1(F_p)`.  If `F` is nonzero, at most `a` fibres in `A` can
vanish identically; every other fibre has at most `b` zeros.  Hence

\[
 |Z(F)|\le ad+(d-a)b=ad+bd-ab.                            \tag{14}
\]

The full products `A^aD^b` span all degree-`a+b` binary forms: on the
diagonal they include the pure powers `A^(a+b)`, and those span because the
degree is below `p+1`.  Therefore any subset larger than (14) still spans.

For degree six the mixed bidegrees are `(4,2),(2,4)`, and

\[
 C-(6d-8)={(p-1)(p-7)\over2}>0.                          \tag{15}
\]

For degree eight they are `(6,2),(4,4),(2,6)`, and

\[
 C-(8d-12)={(p-1)(p-11)\over2}>0.                        \tag{16}
\]

Thus for every branch prime `p>=31`, the clean projective blocks span every
mixed degree-six and degree-eight global binary form.  Combining this with
the independent radial control (10) proves:

\[
 \boxed{\text{Clean four-cycle lattices surject modulo }p
 \text{ onto the full mixed degree-six/eight correction space.}} \tag{17}
\]

This is a source-graph steering theorem respecting the forced fixed word and
all parallel quotas.  It is not another abstract moment-dominance statement:
the moves are explicit symmetric edge-pair exchanges.  Conversely, (17) is
still linear.  Its coefficients need not decompose into a sequence of
conformal toggles that remains in `{0,1}` at every step.

## 6. Exact boundary and the top endpoint

The new reduction retires the possibility of a hidden integrality gap in the
fixed/parallel projection.  It also proves that degree-six/eight mixed forms
are not missing from the legal cycle lattice.  The unresolved data are now
sharply separated:

1. max-flow cut feasibility for the projected target after the actual `U`;
2. conformal availability of the required alternating cycles inside a
   binary flow;
3. the pure radial moment marginals; and
4. the complete nonfixed transverse cell vector, which is stronger than its
   degree-six/eight contractions.

At the top endpoint, if

\[
 |H|-|U|-|a(T_U)|=0,                                     \tag{18}
\]

the divided Hamming equation forces `b=0`.  There is then no completion-side
flow or alternating cycle to optimize.  A failing top-end graph can only be
repaired by changing the antisymmetric Mobius/paired-SDR skeleton itself.
Changing the orientation gauge while leaving `U` fixed cannot do this,
because the forced symmetric pair-total support is unchanged.

The current `p=31` top paired-SDR witness has `|U|=478`, one forced fixed
edge, and `|H|=479`, so it lies exactly in (18).  This observation is a scope
guard, not a finite-prime exclusion theorem.

## 7. Reproduction

The verifier replays all 576 nonfixed symmetric columns at `p=7`.  It checks
that every divided fixed column has precisely the signed-incidence support
in (3), obtains 72 root half-edge columns and 504 two-endpoint arc columns,
and verifies an explicit radial four-cycle with zero fixed projection and
nonzero transverse projection.  This is a fail-when-wrong replay of the
symbolic proof, not a prime or residual-configuration census.

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_symmetric_fixed_parallel_flow.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_symmetric_fixed_parallel_flow.py
```

No parity, Smith, least-norm, coefficient-`l1`, small-prime, or paired-SDR
census is used.  The common graph, residual (ii), E1, and `L=1/2` remain
**OPEN**.
