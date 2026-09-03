# The symmetric image of a Mobius half and the exact overlap gate

**Status:** exact all-prime formulas and a sharp correction to the apparent
one-edge obstruction.  No finite prime census is used.  The formulas reduce
the live problem to a restricted symmetric Boolean fibre; they do not solve
that fibre and do not close residual (ii).

This note continues
`NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md`.  There, for an odd prime
\(p\), independent functionals \(L,M\), and \(j\ne0\), the half

\[
 E=\{e_t=\{u_t,v_t\}:t\in T\},\qquad
 T=\mathbf F_p\setminus\{-1\},
\]

is defined in \((L,M)\)-coordinates by

\[
 u_t=j\left(1,{t\over t+1}\right),\qquad v_t=j(t,t).       \tag{1}
\]

The chain \(z={\bf1}_{-E}-{\bf1}_E\) has ordinary edge--Radon image
\(A_j=S_{-j}-S_j\) in row \(L\) and zero in every other row.  It has
exactly \(p-1\) nonfixed inversion orbits.

## 1. Ordinary image of one half

In row \(L\), \(t=1\) is parallel and the excluded value \(t=-1\) would
have supplied the antipodal star edge.  Therefore

\[
 P_L(E)=1,\qquad K_L(E)=S_j-\{j,-j\}.                     \tag{2}
\]

In row \(M\), the projected endpoints \(jt/(t+1),jt\) are equal only for
\(t=0\), so \(P_M(E)=1\).  Write any remaining row as \(N=L+mM\).  Its
endpoint labels, after division by \(j\), are

\[
 A_t={(1+m)t+1\over t+1},\qquad B_t=(1+m)t.               \tag{3}
\]

For \(m=-1\) they are never equal.  Otherwise parallelism is
\((1+m)t^2=1\).  The value \(t=-1\) is not a solution because here
\(m\ne0\).  Hence, with \(\eta\) the quadratic character,

\[
 \boxed{
 P_N(E)=
 \begin{cases}
 1,&N=L,\\
 1,&N=M,\\
 0,&N=L-M,\\
 1+\eta(1+m),&N=L+mM,\quad m\ne0,-1.
 \end{cases}}                                             \tag{4}
\]

The Mobius involutions from the preceding note pair every transverse cell
in every \(N\ne L\) with its central negative.  Fixed parameters give
centrally fixed cells.  Thus

\[
 K_N(E)=I K_N(E)\quad(N\ne L).                            \tag{5}
\]

Finally, summing (4), or simply assigning each edge its unique parallel
direction, gives

\[
                       \sum_NP_N(E)=p-1.                  \tag{6}
\]

These statements describe \(E\), not the difference \(E-(-E)\).

## 2. The exact Paley pattern

Now assume \(p\equiv3\pmod4\), as in branch C.  Use coordinates in which
the anisotropic norm is \(Q(x,y)=x^2+y^2\).  If \(e_1,e_2\) are the primal
basis dual to \(L,M\), then

\[
 d_t:=u_t-v_t
 =j\left((1-t)e_1-{t^2\over t+1}e_2\right),
\]

and therefore

\[
 {t+1\over j}d_t=e_1-t^2(e_1+e_2).                       \tag{7}
\]

The scalar removed in (7) is a square in \(Q\).  Consequently the Paley
column sign is exactly

\[
 \boxed{\tau_t=\eta\!\left(Q(e_1-t^2(e_1+e_2))\right).}   \tag{8}
\]

In particular \(\tau_t=\tau_{-t}\) whenever both parameters belong to
\(T\).  The only unpaired exceptional value is \(t=1\), because \(-1\)
was deleted.  Moreover

\[
 \tau_1=\eta(Q(e_2))=\epsilon_L.                          \tag{9}
\]

Put \(S_{L,M}=\sum_{t\in T}\tau_t\).  Completing the missing value gives
the exact quartic-character trace

\[
 \boxed{
 S_{L,M}=\sum_{t\in\mathbf F_p}
   \eta\!\left(Q(e_1-t^2(e_1+e_2))\right)-\epsilon_L.}    \tag{10}
\]

This trace depends on the auxiliary direction; it is not replaced here by
a constant or by a finite-prime table.

## 3. The forced symmetric source and its row invariant

The normalized residual map is
\(\widetilde R_N(h)=\epsilon_NR_N(\tau h)\).  The physical ternary
difference realizing \(A_j\) is \(w=\tau z\).  On each orbit of \(E\), its
unique selected graph edge is \(e_t\) when \(\tau_t=-1\), and \(-e_t\)
when \(\tau_t=+1\).  Thus the normalized signed source of that selected
half is

\[
 h_E={1\over2}\left[z+q_E\right],\qquad
 q_E=\sum_{t\in T}\tau_t({\bf e}_{e_t}+{\bf e}_{-e_t}).   \tag{11}
\]

Here \(q_E\) is exactly the forced inversion-symmetric pair-total chain.
Its normalized Radon image is central in every row.  Since a parallel edge
has \(\tau_t=\epsilon_N\), equations (4) and (11) give

\[
 \widetilde P_N(h_E)=P_N(E),\qquad
 \widetilde P_N(q_E)=2P_N(E).                             \tag{12}
\]

Every source edge has one cell in every row.  Therefore the full signed
row sums are the direction-independent trace, up to the row sign:

\[
 \boxed{
 \sum\widetilde R_N(h_E)=\epsilon_NS_{L,M},\qquad
 \sum\widetilde R_N(q_E)=2\epsilon_NS_{L,M}.}             \tag{13}
\]

At the chosen hard row \(\epsilon_L=1\), the transverse part of the
selected half can also be read explicitly.  Starting from
\(-[S_j-\{j,-j\}]\), every parameter with \(\tau_t=+1\) replaces the
negative cell \(\{j,jt\}\) by the positive central mate
\(\{-j,-jt\}\), equivalently adding their central pair.  In all other rows
the selected transverse chain is central.  Equations (12)--(13) are exact
necessary data for any completion by unused double orbits and fixed
antipodal edges.

## 4. Two-trade overlaps

For two signed localized trades, inspect one inversion orbit.  If it occurs
in only one trade, its coefficient remains \(\pm1\).  If the two
orientations agree, the coefficient becomes \(\pm2\), forbidden by the
ternary box.  If they disagree, the two coefficients cancel.  Thus a
pairwise overlap compatible with ternarity must be an opposite-sign
overlap, and each such orbit deletes two nonzero support occurrences.

There is always one exact cancellation for two distinct target directions
\(L_1,L_2\) and arbitrary nonzero centers \(j_1,j_2\).  Put

\[
 X={L_1\over j_1},\qquad Y={L_2\over j_2},\qquad
 M_1=j_1(X+Y),\qquad M_2=j_2(X+Y).                        \tag{14}
\]

These auxiliaries are independent of their corresponding target
directions.  In \((X,Y)\)-coordinates the two halves become

\[
\begin{array}{ll}
 E_1(t):&\{(1,-1/(t+1)),(t,0)\},\\
 E_2(s):&\{(-1/(s+1),1),(0,s)\}.
\end{array}                                               \tag{15}
\]

The \(t=s=0\) edges are \(\{(1,-1),0\}\) and its negative, so their trade
coefficients cancel.  These are the only common inversion orbit: the
second endpoint \((t,0)\) of an \(E_1\)-edge cannot equal, up to sign, the
first \(E_2\)-endpoint because its second coordinate is \(0\ne\pm1\); it
can equal the second endpoint \((0,s)\) only when \(t=s=0\).  Hence the
sum is still ternary, has the two isolated star targets, and uses exactly

\[
                         2(p-1)-2                         \tag{16}
\]

nonzero inversion orbits.

More generally, an origin orbit \(\{0,u\}\) is shared by trades indexed by
\(i\in I\) precisely when signs \(\sigma_i\in\{\pm1\}\) can be chosen with

\[
                         L_i(u)=\sigma_i j_i.              \tag{17}
\]

The coefficient on one orientation is, up to a common sign,
\(\sum_i\sigma_i\), so ternarity is equivalent to
\(|\sum_i\sigma_i|\le1\).  Any two directions admit (17), because two
independent linear equations determine \(u\).  Three or more directions
impose genuine coherence on the centers; arbitrary centers need not meet
it.  This is the exact origin-sharing classification, not a claim that
non-origin overlaps are impossible in other auxiliary choices.

## 5. The apparent one-edge gap is not a contradiction

For \(p=4r+3\), the balanced branch-C ray has

\[
 t_{\min}=2r^2-4r-2,\quad t_{\max}=4r^2-2r-5,
 \quad |H|=4p+2t+1.                                      \tag{18}
\]

At its upper endpoint,

\[
 |H|_{\max}=8r^2+12r+3.                                  \tag{19}
\]

If all \(m=(p+1)/2=2r+2\) hard centers are nonzero, the deliberately
disjoint construction in the preceding note uses

\[
 m(p-1)=8r^2+12r+4=|H|_{\max}+1.                         \tag{20}
\]

Equation (16), however, changes this to

\[
m(p-1)-2=|H|_{\max}-1                                  \tag{21}
\]

for one paired cancellation, before the other trades are chosen disjointly
by the existing greedy avoidance argument.  Therefore (20) is a property
of one convenient lift, not a universal lower bound.  It cannot exclude
the all-nonzero-center case.

There is nevertheless a necessary construction budget across the whole
ray which must not be hidden by the upper-end calculation.  Put

\[
 N=m(p-1),\qquad |U|=N-2\kappa,
\]

where \(\kappa\) is the total opposite-sign cancellation count in a ternary
sum of the \(m\) halves.  Since

\[
 N-|H|=2(t_{\max}-t)+1,
\]

even the first size condition \(|U|\le |H|\) requires

\[
 \boxed{\kappa\ge t_{\max}-t+1.}                         \tag{21a}
\]

At equality the unused edge capacity is exactly one, so the fixed-edge
elimination theorem forces one antipodal fixed edge and no unused double
orbit.  More generally that remaining capacity is
\(2(\kappa-(t_{\max}-t))-1\), and the forced fixed-edge weight is odd.
For \(p=31\), the required cancellation count falls from \(110\) at
\(t=68\) to \(1\) at \(t=177\).  Thus the disjoint lift is not itself
extendable anywhere on this all-active ray.  Equation (16) constructs only
the final one-cancellation endpoint budget; no compatible multi-overlap
family attaining (21a) has been constructed.  This is a corrected
construction gate, not a capacity contradiction for arbitrary
antisymmetric lifts.

## 6. Exact remaining obstruction

The antisymmetric ternary box is solved, but its chosen orbit totals force
the central chain \(q_E\).  A full graph must correct the sum of these
chains to the prescribed hard compact and opposite AE-plus-compact central
rows using only:

* pair total \(0\) or \(2\) on every unused nonfixed inversion orbit;
* the fixed antipodal edges with binary coefficients; and
* any deliberately cancelled or shared trade orbits, subject to the local
  coefficient-imbalance rule above.

The exact parallel vector (4), quartic trace (10), and row-sum invariant
(13) are the current symbolic data for that restricted fibre.  No theorem
here proves that its compact-atom rows are attainable or contradictory.
Residual (ii), E1, \(L=1/2\), and the original MathOverflow problem remain
OPEN.

## Reproduction

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
      tests/test_mobius_half_symmetric.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_mobius_half_symmetric.py
