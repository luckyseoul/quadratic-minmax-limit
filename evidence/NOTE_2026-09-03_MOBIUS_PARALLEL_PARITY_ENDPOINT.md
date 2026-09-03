# The all-active Mobius parallel-parity endpoint obstruction

Date: 2026-09-03

Status: proved for every balanced branch-C prime \(p=4r+3\ge31\), under
the all-active ansatz of one localized Mobius half in every hard direction.
The result strengthens the bare support-size cancellation floor by zero, one,
or two cancellation units according to an exact residue class. It closes
broad \(j=0\) and \(j=1\) endpoint bands. Outside those bands it is only a
necessary condition; it does not solve the divided integral Boolean fibre or
close residual (ii).

## 1. Exact direction parity after every cancellation

Put

\[
 p=4r+3,\qquad m={p+1\over2}=2r+2,
\]

and assume all \(m\) hard centers are nonzero. Index the localized halves by
their distinct hard target directions \(L_i\), and write \(M_i\ne L_i\) for
their auxiliary projective directions.

The exact parallel vector of one half has

\[
 P_{L_i}=P_{M_i}=1,\quad P_{L_i-M_i}=0,\quad
 P_{L_i+qM_i}=1+\eta(1+q)
\]

in the remaining rows. Consequently its parallel vector modulo two is
exactly

\[
                         e_{L_i}+e_{M_i}.                 \tag{1}
\]

Let \(R_D\) count all raw half occurrences of spatial difference direction
\(D\), and let \(u_D\) count the actual nonzero ternary-support orbits in
that direction after every overlap and cancellation. At one orbit \(O\), if
\(n_O\) is its raw multiplicity and \(c_O\in\{-1,0,1\}\) its final
coefficient, then

\[
                         n_O-|c_O|\equiv0\pmod2.
\]

This remains true for triple and higher overlaps. Hence \(u_D\equiv R_D\)
in every direction. If \(c_D\) is the parity of the number of auxiliaries
\(M_i\) in direction \(D\), summing (1) gives

\[
             u\equiv {\bf1}_{\rm hard}+c\pmod2,
             \qquad \operatorname {wt}(c)\le m.          \tag{2}
\]

The inequality uses only that \(c\) is the parity vector of exactly \(m\)
auxiliary choices. In particular \(\operatorname{wt}(c)\) is even because
\(m\) is even.

## 2. Import the actual compact-atom quotas

The balanced hard compact counts \(e_L\) have total

\[
                              E=t+1
\]

and the actual hard parallel quotas are \(P_L=3+e_L\). The opposite quotas
\(Q_L\) are balanced with total

\[
                  10r+6+t=E+5m-5.                       \tag{3}
\]

Thus, with \(P\) denoting the full \(2m\)-direction quota vector, put

\[
 v=P\bmod2+{\bf1}_{\rm hard},\qquad w_0=\operatorname {wt}(v).
                                                                    \tag{4}
\]

On hard directions, \(v_L=e_L\bmod2\); on opposite directions,
\(v_L=Q_L\bmod2\). Write \(E=am+\alpha\), \(0\le\alpha<m\). Directly from
the two balanced allocations and (3),

\[
 w_0=
 \begin{cases}
  m+5,&\alpha\ge5, a\text{ even},\\
  m-5,&\alpha\ge5, a\text{ odd},\\
  m+2\alpha-5,&\alpha<5, a\text{ even},\\
  m+5-2\alpha,&\alpha<5, a\text{ odd}.
 \end{cases}                                             \tag{5}
\]

Equivalently, for

\[
                         s=E\bmod 2m=(t+1)\bmod(p+1),
\]

the complete table is

\[
\begin{array}{c|c}
s&w_0\\ \hline
5\le s\le m&m+5\\
s\in\{4,m+1\}&m+3\\
s\in\{3,m+2\}&m+1\\
s\in\{2,m+3\}&m-1\\
s\in\{1,m+4\}&m-3\\
\text{all remaining residues}&m-5.
\end{array}                                               \tag{6}
\]

This is a closed formula, not a prime or support census.

## 3. A general cancellation-offset lower bound

Let \(a=a(T_U)\) be the forced fixed-edge word, let \(f_D\) count its
selected fixed antipodal edges in spatial direction \(D\), and let \(n_D\)
be the number of selected unused double orbits in that direction. The exact
direction slices are

\[
                         P_D=u_D+f_D+2n_D.                \tag{7}
\]

Equations (2), (4), and (7) imply

\[
 f\bmod2=v+c,
 \qquad |a|=\sum_Df_D\ge\operatorname {wt}(v+c).
\]

The vector \(v\) has odd weight and \(c\) has even weight at most \(m\).
Therefore

\[
                 \boxed{|a(T_U)|\ge\mu(s):=
                        \max\{1,w_0(s)-m\}.}             \tag{8}
\]

The bound is sharp at the aggregate direction-parity level. If \(w_0>m\),
choose \(c\) on any \(m\)-subset of \(\operatorname{supp}v\). If
\(w_0\le m\), choose \(c\) on all but one element of that odd support and
use the remaining auxiliary choices in pairs. The restrictions
\(M_i\ne L_i\) do not improve this aggregate minimum: a derangement handles
the all-hard case, and opposite directions give unrestricted destinations.
This sharpness is only for the direction-parity relaxation, not for actual
Mobius intersections.

Let

\[
 \kappa_0=t_{\max}-t+1,
 \qquad \kappa=\kappa_0+j.
\]

The exact support ledger gives

\[
                         |H|-|U|=1+2j.                    \tag{9}
\]

The Hamming equation then gives \(|a|\le1+2j\). Combining this with (8)
proves the general strengthening

\[
 \boxed{j\ge j_{\rm par}(s):={\mu(s)-1\over2}=
 \begin{cases}
 2,&5\le s\le m,\\
 1,&s\in\{4,m+1\},\\
 0,&\text{otherwise}.
 \end{cases}}                                            \tag{10}
\]

Thus the cancellation floor is really
\(\kappa\ge t_{\max}-t+1+j_{\rm par}(s)\) before any transverse target
coordinate is inspected.

## 4. The first two endpoints

At \(j=0\), the Hamming equation would force one fixed edge and no divided
column. Equations (6)--(10) exclude this uniformly whenever

\[
                         \boxed{4\le s\le m+1.}           \tag{11}
\]

For the remaining residues, parallel parity does not prove
\(T_U=Aa(T_U)\); that equality remains the exact target-dependent test.

At \(j=1\), the two formal Hamming alternatives are three fixed edges and
no divided column, or one fixed edge and one divided column. Equation (10)
gives:

* if \(5\le s\le m\), both alternatives are impossible;
* if \(s\in\{4,m+1\}\), exactly three fixed edges and zero divided columns
  are forced, so the one-column branch is impossible;
* otherwise both alternatives remain target-dependent.

In the one-fixed/one-column branch, (7) determines the unique direction
\(D\) with \(n_D=1\). Its exact remaining test is

\[
                 \widehat T_U=\widehat B_O               \tag{12}
\]

for one retained nonfixed orbit \(O=([a],[\delta])\) parallel to \(D\).
Geometrically, a divided column has one entry in every target direction,
all with its Paley sign. It is \(P_D\) in direction \(D\). In every other
row \(L\), it is the quotient cell with

\[
                  \alpha_L=L(a)\pmod{\pm1},\qquad
                  \beta_L=L(\delta)^2.                   \tag{13}
\]

Thus (12) can be checked as two coherent projective square-evaluation
patterns, not by searching the whole Boolean cube. Equations (12)--(13) are
an exact criterion, not an existence theorem.

## 5. What grouped uncertainty adds at \(j=0\)

The affine-block transform has a useful interpretation for the forced word:
\((M^{\mathsf T}a)|_{B_L}\) is exactly the vector of nonparallel fixed-cell
residues of \(T_U\) in row \(L\). If \(g\) is the number of nonzero such row
groups, grouped uncertainty gives

\[
                              |a|+g\ge p+1=2m.            \tag{14}
\]

Hence a \(j=0\) survivor with \(|a|=1\) must have \(g=2m-1\): its fixed
residue is nonzero in every row except the unique direction parallel to that
fixed antipodal edge. This is exactly the affine-block signature of a
one-point word, so (14) alone cannot exclude it.

The present compact-residual theorem supplies centrality but does not bound
these fixed-cell parities. A compact atom \(K(v,-v;0)\) is already a central
example with an odd fixed-cell coefficient. Consequently no claim that the
actual compact/Mobius residue has two silent groups is currently justified.
Proving that target-specific statement would immediately exclude every
remaining \(j=0\) residue by (14), but it is a new lemma, not a consequence
of centrality.

## 6. Scope and reproduction

The proof applies only to the balanced all-active branch-C Mobius ansatz.
It does not apply unchanged to zero hard centers, another antisymmetric
preimage, unbalanced allocations, nonzero odd forms, or branch B. Passing
(10) does not construct a common graph.

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_mobius_parallel_parity_endpoint.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_mobius_parallel_parity_endpoint.py
```
