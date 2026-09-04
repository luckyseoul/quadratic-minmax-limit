# Exact fractional obstruction to a fixed global semimetric bank

Date: 2026-09-04

**Status:** at the `p=31,t=177,j=0` hard-fixed endpoint, an exact rational
fractional graph obeys the full parallel/fixed/no-double ledger and lies in
the real atom cone in every one of the 32 transverse rows.  Therefore no
fixed sum of row semimetric inequalities separates every graph allowed by
this ledger relaxation.  Adaptive row separation, integrality, and full atom
synchronization remain available; no common graph or residual-(ii) witness
is constructed.

## 1. The hard-fixed fractional graph

Use the frozen labelled top profile

\[
 (P_D)_{D\in\mathbf P^1(\mathbb F_{31})}
 =(15,14,14,15,16,15,16,14,16,14,14,16,16,15,15,14,
   14,16,16,16,16,14,14,16,14,16,16,16,14,14,14,14).
                                                               \tag{1}
\]

Its hard entries are `14^14,15^2`, its opposite entries are `15^3,16^13`,
and its sum is 479.  Freeze the sole fixed antipodal edge in hard direction
`F=1`, whose quota is 14.  Average it uniformly over its 15 magnitudes.
In each direction `D`, average the remaining `P_D-1_{D=F}` selections
uniformly over all 14,400 nonfixed actual edges of that spatial direction.

Equivalently, each of the 7,200 nonfixed inversion-orbit pairs in direction
`D` has total fractional occupancy

\[
                   {P_D-1_{D=F}\over7200}<1.                \tag{2}
\]

The product of the corresponding partition-matroid base polytopes is
integral.  Hence this point is a convex combination of actual 479-edge
simple graphs with profile (1), exactly one fixed edge in `F`, and no doubled
nonfixed orbit.  Independently average each hard literal center uniformly
over `F_31^*`.

This convex average is enough to refute a fixed linear separator: if one
fixed sum of row semimetric inequalities were positive on every allowed
graph and every nonzero center tuple, it would also be positive on their
average.

## 2. The three invariant edge categories

In any projected label row, the average has only three cell values:

\[
 A=\{\{0,a\}:a\ne0\},\quad |A|=30,
\]

\[
 B=\{\{a,-a\}:a\ne0\},\quad |B|=15,
\]

and the 420 remaining nonzero, nonantipodal pairs `C`.

Let `epsilon_D` be the Paley sign, put `f_D=1_{D=F}` and
`r_D=P_D-f_D`, and, for the observed row `N`, set

\[
 U_N=\sum_{D\ne N}\epsilon_Dr_D,
 \qquad V_N=\sum_{D\ne N}\epsilon_Df_D.                  \tag{3}
\]

For `D != N`, each generic projected cell has 31 preimages among all
`D`-parallel edges.  An antipodal cell has 30 nonfixed preimages and one
fixed preimage.  Thus the graph part of the expected row is

\[
 x_A=x_C={31\epsilon_NU_N\over14400},\qquad
 x_B={30\epsilon_NU_N\over14400}+{\epsilon_NV_N\over15}. \tag{4}
\]

For a hard row, subtracting the target literal star means adding `S_j`.
The uniform nonzero center hits an `A`-edge once and a `B`- or `C`-edge
twice, so add respectively `1/30,2/30,2/30` to (4).

There are five resulting row types.  Their exact `(A,B,C)` coefficients are

\[
\begin{array}{c|c}
\text{row type}&(x_A,x_B,x_C)\\ \hline
\text{hard }q=15&(-853/14400,\ 7/160,\ -373/14400)\\
\text{hard }q=14,\ N=F&(-791/14400,\ -3/160,\ -311/14400)\\
\text{hard }q=14,\ N\ne F&(-137/2400,\ 11/240,\ -19/800)\\
\text{opposite }q=16&(31/1200,\ -1/24,\ 31/1200)\\
\text{opposite }q=15&(403/14400,\ -19/480,\ 403/14400).
\end{array}                                                \tag{5}
\]

Their total edge sums are exactly `-b` in a hard row and `18-b` in an
opposite row, with `b=q-3` and `b=q-9`, respectively.

## 3. Exact fractional atom decompositions

For a compact atom `K(a,b;c)=e_ab-e_ac-e_bc`, record the sum of its
coefficients in `(A,B,C)`.  For an all-positive triangle, record the three
positive incidences in the same way.  Uniform averaging over every atom of
one type gives coefficient `v_i/|i|` on each edge in category `i`.
The signed-permutation group of the 15 nonzero antipodal label pairs fixes
zero, preserves each type class, and is transitive on each of `A,B,C`, so
this averaging identity is exact.

The compact type counts are

```text
(0,0,-1): 11760     (0,-1,0): 870      (-2,0,1): 420
(0,1,-2):   420     (-2,1,0):  15
```

and the positive-triangle type counts are

```text
(0,0,3): 3640       (2,0,1): 420       (0,1,2): 420
(2,1,0):   15
```

Here are nonnegative exact weights of these uniform type averages.  Omitted
types have weight zero.

```text
hard q15:
  K(0,-1,0) 10037/1920; K(-2,0,1) 853/960;
  K(0,1,-2) 11297/1920.                         total K = 12

hard q14, N=F:
  K(0,0,-1) 9499/960; K(0,-1,0) 9/32;
  K(-2,0,1) 791/960.                            total K = 11

hard q14, N!=F:
  K(0,-1,0) 1513/320; K(-2,0,1) 137/160;
  K(0,1,-2) 1733/320.                           total K = 11

opposite q16:
  K(0,-1,0) 1201/240; K(-2,1,0) 479/240;        total K = 7
  T(2,1,0) 143/60; T(0,0,3) 217/60.             total T = 6

opposite q15:
  K(0,-1,0) 12493/2880; K(-2,1,0) 4787/2880;   total K = 6
  T(2,1,0) 1499/720; T(0,0,3) 2821/720.         total T = 6
```

Dividing each displayed type coordinate by `(30,15,420)` and summing with
these weights reproduces (5) exactly.  Thus all 32 averaged rows lie in
their exact real atom cones.  By semimetric polarity, every row semimetric
inequality has nonpositive margin there, and so does every fixed nonnegative
sum of them.

## 4. Consequence and scope

This is stronger than observing that one learned bank fails on a second
graph: it excludes every fixed summed semimetric bank on the stated
fractional ledger relaxation.  It does not exclude an adaptive oracle which
chooses a violated row metric after seeing an integral graph, and it does not
address the integer atom semigroup or synchronization of atoms across rows.
The obstruction is deliberately for the still-open hard-fixed case, not the
already closed opposite-fixed `j=0` localized-Mobius design.

Replay with

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
  tests/test_p31_global_semimetric_obstruction.py
```

The executable certificate is
`src/e1_gmin_m4_p31_global_semimetric_obstruction.py`.  Residual (ii), the
hard-fixed integral fibre, and the `j=1` branches remain open.
