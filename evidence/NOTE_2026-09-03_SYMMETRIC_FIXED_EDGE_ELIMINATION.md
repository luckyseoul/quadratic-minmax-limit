# Parity-forced fixed-edge elimination in the symmetric Boolean fibre

Date: 2026-09-03

Status: proved exact all-prime reduction. The fixed antipodal edge variables
are uniquely forced by the fixed-cell parity of the central target, and the
remaining box is an equivalent divided zero-one system on unused nonfixed
inversion orbits. The reduced system is not solved, no common simple graph
is constructed, and residual (ii) remains open.

## 1. The restricted central fibre

Let \(p=2h+1\) be an odd prime, let \(d=p+1\), and let \(J(e)=-e\) be
central inversion on source edges. Let
\(\tau_e\in\{+1,-1\}\) be any source sign with
\(\tau_{Je}=\tau_e\). In the branch-C application this is the Paley column
sign.

Fix a ternary antisymmetric lift and let \(U\) be its **actual support** in
the set of nonfixed \(J\)-orbits: \(O\in U\) exactly when the final ternary
coefficient has absolute value one. Thus cancellations between overlapping
Möbius trades are removed before \(U\) is formed; \(U\) is not the union or
multiset of the individual trade supports. Select the corresponding physical
edge on each orbit and call its signed source \(q_U=\tau x_U\). The exact
central remainder is

\[
                         T_U=Y-Rq_U.                       \tag{1}
\]

Let \(F\) be the fixed antipodal source edges. On every unused nonfixed
orbit \(O=\{e,Je\}\), choose a representative and put

\[
 A_f=R(\tau_f f),\qquad
 B_O=R\bigl(\tau_O(e+Je)\bigr).                            \tag{2}
\]

The restricted symmetric fibre is exactly

\[
 T_U=\sum_{f\in F}a_fA_f+
       \sum_{O\notin U}b_OB_O,
 \qquad a_f,b_O\in\{0,1\}.                                \tag{3}
\]

The used orbits are frozen: changing their pair total would change the
already prescribed antisymmetric coefficient.

## 2. Exact block form

Split the central target coordinates into:

* fixed cells \(P_L\) and \(K_L(0,\beta)=\{s,-s\}\); and
* one common coordinate for every nonfixed pair
  \(K_L(\alpha,\beta),K_L(-\alpha,\beta)\).

This convention records the common value on a nonfixed pair, not the sum of
the two values. In these bases the symmetric map has the integral block form

\[
                 R^+=\begin{pmatrix}A&2B\\0&C\end{pmatrix}. \tag{4}
\]

Indeed, a fixed antipodal edge projects in every row either to \(P_L\) or to
a cell \(\{s,-s\}\), so the lower-left block is zero. If
\(O=\{e,Je\}\) is nonfixed, its two edges hit the two members of a nonfixed
target pair once each, but they hit any fixed target cell together. Since
their \(\tau\)-signs agree, the latter coefficient is even. This proves both
the zero and the factor two in (4) coefficientwise.

## 3. The fixed binary vector is unique and explicit

There are

\[
 |F|=dh,\qquad d(h+1)
\]

fixed source and raw fixed-target coordinates. The latter obey the \(d\)
independent fixed-coordinate compatibility equations from the symmetric
lattice theorem, so their compatible dimension is
\(d(h+1)-d=dh\). That theorem proves that the fixed-edge map modulo two has
rank \(dh\). It is therefore an isomorphism onto the compatible fixed-cell
residue space.

There is also a direct inverse. Index a fixed edge by
\([v]\in(V\setminus\{0\})/\{\pm1\}\). Write \(g_L(0)\) for the \(P_L\)
bit and, for a nonzero square \(\beta\), write \(g_L(\beta)\) for the bit of
the fixed cell \(\{s,-s\}\) with \(s^2=\beta\). Let \(L_v\) be the unique
projective functional annihilating \(v\). Then

\[
 \boxed{
 a_{[v]}=g_{L_v}(0)+\sum_{L\in\mathbf P(V^*)}
                   g_L\bigl(L(v)^2\bigr)\pmod2.}           \tag{5}
\]

To check (5), substitute the fixed-edge image of one class \([u]\). If
\([u]=[v]\), the sum sees that class in all \(p+1\) rows, an even number,
and the extra \(L_v\) term contributes one. If \(u\) is collinear with
\(v\) but \(u\ne\pm v\), only \(L_v\) occurs in the sum and it cancels the
extra term. If \(u,v\) are noncollinear, equality
\(L(u)^2=L(v)^2\) means that \(L\) annihilates \(u-v\) or \(u+v\); these
give two distinct rows, hence zero modulo two. Thus (5) is a left inverse,
and the equal dimensions make it the inverse.

For an integral compatible central \(T_U\), the nonfixed columns in (4)
vanish on fixed cells modulo two. Consequently (3) forces the unique vector

\[
       a(T_U)=\bar A^{-1}\bigl((T_U)_{\rm fix}\bmod2\bigr). \tag{6}
\]

Because \(a_f\) is binary, its residue determines the actual integer choice;
there is no remaining fixed-edge freedom.

## 4. Exact divided zero-one system

Subtract the forced fixed chain. Its fixed-cell remainder is even by (6).
Define

\[
 \widehat T_U=\left(
 { (T_U)_{\rm fix}-A a(T_U)\over2},
 (T_U)_{\rm pair}
 \right),\qquad
 \widehat B_O=\left(
 { (B_O)_{\rm fix}\over2},
 (B_O)_{\rm pair}
 \right).                                                 \tag{7}
\]

The pair coordinate is not divided, because it is the common value on the
two cells. Comparing the two blocks in (4) proves the exact equivalence

\[
 \boxed{
 (3)\ \Longleftrightarrow\
 \exists b\in\{0,1\}^{\Omega\setminus U}:
       \sum_{O\notin U}b_O\widehat B_O=\widehat T_U.}      \tag{8}
\]

Thus the \(dh\) fixed binary variables are eliminated, not relaxed. The
remaining problem is a genuine zero-one fibre on unused nonfixed orbits.

## 5. Exact Hamming and parallel-direction slices

A used orbit contributes one physical edge, a selected fixed coordinate one,
and a selected unused double orbit two. Hence every completion of prescribed
size \(|H|\) satisfies

\[
 \boxed{2\sum_{O\notin U}b_O=|H|-|U|-|a(T_U)|.}           \tag{9}
\]

For the Paley branch \(p\equiv3\pmod4\), the fixed parallel coordinates give
a strictly sharper product of slices. For a projective row \(L\), let
\(u_L\) be the number of orbits in \(U\) whose difference is parallel to
\(\ker L\), let \(f_L\) be the number of forced fixed edges parallel to it,
and let \(n_L\) be the number of selected unused double orbits parallel to
it. In the normalized Paley row, an edge parallel to \(L\) has
\(\tau_e=\epsilon_L\), so \(\epsilon_L\tau_e=1\). Therefore the physical
target quota \(P_L\) obeys

\[
 \boxed{P_L=u_L+f_L+2n_L,\qquad
        n_L={P_L-u_L-f_L\over2}.}                          \tag{10}
\]

There are \(h\) difference classes in \(\ker L\), and each has \(dh\)
nonfixed midpoint orbits. The unused capacity in this row is consequently
\(dh^2-u_L\). For the forced vector (6), parity of
\(P_L-u_L-f_L\) is automatic: it is exactly the \(P_L\) component of the
fixed-cell congruence. The first scalar obstruction exposed by (10) is
therefore nonnegativity; the complete necessary bounds are

\[
              0\le {P_L-u_L-f_L\over2}\le dh^2-u_L.       \tag{11}
\]

The upper capacity is recorded for exactness and may be vacuous on the
branch-C ranges. More importantly, the variables parallel to \(L\) lie on
the constant-weight slice of weight \(n_L\). Summing (10) over the \(d\)
directions recovers (9), so the global Hamming equation is the aggregate of
these directionwise slices.

## 6. The geometric fixed word of one used orbit

Write a nonfixed source orbit in midpoint/difference form

\[
                  e=\{a-\delta,a+\delta\},\qquad a\ne0.   \tag{12}
\]

Apply the inverse (5) to the fixed-cell incidence of one selected edge from
this orbit. Let \(\Phi(a,[\delta])\) be the resulting fixed-edge word. The
parallel row \(L_\delta\) contributes its \(P\)-bit once through the extra
term in (5) and once through the sum, so those two occurrences cancel. If
\(a\) and \(\delta\) are independent, the only other fixed cell occurs in
the row \(L_a\). A fixed class \([v]\) is toggled exactly when

\[
 L_a(v)^2=L_a(\delta)^2,
\]

or equivalently when \(v=\pm\delta+c a\). If \(a\) and \(\delta\) are
parallel, there is no second fixed row. Hence

\[
 \boxed{
 \Phi(a,[\delta])=
 \begin{cases}
 0,&a\parallel\delta,\\[2mm]
 \mathbf1_{\{[\delta+c a]:c\in\mathbf F_p\}},
     &a\not\parallel\delta.
 \end{cases}}                                             \tag{13}
\]

In the second case the \(p\) displayed antipodal classes are distinct.
Replacing \(a\) or \(\delta\) by its negative does not change the set.
Therefore \(\Phi\) depends only on the nonfixed inversion orbit and not on
which physical side of it the antisymmetric lift selects.

Let \(a_Y\) be formula (5) applied to the fixed-cell residue of the full
target. Linearity gives the exact support-only formula

\[
        a(T_U)=a_Y+\sum_{O\in U}\Phi(O)\pmod2.             \tag{14}
\]

Thus the first physical-size optimization attached to the antisymmetric
fibre is the coset weight

\[
 \boxed{
 \min_{\substack{Dw=(1-I)Y\\w\in\{-1,0,1\}^{\Omega}}}
 \left(
 |\operatorname{supp}w|+
 \left|a_Y+\sum_{O\in\operatorname{supp}w}\Phi(O)\right|
 \right).}                                                \tag{15}
\]

A value greater than \(|H|\) excludes a common graph. A value at most
\(|H|\) is only necessary, because the unused double-orbit transverse
equations in (8) remain.

### The paired-affine-line block basis

The nonzero sets in (13) are exactly the paired non-origin affine-line
blocks

\[
 {\cal B}(A,\bar c)=\bigl((c+A)\cup(-c+A)\bigr)/\{\pm1\},
 \qquad c\notin A.                                        \tag{16}
\]

There are \(d\) choices of the spatial direction \(A\) and \(h\) pairs of
nonzero affine cosets for each, hence \(dh=|\Delta|\) block types. Each block
contains \(p\) antipodal point classes. Let \(M\) be their square incidence
matrix over \(\mathbf F_2\), with rows indexed by \(\Delta\) and columns by
the blocks.

One point class lies on one non-origin affine-line pair for every direction
other than its own, hence on \(p\) blocks. Two distinct noncollinear classes
\([v],[w]\) lie together in exactly the two blocks with directions
\(\langle v-w\rangle\) and \(\langle v+w\rangle\). Distinct collinear
classes lie together in no non-origin block. Consequently

\[
                        MM^{\mathsf T}=I,
 \qquad M^{\mathsf T}M=I\quad\text{over }\mathbf F_2.      \tag{17}
\]

The nonzero \(\Phi\)-types therefore form a basis of the fixed-edge parity
space, not merely a spanning family. Among the \(|\Delta|^2\) nonfixed
source orbits, exactly \(|\Delta|h\) have \(a\parallel\delta\) and map to
zero. Every nonzero block type occurs with multiplicity \(ph\).

There is a useful lift of each basis block through the lower paired-cell
map \(C\). Fix one of the \(h\) midpoint magnitude classes \([a]\) in the
direction \(A\), and range \([\delta]\) over the \(p\) members of
\({\cal B}(A,\bar c)\). These \(p\) source columns all have fixed word
\({\cal B}(A,\bar c)\), so their odd sum retains that word. In a row not
annihilating \(a\), the values \(L(\delta)\) run once through
\(\mathbf F_p\): the zero value is fixed, while each nonzero squared value
occurs twice in the paired block and cancels modulo two. In the row
annihilating \(a\), every column is fixed. The \(p\)-column sum is therefore
in \(\ker C\). Varying \([a]\) gives \(h\) pairwise-disjoint
\(p\)-column \(C\)-kernel lifts of every block type.

Finally, one localized Möbius half cannot concentrate heavily in a midpoint
direction. In its normalized \((L,M)\)-coordinates, the midpoint of
parameter \(t\) has slope

\[
 {M(a_t)\over L(a_t)}
 ={t(t+2)\over(t+1)^2}
 =1-{1\over(t+1)^2}.                                     \tag{18}
\]

As \(t+1\) ranges over \(\mathbf F_p^*\), every occurring square has two
roots. The half occupies exactly \(h\) midpoint directions, twice in each.
In particular it deletes at most two columns from the collection with any
fixed midpoint direction. This controls how a sum of Möbius halves can
puncture the disjoint lifts above, but by itself does not prove the required
weight bound in (15).

### The cancellation-aware \(p=31\) parity ladder

There are two different uses of the word *parallel* in this reduction. In
this subsection a **zero-\(\Phi\) parallel orbit** means precisely
\(a\parallel\delta\) in (13), equivalently an orbit whose midpoint is
parallel to its half-difference. This is not the directionwise
parallel-edge statistic \(u_L\) in (10): every edge has one such direction,
whereas only the zero-\(\Phi\) orbits are counted here.

First sum the inverse (5) over all fixed source classes \([v]\). For a fixed
row \(L\), its \(P_L\) bit occurs twice for every \([v]\subset\ker L\), once
in the displayed extra term and once in the sum, and hence vanishes. Each
nonzero fixed-cell bit \(g_L(\beta)\) occurs for the \(p\) antipodal classes
with \(L(v)^2=\beta\), and \(p\) is odd. Therefore

\[
 |a(g)|\equiv\sum_L\sum_{\beta\ne0}g_L(\beta)\pmod2.      \tag{18a}
\]

For the central target of a graph with \(E\) edges, the nonfixed target
cells occur in equal opposite pairs and disappear from each row total
modulo two. Every row total is \(E\), while each source edge contributes to
exactly one \(P_L\), so \(\sum_Lg_L(0)\equiv E\). Since \(p+1\) is even,
(18a) gives

\[
                         |a_Y|\equiv E\pmod2.             \tag{18b}
\]

This uses the full target \(Y\), including any compact fixed-cell word; it
does not set that word to zero or determine its support.

At \(p=31\), the sixteen localized Möbius halves contribute \(480\) raw
orbit occurrences. In each half the parameter \(t=0\) is the unique
zero-\(\Phi\) orbit and the other \(29\) occurrences have odd, \(31\)-point
\(\Phi\) words. Hence the raw zero/nonzero-\(\Phi\) counts are \(16\) and
\(464\). At an orbit \(O\), let \(n_O\) be its raw multiplicity and let
\(c_O\in\{-1,0,1\}\) be its final coefficient in the ternary sum. Define
\(\kappa_O=(n_O-|c_O|)/2\). Splitting the cancellation units into
\(\kappa_0\) on zero-\(\Phi\) orbits and \(\kappa_1\) on nonzero-\(\Phi\)
orbits gives the exact post-cancellation counts

\[
 u_0=16-2\kappa_0,\qquad
 u_{\rm np}=464-2\kappa_1,\qquad
 |U|=480-2(\kappa_0+\kappa_1).                            \tag{18c}
\]

Thus both \(u_0\) and \(u_{\rm np}\) are even for every ternary overlap
pattern, including triple and higher overlaps. Weight parity is linear on
binary words, so cancellations and intersections among their affine-line
supports do not change the conclusion. With \(\kappa=\kappa_0+\kappa_1\)
and the branch-C edge count \(E=125+2t\), equations (14) and (18b)--(18c)
give

\[
 |a(T_U)|\equiv E+u_{\rm np}\equiv1,\qquad
 E-|U|-|a(T_U)|\equiv u_0\equiv0\pmod2.                  \tag{18d}
\]

The Hamming numerator parity is therefore automatic for the sixteen-half
construction: this parity calculation cannot exclude it. What it does
force, conditionally on a completion, is an odd number of fixed antipodal
edges.

Across \(68\le t\le177\), support size first becomes possible at
\(\kappa_{\min}=178-t\), and the remaining edge capacity is exactly

\[
 E-|U|=1+2\{\kappa-(178-t)\}.                             \tag{18e}
\]

At \(\kappa=\kappa_{\min}\), (18d)--(18e) and the exact Hamming equation
(9) force \(|a(T_U)|=1\) and zero selected unused double orbits in every
completion. For larger \(\kappa\), they give only the odd alternatives
\(1,3,\ldots,E-|U|\); the divided transverse equations (8) can still fail.
No compatible Möbius completion is constructed here, and residual (ii)
remains open.

## 7. The next legitimate linear gate

The symmetric-lattice note proved mod-two surjectivity and separately stated
the restricted box. The present theorem combines them: it gives the explicit
inverse, removes every fixed variable, divides the even block, and pins one
constant-weight slice per parallel direction. It retires only a repetition
of the first-layer or unrestricted parity and Smith calculations.

It does **not** retire parity after division and puncturing. The first
legitimate linear test on (8) is the image of the punctured halved joint code

\[
                  \begin{bmatrix}B\bmod2\\ C\bmod2\end{bmatrix},             \tag{19}
\]

equivalently the simultaneous fixed-cell congruences modulo four and paired
cell congruences modulo two after the forced \(a(T_U)\) and used columns have
been removed. Any later exchange or rounding theorem must also preserve all
\(p+1\) slice weights in (10), not merely their global sum.

The all-active theorem may enter through information about the actual support
\(U\), such as \(|U|\ge p-1\); another scalar sharpening of that bound does
not solve (8). What remains is to decide whether (8), on the directionwise
slices (10), is empty or contains a common zero-one lift for the branch-C
target. This note proves neither alternative. Residual (ii), E1,
\(L=1/2\), and the original MathOverflow limit remain OPEN.

## 8. Reproduction

The implementation records the symbolic proof. Its exact matrices at
\(p=3,7\) only check the formulas and are not theorem evidence.

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_symmetric_fixed_edge_elimination.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_fixed_edge_elimination.py
