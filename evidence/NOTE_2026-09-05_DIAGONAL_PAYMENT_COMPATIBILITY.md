# Conjugate payment, exact alternating cycles, and coherent-cross obstruction

Status: proved all-orders conjugate-pair identity, exact real-function
relaxation criterion, and coherent-cross no-improvement theorem (2026-09-05).
The proofs are analytic; finite regression checks are supplementary.
The statewise multiplier-two gate (33) remains OPEN.

This gives an exact characterization for **arbitrary real diagonal-payment
functions with the necessary phase symmetries**. It does not assert that
the resulting functions come from quadratic skew sign matrices. The latter
restriction, and the choice of a suitable cross block for an actually
globally minimal real signing, remain separate unresolved requirements.

## 1. Normalization and finite phase orbits

Let ell,r>=1, let A_L,A_R be real symmetric zero-diagonal matrices, and let
C,R be real ell-by-r matrices. In the signing application their appropriate
entries are signs. Write

\[
 H=\begin{pmatrix}A_L+iS_L&C+iR\\C^T-iR^T&A_R+iS_R\end{pmatrix},
 \qquad S_L^T=-S_L,\quad S_R^T=-S_R.
 \tag{1}
\]

Initially the skew entries may be arbitrary real numbers; actual completion
requires every off-diagonal skew entry to be a sign. Let

\[
 \rho(\xi)=\max(|\operatorname{Re}\xi|,|\operatorname{Im}\xi|),
 \qquad \nu_4(H)=\max_{x\in\mu_4^{\ell+r}}|x^*Hx|,
 \quad \mu_4=\{1,i,-1,-i\}.
\]

Use the finite orbit sets

\[
 V_L=\mu_4^\ell/C_4,\qquad V_R=\mu_4^r/C_4,
 \tag{2}
\]

where each block is independently identified under multiplication by a
common fourth root of unity. The involution kappa on either set is
coordinatewise conjugation. It is well-defined on these orbits.
An orbit is fixed by kappa exactly when it has a real Boolean
representative: normalize its first coordinate to one, after which a fixed
orbit consists of a vector equal to its conjugate.

For z in V_L and w in V_R, define

\[
 \begin{aligned}
 a(z,w)&=z^*A_Lz+w^*A_Rw,\\
 u(z)&=z^*iS_Lz,\qquad v(w)=w^*iS_Rw,\\
 g_\pm(z,w)&=\rho\bigl(z^*(C\pm iR)w\bigr).
 \end{aligned}
 \tag{3}
\]

These are well-defined on the independent phase orbits. The quantities
a,u,v are real. Conjugation obeys

\[
 a(\kappa z,\kappa w)=a(z,w),\quad
 u(\kappa z)=-u(z),\quad v(\kappa w)=-v(w),\quad
 g_+(\kappa z,\kappa w)=g_-(z,w).
 \tag{4}
\]

For the last identity, the conjugated cross scalar is the complex
conjugate of z^*(C-iR)w, and rho is invariant under conjugation.
The same orbit normalization handles the rotated fourth-phase cube used
elsewhere in the repository, since a common rotation leaves Hermitian
quadratics and cross bilinear quantities unchanged.

Rotating w through four phases does not change its diagonal energy and
lets the cross term add with the favorable sign and Cartesian component.
Therefore the exact norm is

\[
 \nu_4(H)=\max_{z,w}\bigl(|a(z,w)+u(z)+v(w)|+2g_+(z,w)\bigr).
 \tag{5}
\]

## 2. Exact conjugate-pair identity

For one orbit pair suppress the arguments and set

\[
 m=g_++g_-,\qquad h=g_+-g_-,
\]

\[
 c=m+\max(|a|,|h|),\qquad
 d_0={|h-a|-|h+a|\over2}.
 \tag{6}
\]

Then, for **every real d**, not merely at its minimizer,

\[
 \boxed{
 \max\{|a+d|+2g_+,\ |a-d|+2g_-\}
       =c+|d-d_0|.}
 \tag{7}
\]

Proof: expand the two absolute values into four affine functions. Their
maximum is

\[
 \max\{d+m+|a+h|,\ -d+m+|a-h|\}.
\]

For real A,B, max(d+A,-d+B)=(A+B)/2+|d-(B-A)/2|.
Now use (|a+h|+|a-h|)/2=max(|a|,|h|), proving (7).

In particular,

\[
 \min_d\max\{|a+d|+2g_+,|a-d|+2g_-\}
 =c=\max\{2\max(g_+,g_-),\ |a|+g_++g_-\}.
 \tag{8}
\]

The optimal transfer d_0 is independent of the target T. It can also be
written -sign(a) times h clipped to [-|a|,|a|], with value zero when a=0.
The permitted transfer interval at target T is exactly

\[
 I_{zw}(T)=[L_{zw}(T),U_{zw}(T)]
 =[d_0-(T-c),\ d_0+(T-c)].
 \tag{9}
\]

If T<c the interval is empty, as it should be. Equivalently,

\[
 \begin{aligned}
 L&=\max(-a-T+2g_+,\ a-T+2g_-),\\
 U&=\min(-a+T-2g_+,\ a+T-2g_-).
 \end{aligned}
\]

Simultaneous conjugation fixes c and reverses d_0:

\[
 c_{\kappa z,\kappa w}=c_{zw},\qquad
 (d_0)_{\kappa z,\kappa w}=-(d_0)_{zw},\qquad
 I_{\kappa z,\kappa w}(T)=-I_{zw}(T).
 \tag{10}
\]

Because (5) ranges over every conjugate pair, physical skew blocks satisfy
the further **exact identity**

\[
 \boxed{\nu_4(H)=\max_{z,w}
       \{c_{zw}+|u(z)+v(w)-(d_0)_{zw}|\}.}
 \tag{11}
\]

Thus the remaining diagonal task is a weighted approximation of the joint
transfer d_0 by the sum of one left and one right payment function. The
symmetry and quadratic/sign realization of those functions must not be
discarded when interpreting this identity.

## 3. Exact criterion for arbitrary real payment functions

Temporarily permit arbitrary functions u:V_L->R and v:V_R->R, subject only
to u(kappa z)=-u(z) and v(kappa w)=-v(w). Define their relaxed norm by the
right side of (11). There exist such functions with relaxed norm <=T if
and only if the cycle inequalities below hold.

Introduce vertex potentials

\[
 t_L(z)=u(z),\qquad t_R(w)=-v(w).
\]

The interval conditions become

\[
 L_{zw}\le t_L(z)-t_R(w)\le U_{zw}.
 \tag{12}
\]

On the complete bipartite graph put a directed arc w->z of weight U_zw
and an arc z->w of weight -L_zw. These are exactly difference constraints.
They have a real solution if and only if every directed cycle has
nonnegative total weight.

Necessity follows by adding the inequalities around a cycle. For
sufficiency, add a source with zero-weight arcs to every vertex. If there
is no negative cycle, shortest-path distances are finite (remove
nonnegative cycles from any path); they satisfy every arc inequality and
therefore supply the required vertex potentials. This is a finite
existence argument, not a proposed phase enumeration.

Every simple directed cycle is either a two-cycle on one bipartite edge
or an alternating cycle of length 2k>=4. With z_1,...,z_k and w_1,...,w_k
and cyclic index w_(k+1)=w_1, its inequality is

\[
 \boxed{\sum_{i=1}^k L_{z_i,w_{i+1}}
                \le\sum_{i=1}^k U_{z_i,w_i}.}
 \tag{13}
\]

The two-cycle is just L_zw<=U_zw, equivalently T>=c_zw. Simple cycles
suffice because every directed closed walk decomposes into them.

### Conjugacy is preserved without an extra obstruction

Suppose arbitrary potentials t solve (12). Define

\[
 t'(q)=-t(\kappa q)
\]

on both parts. By (10), t' also solves (12): its left-minus-right
difference is the negative of a feasible difference at the conjugate
pair. Convexity then makes

\[
 \widetilde t(q)={t(q)-t(\kappa q)\over2}
 \tag{14}
\]

feasible and odd under kappa. Set u=tilde t_L and v=-tilde t_R. These have
the required separate conjugacy oddness and vanish on all fixed real
phase orbits. Independent global fourth-phase invariance was already
built into the vertex sets (2). Thus imposing these symmetries does not
change the difference-constraint feasibility criterion.

### Exact all-cycle threshold

For an undirected simple alternating cycle E of length 2k>=4, let

\[
 C_E=\sum_{e\in E}c_e,\qquad
 D_E=\sum_{i=1}^k\bigl((d_0)_{z_i,w_i}-(d_0)_{z_i,w_{i+1}}\bigr).
\]

Its two orientations in (13) give precisely

\[
 2kT\ge C_E+|D_E|.
\]

Include a repeated single edge as a degenerate two-cycle, with C_E=2c_e
and D_E=0. The exact relaxed optimum, attained by real odd functions, is

\[
 \boxed{T_{\rm rel}(C,R;A_L,A_R)
       =\max_E {C_E+|D_E|\over |E|}.}
 \tag{15}
\]

In particular each four-cycle gives the necessary lower bound

\[
 T\ge {c_{11}+c_{12}+c_{21}+c_{22}
       +|(d_0)_{11}+(d_0)_{22}-(d_0)_{12}-(d_0)_{21}|\over4}.
 \tag{16}
\]

The alternating transfer circulation is an explicit separability penalty,
in addition to the individual conjugate-pair costs. It is not a Gram-only
quantity: a includes both actual real diagonal blocks.

## 4. Four-cycles alone are not a general interval criterion

Consider the 3-by-3 interval system

\[
 I_{11}=I_{22}=I_{33}=\{0\},\quad
 I_{12}=I_{23}=I_{31}=[1,2],\quad
 I_{21}=I_{32}=I_{13}=[-2,2].
 \tag{17}
\]

Every interval is nonempty. Every two-by-two rectangle passes both
four-cycle inequalities:

* If its row and column index sets agree, it has two diagonal entries.
  The sum of their lower endpoints is 0, at most the off-diagonal upper
  sum 4. The off-diagonal lower sum is 1+(-2)=-1, at most the diagonal
  upper sum 0.
* Otherwise it has exactly one diagonal entry. In one inequality the
  lower sum is diagonal zero plus an off-diagonal lower endpoint, hence
  at most 1, while the other upper sum is 4. In the reverse inequality
  the lower sum is at most 1+1=2, while the upper sum is 0+2=2.

Nevertheless the six-cycle with lower edges (1,2),(2,3),(3,1) and upper
edges (1,1),(2,2),(3,3) requires 3<=0. Equivalently, diagonal constraints
u_i+v_i=0 sum to zero, while the three cyclic lower constraints demand
the same sum be at least three.

This is an interval-system counterexample to a four-cycle-only theorem.
It is **not** asserted to arise from physical sign matrices C,R or from a
globally minimal real signing. The complete criterion (15) does not need
such an example or any assumption that physical intervals are arbitrary.

## 5. Quadratic and sign realization remain additional constraints

For an actual real skew matrix S_L, direct expansion gives

\[
 u(z)=-2\sum_{i<j}(S_L)_{ij}\operatorname{Im}(\overline z_i z_j).
 \tag{18}
\]

Thus even arbitrary **real-skew** realizability confines u to the
binomial(ell,2)-dimensional span of these imaginary pair characters.
For the Fourier convention chi_ij(z)=bar z_i z_j, (18) has coefficient
+i(S_L)_ij at chi_ij and -i(S_L)_ij at its inverse, and no other
Fourier characters. Actual skew signing requires (S_L)_ij in {+1,-1},
so these Fourier coefficients must be exactly +i or -i with the stated
inverse relation. The same applies on the right.

The orbit symmetries alone do not imply even real quadratic realization.
For ell>=3 the function

\[
 z\longmapsto\operatorname{Im}(\overline z_1^{\,2}z_2z_3)
 \tag{19}
\]

is invariant under any common phase, is odd under conjugation, and vanishes
on the fixed real phase orbits. On the fourth-phase cube its Fourier
characters are (2,1,1,0,...) and (2,3,3,0,...), neither of which is a
pair-difference character. Character orthogonality therefore proves that
it is outside the space (18).

There is an additional elementary value-lattice restriction in the sign
case. If a phase vector has p coordinates in one Cartesian axis class
and ell-p in the other, exactly p(ell-p) terms in (18) are nonzero.
Consequently u(z)/2 is an integer congruent to p(ell-p) modulo two.
This is necessary, not sufficient, for simultaneous sign realization.

Therefore the cycle theorem is exact only for the stated real-function
relaxation. It supplies necessary conditions for actual skew blocks, not
a rounding theorem into them.

## 6. Uniform obstruction for coherent cross choices R=+C or R=-C

Let xi=z^*Cw=s+it. For R=+C or R=-C,

\[
 g_+=g_-=\rho((1+i)\xi)=\rho((1-i)\xi)=|s|+|t|.
 \tag{20}
\]

Hence h=0, d_0=0, and c=|a|+2(|s|+|t|). Formula (11) becomes

\[
 \boxed{\nu_4(H)=\max_{z,w}
       \{|a(z,w)|+|u(z)+v(w)|
                  +2(|\operatorname{Re}\xi|+|\operatorname{Im}\xi|)\}.}
 \tag{21}
\]

In particular, setting S_L=S_R=0 minimizes this norm over **all real skew
diagonal choices**, and even over all arbitrary real odd payment functions:

\[
 \nu_4(H)\ge
 \nu_4\begin{pmatrix}A_L&(1\pm i)C\\(1\mp i)C^T&A_R\end{pmatrix}.
 \tag{22}
\]

The zero choice is only a relaxation when actual diagonal entries must be
skew signs. The rigorous conclusion is that no diagonal skew filling can
repair a coherent cross block by lowering its norm. This does not prove
that every coherent cross choice exceeds the target for every real A;
it rules out diagonal compensation as a rescue when that cross choice
already fails.

## 7. Relation to the still-open multiplier-two gate

If the real signing A is globally minimal with M=Phi(A), equation (33) in
the cross-rectangle note asks for an actual skew completion satisfying

\[
 \nu_4(H)\le 2\sqrt2 M+o_{\rm Dini}(n^{3/2}).
\]

For each proposed real cross signing R, formula (15) is an exact, generally
stronger necessary test than checking each conjugate pair separately. If
T_rel exceeds the target, no choice of the diagonal skew blocks can work.
If T_rel is below the target, the cycle theorem constructs only real
payment functions, not quadratic skew signs. Actual sign realization
requires the finite Fourier constraints (18) and flat coefficient choices.

Neither a suitable noncoherent R nor a theorem using global minimality to
control these cycle costs and realize the payments has been proved here.
There is no multiplier-two, MathOverflow-limit, or global closure claim.

## Exact reference implementation

`src/original_mo_diagonal_compatibility.py` computes (6) with rational
arithmetic and certifies generic interval feasibility by difference
constraints. It returns either additive potentials or a directed negative
cycle whose input inequalities sum to a strict contradiction. It does not
validate physical phase tables, impose an involution on the returned
potentials, or assert skew realization; (14) is the analytic symmetrization
step when the physical symmetry hypotheses hold.
`tests/test_original_mo_diagonal_compatibility.py` independently replays
the certificates and checks actual coherent Hermitian examples. These
finite checks do not replace the all-orders proofs above.

The complementary actual-energy covariance identities are proved in
`NOTE_2026-09-05_ACTUAL_DIAGONAL_MIXED_MOMENTS.md`. They are averaged
identities, not pointwise cycle or realization bounds.

## Appendix: a physical order-four pair, cycle, and sign separation

This is an analytic example using the definitions above, not an orientation
or signing census. Its exact phase table also has an independent direct-
energy regression check in `tests/test_original_mo_diagonal_compatibility.py`.

The example proves that even global minimality of the real signing does
not make the additive-payment cycle threshold equal the largest individual
pair threshold. It also exhibits a further strict skew-sign realization
gap. It does not refute the asymptotic multiplier-two gate.

### A.1. Physical blocks and the real global minimum

Take two coordinates on each side and set

\[
 A_L=A_R=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
 C=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad
 R=\begin{pmatrix}-1&1\\1&1\end{pmatrix}.
 \tag{A1}
\]

These are genuine sign blocks. With
\(J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}\), they satisfy
\(R=JC=-CJ\). The real symmetric signing

\[
 A=\begin{pmatrix}A_L&C\\C^T&A_R\end{pmatrix}
\]

has exactly one negative edge. Under
\(Q_A(x)=\sum_{i<j}A_{ij}x_ix_j\) and
\(\Phi(A)=\max_{x\in\{\pm1\}^4}|Q_A(x)|\), its norm is four.
Indeed, \(Q_A(\mathbf1)=4\); no switched copy has all six edges positive
or all six edges negative, because A has both a positive and a negative
triangle product. Thus neither score six nor score minus six occurs, and
every score is even. Moreover every order-four signing B satisfies
\(\mathbb E Q_B(x)^2=6\), while every score is even; if \(\Phi(B)<4\),
all scores would have square at most four. Therefore A is an actual global
minimizer, with \(M=4\).

### A.2. Exact type table

Normalize the independent phase orbits as
\(z=(1,\zeta)\), \(w=(1,\omega)\), where
\(\zeta,\omega\in\mu_4\). Write \(\sigma,\tau\in\{\pm1\}\).
Direct multiplication of the matrices in (A1) gives the following table.
The letters R and I designate a real or imaginary second coordinate.

| Type | \(\zeta\) | \(\omega\) | \(a\) | \(g_+\) | \(g_-\) | \(c\) | \(d_0\) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RR | \(\sigma\) | \(\tau\) | \(2(\sigma+\tau)\) | 2 | 2 | \(4+|2(\sigma+\tau)|\) | 0 |
| RI | \(\sigma\) | \(i\tau\) | \(2\sigma\) | \(2(1-\tau)\) | \(2(1+\tau)\) | 8 | \(2\sigma\tau\) |
| IR | \(i\sigma\) | \(\tau\) | \(2\tau\) | \(2(1+\sigma)\) | \(2(1-\sigma)\) | 8 | \(-2\sigma\tau\) |
| II, equal | \(i\sigma\) | \(i\sigma\) | 0 | 0 | 0 | 0 | 0 |
| II, opposite | \(i\sigma\) | \(-i\sigma\) | 0 | \(2(1+\sigma)\) | \(2(1-\sigma)\) | 8 | 0 |

For example, in type IR put \(\xi=z^*Cw\). The identity
\(z^*J=-i\sigma z^*\) yields
\(z^*Rw=-i\sigma\xi\), while \(\rho(\xi)=2\). Hence
\(g_+=2(1+\sigma)\) and \(g_-=2(1-\sigma)\).
For type RI, use \(Jw=-i\tau w\) and \(R=-CJ\) instead.
In type II the cross scalars vanish when \(\sigma=\tau\); when
\(\sigma=-\tau\), the same IR identity gives the last row.
In particular,

\[
 \max_{z,w}c_{zw}=8.
 \tag{A2}
\]

### A.3. Exact relaxed threshold and its four-cycle obstruction

Every conjugacy-odd payment function on either two-coordinate phase orbit
set has the form

\[
 u(\pm1)=v(\pm1)=0,\qquad
 u(i\sigma)=\sigma U,\qquad v(i\tau)=\tau V
 \quad(U,V\in\mathbb R).
 \tag{A3}
\]

Using the exact pair identity in each row of the table gives

\[
 \boxed{N(U,V)=
 \max\{8,\ 10+|U|,\ 10+|V|,\ |U+V|,\ 8+|U-V|\}.}
 \tag{A4}
\]

For example, the IR rows contribute
\(8+\max_{\tau=\pm1}|U+2\tau|=10+|U|\), and the RI rows
contribute \(10+|V|\). Thus

\[
 \boxed{T_{\rm rel}=\min_{U,V\in\mathbb R}N(U,V)=10>8=\max c,}
 \tag{A5}
\]

with the minimum attained at \(U=V=0\).

This is an actual physical four-cycle circulation obstruction, not an
abstract interval-system counterexample. On the rectangle with left
vertices \(\zeta=1,-1\) and right vertices \(\omega=i,-i\), all
four c values are eight and the d_0 table is

\[
 \begin{pmatrix}2&-2\\-2&2\end{pmatrix}.
\]

Its alternating circulation has magnitude eight, so the cycle bound is
\((4\cdot8+8)/4=10\). This also proves the lower bound without
imposing oddness on the payments.

### A.4. Strict skew-sign realization gap

Write
\(S_L=\begin{pmatrix}0&s_L\\-s_L&0\end{pmatrix}\) and
\(S_R=\begin{pmatrix}0&s_R\\-s_R&0\end{pmatrix}\).
Then \(U=-2s_L\) and \(V=-2s_R\). For actual skew-sign blocks,
\(s_L,s_R\in\{\pm1\}\), hence \(U,V\in\{\pm2\}\).
Formula (A4) now yields, for every one of these sign choices,

\[
 \boxed{\nu_4\!\begin{pmatrix}
 A_L+iS_L&C+iR\\C^T-iR^T&A_R+iS_R
 \end{pmatrix}=12>10=T_{\rm rel}>8=\max c.}
 \tag{A6}
\]

There is no additional real-quadratic realization gap in this small
example: any U,V in (A3) are realized by real skew entries
\(s_L=-U/2,s_R=-V/2\). The second strict gap in (A6) comes specifically
from requiring those entries to be signs.

### A.5. Scope and the universal elementary bound

For every physical table,

\[
 |d_0|=\min(|a|,|h|)\le\tfrac12c,
\]

because \(g_++g_-\ge|h|\) and
\(c=g_++g_-+\max(|a|,|h|)\).
The zero payments are admissible in the real-function relaxation, so

\[
 \max c\le T_{\rm rel}
 \le\max_{z,w}(c_{zw}+|(d_0)_{zw}|)
 \le\tfrac32\max c.
 \tag{A7}
\]

The example proves that the first inequality can be strict even when the
real A is globally minimal. It does not assert sharpness of the factor
three-halves. Here the relaxed value ten is below \(2\sqrt2 M=8\sqrt2\),
while every skew-sign diagonal completion for this fixed R has value
twelve above it. Other cross choices R are not excluded, and a fixed
order-four discrepancy does not contradict the asymptotic Dini allowance.
The original multiplier-two gate remains OPEN.
