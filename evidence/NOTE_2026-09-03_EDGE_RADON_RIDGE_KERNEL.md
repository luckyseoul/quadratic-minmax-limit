# An explicit \(p\)-saturated ridge sublattice of the edge--Radon kernel

**Status:** proved all-prime kernel theorem and constructive descent
inequalities; it does not prove signed-Boolean existence or close residual
(ii).

This note continues exactly from the common-graph gate (15.760.10). It does
not repeat a finite-prime census, parity/Smith obstruction, least-norm
calculation, or Parseval bound. Its purpose is to replace the abstract
rational nullspace by explicit integral moves and to state exactly what
those moves still miss.

## 1. Midpoint coordinates and two ridge families

Write every edge of \(V=\mathbf F_p^2\), uniquely, as

\[
 e=(a,[\delta])=\{a-\delta,a+\delta\},
 \qquad a\in V,
 \quad [\delta]\in\Delta:=(V\setminus\{0\})/\{\pm1\}.
                                                               \tag{1}
\]

For a projective nonzero functional \(L\), its edge--Radon row is

\[
 \rho_L(a,[\delta])=
 \begin{cases}
 P_L,&L(\delta)=0,\\
 K_L\bigl(L(a),L(\delta)^2\bigr),&L(\delta)\ne0.
 \end{cases}                                                   \tag{2}
\]

Let \(g:\mathbf F_p\to\mathbf Z\) satisfy
\(\sum_\alpha g(\alpha)=0\).

**Type P (parallel ridge).** Fix \([\delta]\) with \(L(\delta)=0\), and set

\[
 u^{P}_{L,\delta,g}(a,[\epsilon])=
 \begin{cases}
 g(L(a)),&[\epsilon]=[\delta],\\
 0,&\text{otherwise}.
 \end{cases}                                                   \tag{3}
\]

**Type K (same-square transverse ridge).** Fix distinct
\([\delta_1],[\delta_2]\) such that

\[
 L(\delta_1)^2=L(\delta_2)^2=\beta\ne0,
\]

and set

\[
 u^{K}_{L,\delta_1,\delta_2,g}(a,[\epsilon])=
 \begin{cases}
  g(L(a)),&[\epsilon]=[\delta_1],\\
 -g(L(a)),&[\epsilon]=[\delta_2],\\
 0,&\text{otherwise}.
 \end{cases}                                                   \tag{4}
\]

**Lemma 1.** Every move (3)--(4) belongs to
\(\ker_{\mathbf Z}R\).

*Proof.* For (3), the \(P_L\) row receives
\(p\sum_\alpha g(\alpha)=0\). In a different direction \(M\), the
restriction of \(L\) to every affine fibre \(M(a)=c\) is bijective, so each
\(K_M(c,M(\delta)^2)\) row receives \(\sum_\alpha g(\alpha)=0\).

For (4), the two difference classes occupy the same
\(K_L(\alpha,\beta)\) row and cancel there pointwise. In a direction
\(M\ne L\), each difference class cancels separately: if
\(M(\delta_i)=0\), its parallel sum is \(p\sum g=0\); otherwise \(L\) is
bijective on the affine fibre \(M(a)=c\), and its \(K_M\) row again receives
\(\sum g=0\). \(\square\)

For the elementary profile
\(g=\mathbf1_\alpha-\mathbf1_{\alpha'}\), a Type-P move has support and
squared norm \(2p\); a Type-K move has support and squared norm \(4p\).
These supports grow linearly with \(p\), consistently with the lower bound
\(p+1\) in evidence/NOTE_2026-09-02_EDGE_RADON_SWITCHING_DEGREE.md.

**Graver sublemma.** The elementary Type-P moves are primitive circuits,
and hence Graver elements, for every odd prime \(p\). The elementary Type-K
moves are primitive circuits, and hence Graver elements, for \(p\ge5\).

*Proof.* Apply an invertible affine change of coordinates and take
\(L(x,y)=x\). First let \(\delta\) be parallel to \(L\), and suppose that
\(v\in\ker_{\mathbf Z}R\) is supported inside one elementary Type-P support,
the two midpoint lines \(x=\alpha,\alpha'\) for this fixed difference
class. The other projective directions are
\(M_c(x,y)=y+cx\), \(c\in\mathbf F_p\). In every \(M_c\)-row, the unique
supported midpoint on \(x=\alpha\) is paired with the unique one on
\(x=\alpha'\) having the same \(M_c\)-value. Comparing the pairings for two
distinct \(c,c'\) shows that the coefficients on \(x=\alpha\) are invariant
under translation by
\((c-c')(\alpha-\alpha')\ne0\). A nonzero translation is transitive on an
affine \(\mathbf F_p\)-line, so those coefficients are constant, and the
paired coefficients on the other line are their negatives. Thus \(v\) is
a scalar multiple of the displayed Type-P move.

For Type K, change the signs of the two difference representatives so that
\[
 \delta_1=(r,s_1),\qquad \delta_2=(r,s_2),
 \qquad r\ne0,\quad s_1\ne s_2.
\]
In direction \(M_c\), a difference class is parallel precisely at
\(c_i=-s_i/r\), while the two projected squares agree precisely at
\[
 c_0=-{s_1+s_2\over2r}.
\]
The three values \(c_0,c_1,c_2\) are distinct. For \(p\ge5\), choose two
distinct \(c,c'\) outside them. In either direction, the two difference
classes occupy different nonparallel square rows. Therefore the preceding
pairing argument applies separately to each difference class and makes its
coefficients constant on each of the two \(L\)-lines, with opposite
constants between the lines. The \(K_L(\alpha,r^2)\) row then makes the two
difference-class constants opposite. Hence every kernel vector supported
inside the elementary Type-K support is a scalar multiple of that move.

Both displayed moves have coefficients in \(\{0,\pm1\}\), so they are
primitive. The support calculation proves they are circuits. If a
nonzero kernel vector were conformal to and coordinatewise bounded by one
of them, circuit minimality would make it a scalar multiple, and
primitivity would make it the whole move. This is exactly Graver
minimality. \(\square\)

No Type-K Graver assertion is made at \(p=3\), where there are not two
directions outside the three exceptional values. This exception is
irrelevant to the large-prime compact rays but keeps the theorem exact.

## 2. Exact \(p\)-saturation

Let \({\cal K}_{\rm ridge}\) be the integer lattice generated by (3)--(4).

**Theorem 2 (ridge \(p\)-saturation).** For every odd prime \(p\),

\[
 \boxed{
 p\ker_{\mathbf Z}R\ \subseteq\ {\cal K}_{\rm ridge}
 \ \subseteq\ \ker_{\mathbf Z}R.}                            \tag{5}
\]

Consequently

\[
 \ker_{\mathbf Q}R={\cal K}_{\rm ridge}\otimes\mathbf Q,
 \qquad
 \ker_{\mathbf Z}R/{\cal K}_{\rm ridge}
 \text{ is a finite elementary \(p\)-group}.                  \tag{6}
\]

*Proof.* The second inclusion is Lemma 1. Let
\(u\in\ker_{\mathbf Z}R\), and put

\[
 t_\delta=\sum_{a\in V}u(a,[\delta]).
\]

After forgetting midpoints, \(t=(t_\delta)_\delta\) lies in the kernel of
the pure-difference map \(S\) in (15.760.5). Proposition 15.760 proves that
\(S\) is injective over \(\mathbf Z\). Hence

\[
                       t_\delta=0\quad([\delta]\in\Delta).   \tag{7}
\]

For every direction \(L\), difference class \([\delta]\), and
\(\alpha\in\mathbf F_p\), define its fibre sum

\[
 G_{L,\delta}(\alpha)
   :=\sum_{b:L(b)=\alpha}u(b,[\delta]).                        \tag{8}
\]

Equation (7) says that every \(G_{L,\delta}\) has zero total. The elementary
affine-plane incidence identity

\[
 \begin{aligned}
 \sum_LG_{L,\delta}(L(a))
 &=\sum_{b\in V}u(b,[\delta])
       \#\{L:L(b)=L(a)\}\\
 &=(p+1)u(a,[\delta])+
       \sum_{b\ne a}u(b,[\delta])\\
 &=p\,u(a,[\delta])                                           \tag{9}
 \end{aligned}
\]

uses the fact that a nonzero vector \(b-a\) is annihilated by exactly one
projective functional.

Fix \(L\). If \(L(\delta)=0\), the summand
\(G_{L,\delta}(L(a))\) is exactly a Type-P move. If
\(L(\delta)^2=\beta\ne0\), the kernel equation in the row
\(K_L(\alpha,\beta)\) gives

\[
 \sum_{[\delta]:L(\delta)^2=\beta}
          G_{L,\delta}(\alpha)=0                              \tag{10}
\]

for every \(\alpha\). There are exactly \(p\) such difference classes.
Choose one as a base. Equation (10) writes this entire same-square summand
as the integral sum, over the other \(p-1\) classes, of Type-K moves with
profile \(G_{L,\delta}\). Thus the right side of (9), summed over all
\([\delta]\), lies in \({\cal K}_{\rm ridge}\). It equals \(pu\), proving
the first inclusion in (5).

Both lattices have the same rank by (5), so their quotient is finite; (5)
also kills every quotient element by \(p\). This proves (6). \(\square\)

This proof is constructive: equations (8)--(10) explicitly decompose
\(pu\) into ridge moves. It does **not** divide that decomposition by \(p\)
inside the integer lattice.

## 3. Rank audit

Put \(m=(p-1)/2\) and \(d=p+1\). There are \(dm\) difference classes and
\(p^2dm\) edge coordinates. Proposition 15.760 gives
\(\operatorname {rank}R=dpm\), so

\[
 \operatorname {rank}\ker R=dmp(p-1).                         \tag{11}
\]

Choose the \(p-1\) profiles
\(g_\alpha=\mathbf1_\alpha-\mathbf1_0\), each of the \(m\) parallel
difference classes in every direction, and a base among the \(p\) classes
in each of its \(m\) nonzero-square fibres. The resulting counts are

\[
 \begin{array}{c|c}
 \text{family}&\text{number}\\ \hline
 \text{Type P}&dm(p-1)\\
 \text{Type K}&dm(p-1)^2\\ \hline
 \text{total}&dmp(p-1).
 \end{array}                                                   \tag{12}
\]

They generate every profile used in the proof of Theorem 2; an arbitrary
Type-K pair is the difference of the two corresponding base-star moves.
Hence they generate \({\cal K}_{\rm ridge}\); (5), (11), and the matching
count in (12) show that they are a rational basis of the full kernel and an
integer basis of the ridge sublattice. They need not be a basis of the
saturated integer kernel.

### Exact one-step saturation invariant

Let \(B_p\) be the integer edge-by-ridge matrix whose columns are the
canonical basis in (12), so
\({\cal K}_{\rm ridge}=B_p\mathbf Z^{dmp(p-1)}\). Theorem 2 gives the
stronger exact identity

\[
 \ker_{\mathbf Z}R
 =\{v\in\mathbf Z^E:p v\in{\cal K}_{\rm ridge}\}.             \tag{12a}
\]

Indeed, one inclusion is (5). Conversely, if \(pv\in{\cal K}_{\rm ridge}\),
then \(pRv=0\); the integral target lattice is torsion-free, so \(Rv=0\).
Thus only one \(p\)-saturation step is required.

The remaining quotient has a completely explicit linear invariant:

\[
 \boxed{\quad
 \ker_{\mathbf Z}R/{\cal K}_{\rm ridge}
 \ \cong\ \ker\!\left(
   \overline {B_p}:\mathbf F_p^{dmp(p-1)}\to\mathbf F_p^E
 \right).
 \quad}                                                       \tag{12b}
\]

For \(\bar c\in\ker_{\mathbf F_p}\overline {B_p}\), choose an integral lift
\(c\). Then \(B_pc\) is coordinatewise divisible by \(p\), and the
isomorphism is

\[
             \bar c\longmapsto {B_pc\over p}
                       \pmod{{\cal K}_{\rm ridge}}.            \tag{12c}
\]

Changing the lift adds a ridge vector. If (12c) is zero, then
\(B_p(c-pq)=0\) for some integral \(q\); the columns of \(B_p\) are
rationally independent by (11)--(12), so \(\bar c=0\). Surjectivity follows
from (12a): for \(v\in\ker_{\mathbf Z}R\), write \(pv=B_pc\), and then
\(\bar c\) is a dependency of \(\overline {B_p}\) mapping to \(v\).

Consequently the exact new residual invariant is

\[
 \nu_p:=\dim_{\mathbf F_p}\ker\overline {B_p},\qquad
 [\ker_{\mathbf Z}R:{\cal K}_{\rm ridge}]=p^{\nu_p}.          \tag{12d}
\]

This invariant has a closed form. Let
\[
 {\cal U}:=\bigoplus_{L,[\delta]}
 \{g:\mathbf F_p\to\mathbf Z:\sum_\alpha g(\alpha)=0\},
\]
and define the unrestricted ridge-synthesis map
\[
 {\cal T}:{\cal U}\longrightarrow E_0,\qquad
 ({\cal T}g)(a,[\delta])=\sum_Lg_{L,\delta}(L(a)).
                                                               \tag{12e}
\]
Here \(E_0\) is the lattice of source arrays with zero midpoint total in
each difference class. Domain and target both have rank
\(dm(p^2-1)\). Equation (9) gives \(pE_0\subseteq\operatorname {im}{\cal T}\),
so \({\cal T}\) is injective and its cokernel is elementary \(p\)-torsion.

Modulo \(p\), a function on \(\mathbf F_p\) has sum zero exactly when its
interpolation polynomial has degree at most \(p-2\): the only coefficient
detected by summation is that of \(x^{p-1}\). For fixed degree
\(0\le i\le p-2\), the projective powers \(L(a)^i\) span the \(i+1\)
binary forms of degree \(i\). Different degrees are independent as
functions on \(\mathbf F_p^2\). Hence the reduction of \({\cal T}\), for
each fixed difference class, has rank
\[
                    \sum_{i=0}^{p-2}(i+1)=pm.
\]
It follows from the elementary-\(p\) cokernel that
\[
 [E_0:{\cal T}{\cal U}]=p^A,\qquad
 A=dm\bigl((p^2-1)-pm\bigr).                                  \tag{12f}
\]

Let \(M\) be Proposition 15.760's midpoint target lattice, of rank
\[
                         r_M=dm(p-1),
\]
and define \({\cal D}:{\cal U}\to M\) by
\[
 ({\cal D}g)_{L,\beta}(\alpha)
   =\sum_{[\delta]:L(\delta)^2=\beta}g_{L,\delta}(\alpha).
                                                               \tag{12g}
\]
Parallel profiles map to zero. The map \({\cal D}\) is onto, its kernel is
exactly the Type-P/Type-K profile lattice, and direct row evaluation gives
\[
                         R_0{\cal T}=p{\cal D}.                \tag{12h}
\]
Since \({\cal T}\) is injective and \(M\) is torsion-free,
\[
 {\cal T}{\cal U}\cap\ker R_0
   ={\cal T}(\ker{\cal D})={\cal K}_{\rm ridge}.               \tag{12i}
\]

Put
\[
 S_0={m(m-1)(4m+1)\over6}.
\]
Proposition 15.760 proves
\(\operatorname {coker}R_0\cong(\mathbf Z/p\mathbf Z)^{S_0}\).
On the quotients
\({\cal U}/\ker{\cal D}\cong M\) and
\(E_0/\ker R_0\cong\operatorname {im}R_0\), equation (12h) says that the
map induced by \({\cal T}\) has image \(pM\). Therefore
\[
 [\operatorname {im}R_0:pM]=p^{r_M-S_0}.
\]
Factoring the index in (12f) through (12i) now gives
\[
 A=\nu_p+r_M-S_0.
\]
Using \(p-1=2m\) simplifies this to the exact formula
\[
 \boxed{\quad
 \nu_p=dpm^2+{m(m-1)(4m+1)\over6},\qquad
 \ker_{\mathbf Z}R/{\cal K}_{\rm ridge}
 \cong(\mathbf Z/p\mathbf Z)^{\nu_p}.
 \quad}                                                       \tag{12j}
\]

In particular, the ridge lattice is proper for every odd prime. A basis of
the mod-\(p\) dependency space produces explicit saturating moves
\(B_pc/p\); together with the ridge columns, those moves generate the
entire integral kernel. This is an exact finite-linear-algebra reduction,
not a claim that the resulting full-kernel generating set is a Graver
basis or is practical to materialize at large \(p\).

More explicitly, choose dependency lifts \(c_1,\ldots,c_{\nu_p}\) whose
reductions form a basis of \(\ker\overline {B_p}\), and put
\(v_j=B_pc_j/p\). For any one integral lift \(z_0\), the Boolean gate is
exactly
\[
 \begin{split}
 &(z_0+\ker_{\mathbf Z}R)\cap\prod_e\{0,\tau_e\}\ne\varnothing\\
 &\quad\Longleftrightarrow\quad
 \exists\,a\in\{0,\ldots,p-1\}^{\nu_p},\
 q\in\mathbf Z^{dmp(p-1)}:
 z_0+\sum_ja_jv_j+B_pq\in\prod_e\{0,\tau_e\}.
                                                               \tag{12k}
 \end{split}
\]
Thus (12j)--(12k), rather than an unspecified integer nullspace, are the
new explicit residual invariant and exact fibre parametrization.

There is also a sharp warning for the full Graver route. Graver elements
generate \(\ker_{\mathbf Z}R\), so their images must generate the quotient
in (12j). Since the Graver basis is symmetric and \(p\) is odd, it must
contain at least \(\nu_p\) non-ridge sign pairs, hence at least
\(2\nu_p\) elements outside \({\cal K}_{\rm ridge}\). The displayed
linearly supported circuits are therefore only a rigorously identified
part of the complete Graver system.

## 4. Constructive defect descent inequalities

Use the signed source and defect from
evidence/NOTE_2026-09-03_EDGE_RADON_SIGNED_BOOLEAN_DEFECT.md:

\[
 \beta(z)={1\over2}\sum_ez_e(z_e-\tau_e).
\]

Since the ridge moves are kernel moves, \(\tau\cdot u=0\), and

\[
             \beta(z+u)-\beta(z)=z\cdot u+{\|u\|^2\over2}.   \tag{13}
\]

For fixed \(L,[\delta]\), put

\[
 S_\delta(\alpha)=\sum_{a:L(a)=\alpha}z(a,[\delta]).          \tag{14}
\]

If \(z\) is a global defect minimizer in its integral fibre, applying
(13) to both signs of the elementary Type-P and Type-K moves gives the
necessary inequalities

\[
 |S_\delta(\alpha)-S_\delta(\alpha')|\le p
 \quad\bigl(L(\delta)=0\bigr),                                \tag{15}
\]

and

\[
 \left|
   S_{\delta_1}(\alpha)-S_{\delta_1}(\alpha')
  -S_{\delta_2}(\alpha)+S_{\delta_2}(\alpha')
 \right|\le2p                                                 \tag{16}
\]

whenever the two difference classes are distinct and have the same
nonzero projected square. If (15) or (16) fails, the appropriate sign of
the displayed ridge is an explicit kernel move that decreases the
nonnegative integer \(\beta\) by at least one.

These are genuine Graver inequalities for Type P, and also for Type K when
\(p\ge5\). They are new necessary conditions on a closest integral lift,
but they are not sufficient for global optimality: the displayed Graver
elements have not been proved to be the complete Graver basis, and the
ridge lattice's elementary \(p\)-torsion quotient may be nontrivial. In
particular, (5) does not turn an unrestricted integral lift into a signed
Boolean lift.

## Exact scope

The theorem supplies an all-prime, linearly supported, constructive kernel
system; an exact \(p\)-saturation statement; and explicit descent tests. It
does **not** prove
the complete Graver basis beyond the displayed circuit families, show that
the compact/all-equal target has defect zero, or close residual (ii). In
fact, (12j) proves that the ridge lattice is strictly smaller than the
integer kernel and hence that the displayed ridge families cannot be the
complete Graver basis.

Executable exact-arithmetic checks are in
src/e1_gmin_m4_ridge_kernel.py and tests/test_ridge_kernel.py.
