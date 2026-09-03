# Exact signed-Boolean defect and Graver alternative

**Status:** proved exact nonlinear reduction; the compact/all-equal target
still has unknown defect, so residual (ii) remains open.

This note starts precisely after Proposition 15.760. Let $R$ denote its
unsigned integral edge--Radon map and retain the signed source convention in
(15.760.10): a graph $H$ has source vector

$$
                    z_e=\tau_e\mathbf 1_{e\in H},
$$

and the desired box is $z_e\in\{0,\tau_e\}$. This is not another Smith,
parity, real least-norm, or local-switch calculation.

## 1. The missing linear term is already a Radon row functional

Every affine-plane edge $e=\{a-d,a+d\}$ is parallel to the unique
projective direction $L_e$ satisfying $L_e(d)=0$. If
$\epsilon_L\in\{\pm1\}$ is the Paley sign of direction $L$, then the edge
sign is

$$
                         \tau_e=\epsilon_{L_e}.              \tag{1}
$$

This is also the identity behind Proposition 15.759's normalized parallel
entry $\epsilon_L\tau_e n_e=n_e$.

Give the parallel target row $P_L$ weight $\epsilon_L$, and give every
off-diagonal row weight zero. Each edge column meets exactly one parallel
row, so (1) gives the integral row identity

$$
                  \tau^t=\sum_L\epsilon_L P_L^t.             \tag{2}
$$

Consequently, on every integral fibre $Rz=y$,

$$
 H_y:=\sum_L\epsilon_LP_L(y)=\tau\mathbin\cdot z             \tag{3}
$$

is fixed. For a genuine graph source, (1) shows that every selected edge
contributes one to the right side, and hence $H_y=|H|$.

## 2. An exact integer quadratic defect

For any integer $n$ and sign $s$, put $w=sn$. Then

$$
 {n(n-s)\over2}={w(w-1)\over2}\ge0,                          \tag{4}
$$

because two consecutive integers have nonnegative product. Equality in (4)
holds exactly when $w\in\{0,1\}$, or equivalently
$n\in\{0,s\}$.

Apply (4) coordinatewise and use (3). For an integral lift $Rz=y$, define

$$
 \beta(z):={1\over2}\sum_ez_e(z_e-\tau_e)
           ={\|z\|_2^2-H_y\over2}.                           \tag{5}
$$

Thus $\beta(z)$ is a nonnegative integer and

$$
 \boxed{\quad z\in\prod_e\{0,\tau_e\}
       \quad\Longleftrightarrow\quad\beta(z)=0.\quad}        \tag{6}
$$

The exact invariant of the fibre is therefore

$$
 \boxed{\quad
 \beta_R(y)=\min_{z\in\mathbf Z^E,\ Rz=y}\beta(z)
            ={1\over2}\left(
              \min_{Rz=y}\|z\|_2^2-H_y\right).
 \quad}                                                       \tag{7}
$$

The minimum exists because the squared norm is coercive on the discrete
affine lattice. Equations (6)--(7) prove

$$
 (z_0+\ker_{\mathbf Z}R)\cap\prod_e\{0,\tau_e\}\ne\varnothing
 \quad\Longleftrightarrow\quad \beta_R(y)=0.                 \tag{8}
$$

If the Boolean fibre is empty, integrality supplies the sharp next gap

$$
                    \min_{Rz=y}\|z\|_2^2\ge H_y+2.           \tag{9}
$$

Proposition 15.761 minimized the norm over the *real* fibre. Its exact
Moore--Penrose value can lie strictly below $H_y$ on the compact rays.
The new quantity in (7) is the closest-vector norm in the integral fibre;
(9) identifies the exact threshold that a successful integer argument must
cross. It does not claim that the threshold is crossed.

## 3. Exact Graver/Voronoi alternative

For integer vectors write $g\sqsubseteq d$ when they are in the same
closed orthant and $|g_e|\le|d_e|$ in every coordinate. The Graver basis
$\mathcal G(R)$ consists of the nonzero $\sqsubseteq$-minimal vectors of
$\ker_{\mathbf Z}R$. It is finite and symmetric. Repeatedly splitting a
nonminimal kernel vector proves directly that every
$d\in\ker_{\mathbf Z}R$ is a conformal sum

$$
                         d=g_1+\cdots+g_k,
 \qquad g_i\in\mathcal G(R).                                \tag{10}
$$

Equation (2) implies $\tau\cdot g=0$ for every kernel move. Expanding (5)
therefore gives

$$
 \beta(z+g)-\beta(z)=z\mathbin\cdot g+{\|g\|_2^2\over2}.    \tag{11}
$$

For the conformal sum (10), all cross-products
$g_i\cdot g_j$ are nonnegative, and hence

$$
 \begin{aligned}
 \beta(z+d)-\beta(z)
 &=z\cdot d+{\|d\|^2\over2}\\
 &\ge\sum_i\left(z\cdot g_i+{\|g_i\|^2\over2}\right)
  =\sum_i\bigl(\beta(z+g_i)-\beta(z)\bigr).                 \tag{12}
 \end{aligned}
$$

If any kernel move improves $z$, (12) forces one Graver summand to improve
it. Conversely, a global minimizer plainly admits no improving Graver move.
Thus

$$
 z\text{ minimizes (7)}
 \quad\Longleftrightarrow\quad
 \beta(z+g)\ge\beta(z)\quad(g\in\mathcal G(R)).              \tag{13}
$$

Since the basis contains both $g$ and $-g$, (11) rewrites (13) as the
Voronoi inequalities

$$
                  \boxed{|2z\cdot g|\le\|g\|_2^2
                         \quad(g\in\mathcal G(R)).}           \tag{14}
$$

Starting from Proposition 15.760's integral lift, any strict Graver descent
decreases the nonnegative integer $\beta$, so it terminates. The terminal
alternative is exact:

- $\beta=0$, in which case (6) constructs the common graph;
- $\beta>0$ together with all inequalities (14), in which case (13)
  certifies that the Boolean fibre is empty.

This is a mathematical finite alternative, not a claim that a practical
complete Graver basis of the full matrix has been generated. The switch
support theorem in
evidence/NOTE_2026-09-02_EDGE_RADON_SWITCHING_DEGREE.md proves that every
nonzero move already has support at least $p+1$, so no bounded local move
list can implement (14).

## 4. The binary theorem removes parity but not the defect

Proposition 15.760 gives

$$
 \operatorname {coker}R\cong(\mathbf Z/p\mathbf Z)^{S(p)},
$$

with $p$ odd. Therefore this cokernel has neither 2-torsion nor a
nonzero quotient modulo two. Reducing the exact lattice sequence modulo two
shows

$$
 \ker_{\mathbf Z}R\longrightarrow\ker_{\mathbf F_2}R_2
 \quad\hbox{is onto}.                                       \tag{15}
$$

Indeed, factor $R$ through its image. Tensoring
$0\to\ker R\to\mathbf Z^E\to\operatorname {im}R\to0$ and
$0\to\operatorname {im}R\to\mathcal A\to\operatorname {coker}R\to0$
with $\mathbf F_2$ is exact at the required spots because both
$\operatorname {Tor}_1(\operatorname {coker}R,\mathbf F_2)$ and
$\operatorname {coker}R\otimes\mathbf F_2$ vanish.

Hence any binary preimage furnished by Proposition 15.757 can be matched by
an integral lift with the same coordinate parities. What remains is exactly
the sign-and-magnitude defect (7), not another mod-two equation.

## 5. Exact pointwise formula for the canonical fractional lift

Proposition 15.761 gives the three spectral blocks of $RR^t$. They also
give a coordinate formula which its scalar norm does not display. Use the
unsigned target convention

$$
 y=(\widetilde P_L,\widetilde K_L),\qquad
 \widetilde P_L+\sum_c\widetilde K_L(c)=T,
 \qquad \sum_L\widetilde P_L=T.
$$

For an edge $e$, let $L_e$ be its parallel direction and let $c_L(e)$ be
its off-diagonal cell in direction $L\ne L_e$. Put
$C=\binom p2$, $d=p+1$, and

$$
 k_L={T-\widetilde P_L\over C},\qquad
 q_L=\widetilde P_L-{T\over d},\qquad
 w_L=\widetilde K_L-k_L\mathbf1.
$$

The within-direction block, directional block, and uniform block contribute
to $z^*=R^+y$, respectively,

$$
 {1\over p^2}\sum_{L\ne L_e}w_L(c_L(e)),\qquad
 {q_{L_e}\over pC},\qquad
 {T\over\binom{p^2}{2}}.
$$

Substitution of $k_L,q_L$, using
$\sum_L\widetilde P_L=T$ and
$\binom{p^2}{2}=dpC$, simplifies all constant terms to

$$
 \boxed{\quad
 z_e^*={2(\widetilde P_{L_e}-T)\over p^3}
       +{1\over p^2}\sum_{L\ne L_e}
          \widetilde K_L(c_L(e)).
 \quad}                                                       \tag{16}
$$

In residual notation,

$$
 \widetilde P_L=\epsilon_LP_L,\qquad
 \widetilde K_L=\epsilon_LW_L,\qquad
 \tau_e=\epsilon_{L_e}.
$$

Define the pointwise signed backprojection

$$
 B_e:=\sum_{L\ne L_e}\tau_e\epsilon_LW_L(c_L(e)).             \tag{17}
$$

Multiplying (16) by $\tau_e$ gives the physical fractional coordinate

$$
 \boxed{\quad
 h_e^*:=\tau_e z_e^*
   ={2(P_{L_e}-\tau_eT)+pB_e\over p^3}.
 \quad}                                                       \tag{18}
$$

Since $R^+y$ is an exact real lift for every compatible target, (18) yields
the structural sufficient theorem

$$
 0\le2(P_{L_e}-\tau_eT)+pB_e\le p^3
 \quad\hbox{for every edge }e
 \quad\Longrightarrow\quad
 y\in R\prod_e[0,\tau_e].                                    \tag{19}
$$

This is a genuine pointwise invariant beyond Proposition 15.761's aggregate
$\ell^2$ norm. It is sufficient, not necessary: if the Moore--Penrose lift
leaves the box, another real lift obtained by adding a kernel vector may
still enter it. For the compact/all-equal target, the signs of $B_e$ are not
controlled by the current atom-count or $\ell^1$ bounds. Thus (19) isolates
a concrete next calculation but does not yet settle its fractional box.

## 6. The periodic dual collapses exactly to the fractional-box gate

There is also a one-sided dual certificate that retains integer rounding.
For a real target multiplier $\lambda$, put $c=R^t\lambda$. From (5) and
$Rz=y$, coordinatewise minimization over $z_e\in\mathbf Z$ gives

$$
 \beta_R(y)\ge D_y(\lambda):=
 \lambda\cdot y+{1\over2}\sum_e\left[
 \operatorname {dist}\left(c_e+{\tau_e\over2},\mathbf Z\right)^2
 -\left(c_e+{\tau_e\over2}\right)^2\right].                 \tag{20}
$$

Thus any $\lambda$ with $D_y(\lambda)>0$ rigorously excludes a Boolean
lift. It does **not**, however, see an integer obstruction after the
fractional signed box becomes feasible.

Indeed, the coordinatewise minimum defining (20) includes the two choices
$z_e=0,\tau_e$, both of defect zero. Therefore

$$
 D_y(\lambda)\le
 F_y(\lambda):=\lambda\cdot y
   -\sum_e\max\{0,\tau_e(R^t\lambda)_e\}.                    \tag{21}
$$

The right side is precisely the signed-box Farkas slack. If
$y\in R\prod_e[0,\tau_e]$, then $F_y(\lambda)\le0$ for every $\lambda$,
so $D_y(\lambda)\le0$ as well; $\lambda=0$ makes both values zero.

Conversely, if the fractional signed box is empty, separation supplies a
$\lambda$ with $F_y(\lambda)>0$. Scale it by a sufficiently small positive
$s$ so that

$$
               |s(R^t\lambda)_e|<1\qquad\hbox{for every }e.
$$

For $|c_e|<1$, the integer minimum in (20) is attained by one of
$0,\tau_e$, and its value is exactly
$-\max(0,\tau_ec_e)$. Consequently

$$
                       D_y(s\lambda)=F_y(s\lambda)>0.         \tag{22}
$$

Equations (21)--(22) prove the exact method verdict

$$
 \sup_\lambda D_y(\lambda)>0
 \quad\Longleftrightarrow\quad
 R\prod_e[0,\tau_e]\ \hbox{does not contain }y.              \tag{23}
$$

Thus the periodic dual is another expression of the zonotope/Farkas gate
already isolated in the cohomology/topology audit. It can exclude a target
whose *fractional* box is empty, but it cannot certify a positive
$\beta_R(y)$ caused only by integrality. For example, the executable
rank-one-kernel test has $\beta_R(y)=1$ and no Boolean lift, while
$(1/2,1/2,0)$ is a fractional lift; every value in (20) is then nonpositive.
The Graver alternative (13)--(14), not (20), is the genuinely integer part
of this note.

## Exact scope

This note proves an exact nonlinear reformulation and a complete Graver
alternative. It does **not** prove
$\beta_R(y)=0$ for the compact/all-equal atoms, decide whether their
fractional zonotope fibre is nonempty, construct the full Graver basis, or
close residual (ii). The live integer target is now precise: after the
moment congruences and fractional-box inequalities pass, prove defect zero
or prove a positive integer defect using (14).

Executable identity checks are in
src/e1_gmin_m4_signed_boolean_defect.py and
tests/test_signed_boolean_defect.py.
