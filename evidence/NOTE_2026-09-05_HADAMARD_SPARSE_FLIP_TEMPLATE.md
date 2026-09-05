# Actual sign matrices with a finite template on a Hadamard background

2026-09-05. Analytic construction. This realizes a flat singular-value
bulk, finite operator outliers, and an asymptotically optimal scalar
SDP dual by actual complete cross sign matrices. It leaves a precisely
stated Boolean-norm gap between a lower bound and a magnetization
completion upper. It does not construct an original conditional
optimizer or prove compatibility with a particular original covariance.

No matrix was sampled, enumerated, or numerically optimized to obtain
this existence theorem. The only probabilistic input is the published
matrix Bernstein theorem stated explicitly below.

## 1. Fixed template and regular Hadamard background

Fix p=4^ell with ell>=0, and a real p by p matrix C. Write

\[
 q=\|C\|_{\rm op},\qquad
 \beta(C)=\max_{a,b\in\{-1,1\}^p}|a^TCb|,
\]
\[
 \tau(C)=\max_{\|u_i\|=\|v_j\|=1}
                   \sum_{ij}C_{ij}\langle u_i,v_j\rangle.
                                                               \tag{1}
\]

The construction itself permits arbitrary C. The scalar-SDP conclusion
will assume q>=1 and the actual finite-template equality

\[
                              \tau(C)=pq.                 \tag{2}
\]

No Boolean/SDP ratio or sign condition on the entries of C is assumed.
Condition (2) must be proved for a proposed template; it is not inferred
from its operator norm.

Let F=J_4-2I_4, where J_4 is the all-ones matrix. For r=4^j put

\[
 H_r=F^{\otimes j},\qquad O_r=H_r/\sqrt r,
 \qquad H_1=O_1=(1).
\]

Direct multiplication gives FF^T=4I_4 and F1_4=2 1_4. Therefore
each H_r is a real sign matrix and

\[
 O_rO_r^T=I_r,\quad O_r^T=O_r,\quad O_r1_r=1_r.          \tag{3}
\]

Let m=4^k with k>=1 and n=pm. The order-n background is
H_n=H_p tensor H_m and O_n=O_p tensor O_m. Define the block-constant
isometry J_n from R^p to R^n by

\[
                  (J_n a)_{(i,\alpha)}=a_i/\sqrt m.
\]

Its range S_n is the subspace constant on each of the p blocks.
Equations (3) imply O_n J_n=J_n O_p; S_n and its orthogonal complement
are invariant for O_n. Define the deterministic normalized model

\[
              T_n=O_n+J_n(C-O_p)J_n^T.                   \tag{4}
\]

Relative to S_n direct-sum S_n-perp, this is exactly C direct-sum
the restriction of the orthogonal O_n to S_n-perp. Hence its singular
multiset is

\[
       \{\sigma_1(C),\ldots,\sigma_p(C)\}
                     \ \cup\ \{1\text{ repeated }n-p\text{ times}\},
 \qquad \|T_n\|_{\rm op}=\max(1,q).                       \tag{5}
\]

## 2. Sparse independent flips and their exact mean

Set

\[
 M=p(C-O_p),\qquad \gamma=\max(1,\max_{ij}|M_{ij}|).
\]

For n>=gamma^2, independently flip the background entry
`H_n[(i,alpha),(j,beta)]` with probability

\[
 \pi_{i\alpha,j\beta}
 ={\gamma-H_n[(i,\alpha),(j,\beta)]M_{ij}\over2\sqrt n}.
                                                               \tag{6}
\]

Every probability lies in [0,gamma/sqrt(n)] and hence in [0,1].
Let B_n be the resulting actual n by n sign matrix. Its mean is

\[
 \mathbb E B_n
 =(1-\gamma/\sqrt n)H_n+{M\otimes J_m\over\sqrt n},
 \qquad
 {\mathbb E B_n\over\sqrt n}
                         =T_n-{\gamma\over\sqrt n}O_n.   \tag{7}
\]

Here J_m in the tensor product is the all-ones m by m matrix, not the
isometry J_n. Equation (7) follows entry by entry from the sign flip,
using `H_n[entry]^2=1`, and from n=pm.

Let R_n=B_n-E B_n. Its independent entries are centered, have absolute
value at most two, and satisfy

\[
 \mathbb E(R_n)_{ab}^2
       =4\pi_{ab}(1-\pi_{ab})\le4\gamma/\sqrt n.
                                                               \tag{8}
\]

In particular both maximum row variance sum and maximum column
variance sum are at most 4 gamma sqrt(n).

## 3. Published concentration input and deterministic existence

We use the following precise version of matrix Bernstein. If Z_j are
independent centered real d_1 by d_2 random matrices, ||Z_j||op<=R
almost surely, and

\[
 v=\max\left(\left\|\sum_j\mathbb E Z_jZ_j^T\right\|_{\rm op},
             \left\|\sum_j\mathbb E Z_j^TZ_j\right\|_{\rm op}\right),
\]

then, for every t>=0,

\[
 \Pr\left\{\left\|\sum_j Z_j\right\|_{\rm op}\ge t\right\}
       \le(d_1+d_2)\exp\left[-{t^2\over2(v+Rt/3)}\right]. \tag{9}
\]

This is Theorem 1.6 in Joel Tropp's author-hosted
[User-Friendly Tail Bounds for Sums of Random Matrices](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf),
published in Foundations of Computational Mathematics 12 (2012),
389--434, [DOI 10.1007/s10208-011-9099-z](https://doi.org/10.1007/s10208-011-9099-z).
The rectangular theorem's entire statement and its variance convention
were checked in the primary source; no unproved iid-variance surrogate
is substituted for v.

Apply (9) to the individual centered-entry matrices of R_n. Their
norm bound is R=2, their dimensions are n by n, and their variance
parameter is at most 4 gamma sqrt(n), by (8). Thus

\[
 \Pr\{\|R_n\|_{\rm op}\ge n^{3/8}\}
 \le 2n\exp\left[-{n^{3/4}\over
                    2(4\gamma\sqrt n+(2/3)n^{3/8})}\right]
 \le 2n\exp\left[-{n^{1/4}\over8\gamma+4/3}\right]
                         =:b_n\longrightarrow0.         \tag{10}
\]

Let F_n count the actual flips. Equation (6) gives
E F_n<=gamma n^(3/2), so Markov's inequality gives

\[
                  \Pr\{F_n>4\gamma n^{3/2}\}\le1/4.     \tag{11}
\]

For every sufficiently large n in the indicated sequence, b_n<1/4.
The events `||R_n||op<=n^(3/8)` and `F_n<=4 gamma n^(3/2)` therefore
intersect with probability at least 1/2. Choose any realization in
their intersection. This supplies a deterministic sequence of complete
sign matrices B_n with only O_C(n^(3/2)) changed entries and

\[
 \boxed{\displaystyle
 \left\|{B_n\over\sqrt n}-T_n\right\|_{\rm op}
 \le\varepsilon_n:=n^{-1/8}+\gamma n^{-1/2}\longrightarrow0.}
                                                               \tag{12}
\]

The existence assertion does not prescribe an efficient deterministic
search for the selected realizations, and no such search was run.

## 4. Actual singular bulk and asymptotically optimal scalar dual

Weyl's singular-value perturbation inequality and (5), (12) identify
the singular spectrum of B_n/sqrt(n) to uniform error epsilon_n.
Consequently its empirical singular-value law converges to delta_1.
Its finitely many non-bulk limiting singular values are those of the
fixed template C; values equal to one are of course indistinguishable
from the bulk. In particular, when q>=1,

\[
       q-\varepsilon_n\le\|B_n\|_{\rm op}/\sqrt n
                                      \le q+\varepsilon_n. \tag{13}
\]

Assume now (2). For any n by n real matrix W,

\[
                            \tau(W)\le n\|W\|_{\rm op},
\]

by Cauchy--Schwarz on its two unit-row vector families. Moreover
`|tau(W)-tau(W')|<=n||W-W'||op` by the same argument.

The deterministic matrix sqrt(n) T_n has SDP value exactly
q n^(3/2). The operator bound proves the upper. For the matching
lower, take optimal unit vectors u_i,v_j for C and repeat each
vector throughout its block. The sum of entries of T_n in block
(i,j) is m C_ij, so this primal has value

\[
             \sqrt n\,m\,\tau(C)=q n^{3/2}.
\]

Equation (12) therefore yields

\[
 \boxed{\displaystyle
 q-\varepsilon_n\le{\tau(B_n)\over n^{3/2}}
                                   \le q+\varepsilon_n.} \tag{14}
\]

The scalar diagonal `D_n=||B_n||op I_(2n)` is feasible in the
bipartite SDP dual, since `D_n +- H_(B_n)>=0`. Equations (13)--(14)
give

\[
 0\le n\|B_n\|_{\rm op}-\tau(B_n)
                      \le2\varepsilon_n n^{3/2},\qquad
 {\tau(B_n)\over n\|B_n\|_{\rm op}}\longrightarrow1.     \tag{15}
\]

This is ASYMPTOTIC scalar-dual optimality. It does not assert that an
exact optimal diagonal of each realized B_n is scalar, or that the
exact finite-n optimal Gram has a prescribed rank or eigenspace.

## 5. The correctly directed Boolean-norm sandwich

Define the finite-template magnetization completion

\[
 \Gamma(C)=\max_{a,b\in[-1,1]^p}\left[
 {a^TCb\over p}
 +\sqrt{\left(1-{\|a\|_2^2\over p}\right)
        \left(1-{\|b\|_2^2\over p}\right)}\right].        \tag{16}
\]

The maximum exists by compactness. Changing a to -a shows that the
bilinear term can equivalently be replaced by its absolute value.
For q>=1, direct Cauchy--Schwarz gives

\[
                  \max(1,\beta(C)/p)\le\Gamma(C)\le q.   \tag{17}
\]

For the last bound put r=||a||/sqrt(p), s=||b||/sqrt(p). The
expression is at most `q r s+sqrt((1-r^2)(1-s^2))`, which is at
most `q r s+1-rs<=q`.

The actual complete sign matrices constructed above satisfy

\[
 \boxed{\displaystyle
 \max(1,\beta(C)/p)-\varepsilon_n
 \le{\beta(B_n)\over n^{3/2}}
 \le\Gamma(C)+\varepsilon_n.}                            \tag{18}
\]

To prove the upper, let x,y be any actual sign vectors and let a_i,
b_j be their block means. Decompose them into their block-constant
parts and perpendicular parts x_perp,y_perp. Equation (4) gives

\[
 {x^TT_n y\over n}
 ={a^TCb\over p}+{x_\perp^TO_n y_\perp\over n},
\]
\[
 \|x_\perp\|_2^2=n\left(1-{\|a\|_2^2\over p}\right),
 \qquad
 \|y_\perp\|_2^2=n\left(1-{\|b\|_2^2\over p}\right).
\]

Orthogonality of O_n bounds the perpendicular term by the square
root in (16). Taking the absolute maximum and then using (12)
proves the upper in (18).

For the lower term beta(C)/p, use signs constant throughout each
block and choose their p coordinates to attain beta(C). Their
perpendicular components vanish, giving exactly beta(C)/p in the
normalized deterministic model.

The background lower term one requires an actual balanced witness,
not an assumption about the background maximum. Put

\[
 w_4=(1,1,-1,-1)^T,\qquad
 w_m=w_4\otimes1_{m/4}.
\]

Directly `1_m^T w_m=0` and `H_m w_m=-sqrt(m) w_m`, using
`F w_4=-2 w_4` and regularity of the other tensor factors. Take

\[
                 y=1_p\otimes w_m,\qquad x=-y.
\]

These are actual sign vectors with zero mean in every block.
Their template components vanish. Since O_p 1_p=1_p and
O_m w_m=-w_m, one has T_n y=O_n y=x and therefore
`x^T T_n y=n`. This proves the exact normalized background value one.
Finally (12) changes any sign-pair value by at most epsilon_n after
normalization, proving both lower terms in (18).

Equivalently, (18) implies the liminf/limsup sandwich

\[
 \max(1,\beta(C)/p)
 \le\liminf_n{\beta(B_n)\over n^{3/2}}
 \le\limsup_n{\beta(B_n)\over n^{3/2}}\le\Gamma(C).        \tag{19}
\]

No equality between the actual Boolean limit and Gamma(C) is claimed.
The square-root completion allowed arbitrary perpendicular vectors;
actual signs need not realize equality in that Cauchy--Schwarz step.

## 6. What this realizes and what still requires a theorem

This construction shows that a delta_1 bulk, finitely many bounded
sqrt(n)-scale outliers, and asymptotically scalar SDP optimality can
coexist in actual complete sign matrices. They are not by themselves
an obstruction to attainability.

For a proposed conditional-cap counterprofile at target f, the
sufficient finite-template conditions
`||C||op=q>=1`, `tau(C)=pq`, and `Gamma(C)<=f` would give actual
sign matrices with asymptotic SDP value q n^(3/2) and Boolean norm
at most f n^(3/2). No template satisfying the problematic sharp
parameters is supplied here. A small value of beta(C)/p alone is
not sufficient because Gamma(C) may be strictly larger.

Conversely, proving Gamma(C)>f for a class of templates only defeats
this PARTICULAR completion-based upper certificate. It does not by
itself prove that the actual B_n have Boolean norm above f n^(3/2),
because Gamma is an upper, not a lower, in (18).

The original internal signing, actual conditional minimization,
intrinsic source normalization, and joint-shell compatibility remain
separate requirements. No assertion about them, the final Gaussian
comparison, or the original MO limit follows from this construction.
