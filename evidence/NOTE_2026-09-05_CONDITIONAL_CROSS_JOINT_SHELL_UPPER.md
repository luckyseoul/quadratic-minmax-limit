# Conditional cross-covariance: a joint-shell Gaussian upper

2026-09-05. This is an upper bound for the actual conditional cross
proposal. It retains the exchange covariance and the independent Gaussian
cushion. The leading comparison with the smaller-order optimum is not
proved here.

## 1. The source, conditional optimizer, and covariance

Let n>=2. Let A be a complete symmetric zero-diagonal order-n signing,
and let B be any n by n sign matrix. Set

\[
 K=\begin{pmatrix}A&B\\B^T&-A\end{pmatrix},\qquad
 L=\|K\|_{\rm op},\qquad D=L^2-1.
\]

Thus `L^2>=2n-1>=3`, and D is positive. In cross coordinates define

\[
 (\mathcal S_B)_{ij,kl}=B_{il}B_{kj},\qquad
 R={L^2 I+A\otimes A-\mathcal S_B\over D}.             \tag{1}
\]

The tensor convention in (1) means that its `ij,kl` entry is
`A_ik A_jl`. Equivalently `S_B(X)=B X^T B`.

One has `0<=R<=3I` and `diag R=1`. To see positivity directly, compress
the self-adjoint operator `X -> K X K` to the orthonormal symmetric
cross-edge matrices `(E_(i,n+j)+E_(n+j,i))/sqrt(2)`. Its compression has
entries `-A_ik A_jl+B_il B_kj`, and its operator norm is at most L^2.
Consequently R is PSD and has norm at most `2L^2/D<=3`. The compression
has diagonal one, proving the diagonal assertion.

Let k,v>=0 and let Z be a centered Gaussian n by n matrix with

\[
 \operatorname{Cov}(\operatorname{vec}Z)=C=kR+vI.        \tag{2}
\]

The actual shifted-sign linearized proposal has
`k=4 phi(h)^2`, `v=1-s_h^2-k`, and deterministic cross mean `s_h B`.
The proof below allows every nonnegative k and v. In particular, the
independent cushion is not removed.

For a nonempty joint shell put

\[
 \mathcal T_{p,q,c}=
 \{(x,y)\in\{-1,1\}^{2n}:x^TAx=p,\ y^TAy=q,\ x^TBy=c\}.
                                                               \tag{3}
\]

No optimization assumption on B is needed for the upper bound. The
conditional minimizing property of B may subsequently be combined with
the separately proved Gaussian-refill floor.

## 2. A Gaussian comparison retaining the full exchange term

Write

\[
 a_0=n\left({kL^2\over D}+v\right),\qquad t={k\over D}.
                                                               \tag{4}
\]

Let `(xi,eta)` be a centered Gaussian vector in R^(2n) with covariance

\[
 M_{p,q,c}=
 \begin{pmatrix}
 a_0 I+tqA&-tcB\\
 -tcB^T&a_0 I+tpA
 \end{pmatrix}.                                               \tag{5}
\]

This is PSD. It suffices to check k=1,v=0. The compression of K to the
orthogonal vectors `(x,0)/sqrt(n)` and `(0,y)/sqrt(n)` is

\[
 {1\over n}\begin{pmatrix}p&c\\c&-q\end{pmatrix},
\]

so `H=[[q,-c],[-c,-p]]` has norm at most nL: it is minus the orthogonal
conjugate of the displayed unnormalized 2 by 2 matrix by the coordinate
swap `J=[[0,1],[1,0]]`. The nonconstant block
matrix in `D M` is a coordinate compression of `H tensor K`, hence
has norm at most nL^2. Adding `nL^2 I` proves positivity. For general
k,v, add the independent covariance `vn I_(2n)` and multiply the first
covariance by k.

Then the following is a genuine upper comparison:

\[
 \boxed{\displaystyle
 \mathbb E\sup_{(x,y)\in\mathcal T_{p,q,c}}x^TZy
 \le
 w_{p,q,c}:=
 \mathbb E\sup_{(x,y)\in\mathcal T_{p,q,c}}
                  (\xi^Tx+\eta^Ty).}                         \tag{6}
\]

Here are the complete increment calculations. First take k=1,v=0 and
two states in the same shell. Abbreviate

\[
 r_x=x^Tx',\quad r_y=y^Ty',\quad
 a=x^TAx',\quad b=y^TAy',\quad
 d=x^TBy',\quad e=x'^TBy.
\]

For `X_(x,y)=x^TZy` and `Y_(x,y)=xi^Tx+eta^Ty`,

\[
 \operatorname{Cov}(X,X')={L^2r_xr_y+ab-de\over D},
 \qquad \operatorname{Var}X={L^2n^2+pq-c^2\over D},             \tag{7}
\]

\[
 \operatorname{Cov}(Y,Y')=
 {L^2n(r_x+r_y)+qa+pb-c(d+e)\over D},\qquad
 \operatorname{Var}Y=2\operatorname{Var}X.                     \tag{8}
\]

In particular the exchange term in (7) is `-de`, not `-c c'`.
If `Delta=Var(Y-Y')-Var(X-X')`, direct subtraction gives

\[
 {\Delta\over2}
 ={L^2(n-r_x)(n-r_y)+(p-a)(q-b)-(c-d)(c-e)\over D}
                                                               \tag{9}
\]

\[
 ={D\,\langle\delta x\otimes\delta y,
                   R(\delta x\otimes\delta y)\rangle
                  +(d-e)^2\over4D}\ \ge 0,
 \qquad \delta x=x-x',\quad\delta y=y-y'.                     \tag{10}
\]

For the independent covariance vI, the corresponding linear fields
have covariance `vn I_(2n)` and their increment excess divided by two
is `v(n-r_x)(n-r_y)>=0`. Adding this to k times (10) proves increment
domination for (2). The finite Gaussian comparison theorem, or its
log-sum-exp interpolation proof followed by the zero-temperature limit,
now proves (6). Singularity of either covariance causes no difficulty:
one may use Gaussian factor representations or positive diagonal
regularizations and then pass to the limit.

## 3. Exact reference energies of the maximizing linear-field signs

Suppose a_0>0; if a_0=0 all Gaussian widths in this note are zero. Put

\[
 x^0_i=\operatorname{sign}\xi_i,\qquad
 y^0_j=\operatorname{sign}\eta_j.
\]

Every field coordinate has variance a_0. Gaussian arcsine identities
therefore give the exact three expectations

\[
 \begin{aligned}
 m_p:=\mathbb E[(x^0)^TAx^0]
    &={2\over\pi}n(n-1)\arcsin(tq/a_0),\\
 m_q:=\mathbb E[(y^0)^TAy^0]
    &={2\over\pi}n(n-1)\arcsin(tp/a_0),\\
 m_c:=\mathbb E[(x^0)^TBy^0]
    &=-{2\over\pi}n^2\arcsin(tc/a_0).
 \end{aligned}                                               \tag{11}
\]

In particular m_c has the opposite sign to c when k>0 and c is nonzero.
The actual cushion appears in the argument

\[
 {tc\over a_0}={kc\over n(kL^2+vD)}.                          \tag{12}
\]

No field-coordinate independence is asserted or used.

## 4. A joint mismatch deficit

Define the mismatch vector in quadratic-energy units by

\[
 b=\left({p-m_p\over2},{q-m_q\over2},c-m_c\right).
                                                               \tag{13}
\]

For `lambda=(lambda_1,lambda_2,lambda_3)` define

\[
 H_\lambda=
 \begin{pmatrix}\lambda_1A&\lambda_3B\\
                 \lambda_3B^T&\lambda_2A\end{pmatrix},
 \qquad
 R_\lambda={\lambda_{\max}(H_\lambda)-
                  \lambda_{\min}(H_\lambda)\over2},
\]

\[
 r=\sup_{\lambda\ne0}{|\lambda\cdot b|\over2nR_\lambda},
 \qquad \delta(r)={1-\sqrt{1-r^2}\over2}.                     \tag{14}
\]

For n>=2, nonzero lambda gives nonzero H_lambda and R_lambda>0.
The proof below also shows `0<=r<=1` on every nonempty shell.

Let Phi_G denote the standard Gaussian distribution function and put

\[
 z_d=\Phi_G^{-1}((1+d)/2),\qquad
 \ell(d)=1-e^{-z_d^2/2},\quad 0\le d\le1,
                                                               \tag{15}
\]

with `ell(1)=1`. Then

\[
 \boxed{\displaystyle
 w_{p,q,c}\le
 2n\sqrt{2a_0/\pi}\,[1-2\ell(\delta(r))]
 =2n\sqrt{2a_0/\pi}\,[2e^{-z_{\delta(r)}^2/2}-1].}            \tag{16}
\]

In particular,

\[
 w_{p,q,c}\le 2n\sqrt{2a_0/\pi}
                   -2n\sqrt{\pi a_0/2}\,\delta(r)^2.         \tag{17}
\]

Here is a proof with the joint selection dependence retained. Choose a
measurable maximizing state `z*=(x*,y*)` in the finite shell and let
`z0=(x0,y0)`. Let D_H be their Hamming distance and N=2n. For any
symmetric H, subtracting the midpoint of its two extreme eigenvalues
does not affect the difference of its quadratic forms on Boolean
vectors. The vectors `z*+z0` and `z*-z0` are orthogonal and have norms
`2sqrt(N-D_H)` and `2sqrt(D_H)`. Consequently

\[
 |Q_{H_\lambda}(z^*)-Q_{H_\lambda}(z^0)|
 \le2R_\lambda\sqrt{D_H(N-D_H)}.
\]

Taking expectations and applying concavity of `sqrt(d(N-d))`, with
`dbar=E D_H/N`, gives

\[
 |\lambda\cdot b|\le2NR_\lambda
              \sqrt{\bar d(1-\bar d)},\qquad
 r\le2\sqrt{\bar d(1-\bar d)}\le1.
\]

Thus `dbar>=delta(r)`. For any flip indicators f_i, including indicators
depending on ALL the fields, and any threshold u>=0,

\[
 \sum_i\mathbb E[|g_i|f_i]
 \ge u\sum_i\mathbb E f_i-\sum_i\mathbb E(u-|g_i|)_+,
 \qquad g=(\xi,\eta).
\]

Optimizing the right side at `u=sqrt(a_0) z_dbar` uses only the common
one-coordinate normal marginal and yields

\[
 \sum_i\mathbb E[|g_i|f_i]
 \ge N\sqrt{2a_0/\pi}\,\ell(\bar d).
\]

Since `g dot z*=sum_i |g_i|-2sum_i |g_i|f_i`, (16) follows.
Finally `ell'(d)=sqrt(pi/2) z_d` and
`z_d>=sqrt(pi/2)d`, so `ell(d)>=pi d^2/4`, proving (17).

For example, taking `lambda=(1,-1,1)` in (14) retains the JOINT
internal-plus-cross mismatch and uses `R_lambda<=L`. Taking only
lambda_3 nonzero gives the separate cross mismatch with
`R_lambda=|lambda_3| ||B||op`. Neither is replaced by the source
unperturbed Gibbs energy.

## 5. An optional sharper two-block Hamming bound

Let `R_A=(lambda_max(A)-lambda_min(A))/2` and `L_B=||B||op`.
Consider all `(d_x,d_y)` in [0,1]^2 satisfying

\[
 \begin{aligned}
 |b_1|&\le2nR_A\sqrt{d_x(1-d_x)},\\
 |b_2|&\le2nR_A\sqrt{d_y(1-d_y)},\\
 |b_3|&\le2nL_B\bigl[
          \sqrt{(1-d_x)d_y}+\sqrt{d_x(1-d_y)}\bigr],\\
 r&\le2\sqrt{\bar d(1-\bar d)},\qquad
                  \bar d=(d_x+d_y)/2.
 \end{aligned}                                               \tag{18}
\]

This compact feasible set is nonempty, and

\[
 \boxed{\displaystyle
 w_{p,q,c}\le2n\sqrt{2a_0/\pi}
     \left[1-\min_{(d_x,d_y)\ {m satisfying}\ (18)}
                          (\ell(d_x)+\ell(d_y))\right].}      \tag{19}
\]

Indeed, use `d_x=E Ham(x*,x0)/n` and the analogous d_y. The first two
constraints follow from the same one-block quadratic argument. For the
third, the exact polarization

\[
 x^{*T}By^*-x^{0T}By^0
 ={1\over2}\bigl[(x^*+x^0)^TB(y^*-y^0)
                 +(x^*-x^0)^TB(y^*+y^0)\bigr]
\]

and the two operator-norm bounds give the stated inequality after
expectation: the geometric mean is jointly concave on the nonnegative
quadrant. Apply the threshold argument separately to the two blocks to
prove (19). Convexity of ell shows (19) is at least as strong as (16).

## 6. Returning to the actual Gaussian maximum

There are at most `J=(2n^2+1)^3` nonempty shells. For every deterministic
real s, (6), Gaussian concentration, and the two absolute-value phases
give

\[
 \boxed{\displaystyle
 \mathbb E\Phi\!\begin{pmatrix}A&sB+Z\\sB^T+Z^T&-A\end{pmatrix}
 \le\max_{\mathcal T_{p,q,c}\ne\varnothing}
       \left[\left|{p-q\over2}+sc\right|+w_{p,q,c}\right]
       +n\sqrt{2\|C\|_{\rm op}\log(2J)}.}                   \tag{20}
\]

Any right side of (16), (17), or (19) may replace w in (20).
The concentration term follows because each shell maximum, including
each of its two phases, is n-Lipschitz in the cross Gaussian vector,
whose covariance is C. Its centered logarithmic moment generating
function is at most `lambda^2 n^2 ||C||op/2`; the finite maximum bound
then gives the displayed remainder. Also `||C||op<=3k+v`.

If B is an ACTUAL conditional cross optimizer, then every nonempty shell
additionally satisfies

\[
 \left|{p-q\over2}\right|+|c|\le
 F_A^*:=\min_{B'\in\{-1,1\}^{n\times n}}
       \max_{x,y}\left( |Q_A(x)-Q_A(y)|+|x^TB'y|\right).
                                                               \tag{21}
\]

The masked-cross comparison proved next combines with (20) at
`s=s_h`, `k=4phi(h)^2`, `v=1-s_h^2-k`. This gives a genuine
finite-dimensional self-consistency upper involving the actual
optimizer and its attainable shells. No evaluation proving
`F_A^*<=2sqrt(2) Phi(A)+o(n^(3/2))` is asserted.

## 7. The required masked-cross Gaussian-refill floor

For precision, this step does NOT apply a theorem stated only for a
whole quadratic observable to a silently masked observable. Let
`G_R~N(0,R)` on the n^2 cross coordinates in (1), let h be any fixed
real threshold, and define the ACTUAL cross signing

\[
 \widehat B_{h,ij}=\operatorname{sign}(G_{R,ij}+hB_{ij}),
 \qquad s=2\Phi_G(h)-1,\quad k=4\phi(h)^2,
 \quad v=1-s^2-k.
\]

Keep A and -A unchanged. For any deterministic internal energy I(x,y),
put `M_I(H)=max_(x,y)|I(x,y)+x^T H y|`. If

\[
 Z_h=sB+\sqrt{k}\,G_R+\sqrt v\,W,
\]

where W has independent standard cross entries and is independent of
G_R, then an absolute constant C_* satisfies

\[
 \boxed{\displaystyle
 |\mathbb E M_I(\widehat B_h)-\mathbb E M_I(Z_h)|
                 \le C_* n^{16/11}.}                       \tag{22}
\]

Here is the masked proof and its provenance. Apply the arbitrary finite
observable version, including Section 6, of the
[shifted-sign theorem](NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md)
to the n^2 cross coordinates, with observable `sigma vec(xy^T)` and
prior proportional to `exp(beta sigma I(x,y))`. There are at most
`2^(2n+1)` augmented states and `||R||op<=3`. Equivalently, on all
unordered edges of K these observables have ZERO internal-edge
coordinates; the generic theorem permits this mask. With Y_h Gaussian
of mean sB and the EXACT centered covariance of the cross signing,
that theorem gives

\[
 |\mathbb E M_I(\widehat B_h)-\mathbb E M_I(Y_h)|
                       \le C_1 n^{16/11}.                 \tag{23}
\]

To replace the exact covariance, use Sections 3--5, in particular
equations (8) and (13), of the
[whole-edge remainder theorem](NOTE_2026-09-05_WHOLE_EDGE_SOURCE_PRESERVING_GAUSSIAN_REDUCTION.md)
on the complete source K of order N=2n. Its principal cross restriction
is exactly R. Principal restriction preserves the odd/even covariance
identities and does not increase operator norm. Thus, with
`b=vec(B)`, the exact cross covariance is

\[
 C_h=(kR+vI+wbb^T)+E_{\rm cross},\qquad
 0\le w\le {2\over D^2},\qquad
 \|E_{\rm cross}\|_{\rm op}\le{40\over\sqrt{2n}}.          \tag{24}
\]

Both C_h and `kR+vI+wbb^T` are PSD. The Gaussian finite-maximum
comparison with arbitrary deterministic offsets, applied to the
n^2-dimensional cross coefficient vectors of norm n and at most
`2^(2n+1)` states, shows that the operator error in (24) costs
`O(n^(5/4))` in the expected maximum. The retained rank-one Gaussian
is precisely `sqrt(w) zeta B`, with an independent scalar standard
Gaussian zeta. Its expected cross-norm cost is at most

\[
 \sqrt w\,\mathbb E|\zeta|\,
       \max_{x,y}|x^TBy|
 \le {2n^2\over\sqrt\pi D}=O(n),\qquad D\ge2n-2.           \tag{25}
\]

The same bound holds with every deterministic internal energy, by the
pointwise maximum Lipschitz inequality. Consequently
`|E M_I(Y_h)-E M_I(Z_h)|<=C_2 n^(5/4)`. Combining this with (23)
proves (22), uniformly in h and I. No full-order optimality of K and
no source operator cap are used.

Finally, if B is the actual conditional cross minimizer, every outcome
of the cross signing in (22) is an admissible competitor with the SAME
fixed internals. Taking `I=Q_A(x)-Q_A(y)` proves the floor

\[
 \boxed{\displaystyle
 F_A^*\le
 \mathbb E\Phi\!\begin{pmatrix}A&Z_h\\Z_h^T&-A\end{pmatrix}
                         +C_* n^{16/11}.}                \tag{26}
\]

This is the conditional-optimality inequality used with (20), not an
assumption that K is a full order-2n optimizer.
