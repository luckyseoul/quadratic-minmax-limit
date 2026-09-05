# A weighted fourth-moment budget along the actual row tilt

2026-09-05. **All-orders actual-law estimate; no transport conclusion.**
This is an integral in the amplitude of one fixed row, not an integral
over independently reoptimized temperatures. No computation is used.

## 1. The actual deleted law and row interpolation

Fix c,lambda>0 and an ACTUAL global minimizer of the balanced quartic
profile. Use the definitions and normalization of
`NOTE_2026-09-05_QUARTIC_PROFILE_ROW_RESET.md`. For a fixed vertex write

\[
 M=\begin{pmatrix}0&w^T\\w&B\end{pmatrix},\qquad
 d=\|w\|_2^2\le1,\qquad h(x_T)=c w^Tx_T,
 \quad H=|h|.
\]

Let mu_B be the ACTUAL augmented symmetric Gibbs law of B on
(sigma,x_T). It includes the original temperature and normalization;
it is not a smaller-order minimizing law. Let nu be mu_B times an
independent uniform removed spin x_i. For 0<=alpha<=1 put

\[
 Z_\alpha=\mathbb E_{\mu_B}\cosh(\alpha h),\qquad
 f(\alpha)=\log Z_\alpha,\qquad
 \frac{d\nu_\alpha}{d\nu}
     =\frac{\exp(\alpha\sigma x_i h)}{Z_\alpha}.
                                                               \tag{1}
\]

Thus nu_1 is exactly the full Gibbs law of M. At other alpha this is
the actual Gibbs law of the same matrix with just this row and column
multiplied by alpha. These intermediate matrices need not be sign
optimizers. For any function of H, the remaining-spin marginal gives

\[
 \mathbb E_{\nu_\alpha}\psi(H)
 =\frac{\mathbb E_{\mu_B}[\psi(H)\cosh(\alpha H)]}{Z_\alpha}.
                                                               \tag{2}
\]

The independently proved whole-row and edge-optimality budgets give

\[
 0\le p_i:=f(1)\le c^2/2+4\lambda,
 \qquad
 0\le\mathcal E_i:=f'(1)
      =\mathbb E_{\nu_1}[H\tanh H]\le c^2+8\lambda.       \tag{3}
\]

These inequalities hold separately for every row of every active
optimizer at every profile parameter, including zero cross weights.

## 2. Exact weighted identity and its consequence

For every H>=0, with the value at H=0 interpreted by continuity,

\[
 \int_0^1(1-\alpha)^2\cosh(\alpha H)\,d\alpha
       =\frac{2(\sinh H-H)}{H^3}.                         \tag{4}
\]

All terms below are nonnegative. Using (2), Tonelli, and (4) yields
the EXACT unnormalized identity

\[
 \begin{aligned}
 &\int_0^1(1-\alpha)^2Z_\alpha
               \mathbb E_{\nu_\alpha}H^4\,d\alpha\\
 &\quad=2\mathbb E_{\mu_B}[H(\sinh H-H)]\\
 &\quad=2Z_1\left(\mathcal E_i-
                      \mathbb E_{\nu_1}[H^2\operatorname{sech}H]
              \right).
 \end{aligned}                                              \tag{5}
\]

Since Z_alpha>=1, equations (3)--(5) prove the normalized bound

\[
 \boxed{\quad
 \int_0^1(1-\alpha)^2\mathbb E_{\nu_\alpha}|h|^4\,d\alpha
 \le2e^{p_i}\mathcal E_i
 \le2e^{c^2/2+4\lambda}(c^2+8\lambda).
 \quad}                                                       \tag{6}
\]

The same fixed cavity law occurs throughout (4)--(6). In particular,
no full-law exponential tail was inferred from a cavity-law tail.
Summing (6) over rows gives O_(c,lambda)(N), not o(N).

## 3. What the budget does and does not bound

Formula (6) is a weighted interpolation estimate. It is not an
unweighted fourth moment at alpha=1 or a uniform-integrability claim
for the endpoint actual laws. Over the full interval it contains no
N^(-1/2) factor. A fourth-order Taylor estimate for a full row based
only on (6) therefore does not supply an O(N^(-1/2)) per-row error.

There is a different, coordinatewise reason that a biased row refill
has a small conditional averaging error. If xi_j are independent signs
with mean alpha and the refilled row is w'_j=xi_j w_j, then for each
fixed remaining spin configuration and each s in {-1,1},

\[
 \begin{aligned}
 \log\mathbb E_\xi\exp(sc\sum_j\xi_jw_jx_j)
 &=\sum_j\log[\cosh(cw_jx_j)+s\alpha\sinh(cw_jx_j)]\\
 &=s\alpha h+\tfrac12(1-\alpha^2)c^2d+\varepsilon_s(x_T),
 \end{aligned}                                               \tag{7}
\]

where, uniformly in alpha, s, and x_T,

\[
 |\varepsilon_s(x_T)|
 \le\tfrac43c^3\sum_j|w_j|^3
 \le\tfrac{4\sqrt2}{3}c^3d\,N^{-1/2}.                    \tag{8}
\]

Indeed, the third derivative of the scalar log moment-generating
function of a sign is its third centered moment, whose absolute value
is at most 8; Taylor's theorem gives (8). The profile bound is
`max_j |w_j|<=sqrt(2/N)`.

After summing both removed-spin values and integrating against mu_B,
(7)--(8) give

\[
 \left|\log\mathbb E_\xi e^{F(M')}-
       \left[F(B)+f(\alpha)+\tfrac12(1-\alpha^2)c^2d\right]
 \right|
 \le\tfrac{4\sqrt2}{3}c^3d\,N^{-1/2}.                    \tag{9}
\]

This bound is uniform for every fixed host, so it may be applied to
whichever optimizer is active. It does not differentiate an optimizer
selector or identify log E exp(F) with E F. In a global-optimality
comparison, Jensen supplies only the direction
`E_xi F(M')<=log E_xi exp(F(M'))`. The annealed/quenched gap has not
been bounded by (9). The weighted fourth-moment budget (6) is not
needed to prove (9).

The remaining balanced-profile signed gap integral is unchanged.
