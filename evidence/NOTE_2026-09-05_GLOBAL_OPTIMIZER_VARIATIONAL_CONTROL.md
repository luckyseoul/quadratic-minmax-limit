# Variational control of actual global optimizers

2026-09-05. **Proved structural results; original convergence OPEN.**
These arguments use complete signings, not a Paley or other prescribed
family. The norm results require global norm optimality where explicitly
stated. The pressure results concern pressure minima, which need not be
norm minima. The latter results often need only single-edge local
optimality. No comparison of normalized optima at different orders is
proved by this note.

## 1. Uniform sparse rounding near the sign domain

Write `E` for the unordered edges of `K_N`, `k=binom(N,2)`, and

\[
 Q_H(x)=\sum_{e=\{i,j\}}H_e x_i x_j,\quad
 \Phi(H)=\max_x|Q_H(x)|,\quad
 m_N=\min_{A_e\in\{-1,1\}}\Phi(A),\quad
 L=(N+1)\log2+1.
\]

For independent centered edge variables bounded in absolute value by two
and with total variance at most `V`, Bernstein's inequality and a union
bound over the `2^N` spin states give a realization with

\[
 \Phi(\xi)\le T(V):=\sqrt{2VL}+\frac43L.                 \tag{1}
\]

Indeed the fixed-state two-sided tail is at most
`2 exp(-t^2/[2(V+2t/3)])`. At the displayed threshold its exponent is
at least `L`, so the union bound is `2^(N+1) exp(-L)=exp(-1)<1`.
Independence between state energies is not required.

Fix `C>0`, any signing `A`, and any edge matrix `B` with `|B_e|<=C`.
Put `epsilon=C/sqrt(N)` and `H=A+B/sqrt(N)`. Round `H/(1+epsilon)`
independently to signs with its prescribed expectation. Relative to `A`
this flips edge `e` with probability

\[
 p_e=\frac{C-A_eB_e}{2(\sqrt N+C)}
 \in\left[0,\frac{C}{\sqrt N+C}\right].                    \tag{2}
\]

Thus the total centered variance is at most `4Ck/(sqrt(N)+C)`. Define

\[
 T_C=\sqrt{\frac{8CkL}{\sqrt N+C}}+\frac43L,\qquad
 E_C=(1+C/\sqrt N)T_C.
\]

Equation (1) supplies a complete signing `A'` such that
`Phi(A'-H/(1+epsilon))<=T_C`. Consequently, uniformly over all `A,B`,

\[
 \boxed{\Phi(A+B/\sqrt N)
 \ge(1+C/\sqrt N)(m_N-T_C)\ge m_N-E_C.}                    \tag{3}
\]

No norm cap on the underlying `A` was assumed. For `0<C<=sqrt(N)`,
`T_C=O(sqrt(C) N^(5/4)+N)`. In particular `C=o(sqrt(N))` gives
`E_C=o(N^(3/2))`.

This is not arbitrary dense fractional rounding: the per-edge variance
in (2) is `o(1)`. Coefficients in the interior of the cube, including
zero cross blocks, do not automatically have this property.

## 2. One common near-maximizing ensemble for a norm minimizer

Now assume `Phi(A)=m_N=M`. A signed state is `s=(sigma,x)`, where
`sigma in {-1,1}`. Put `q_s=sigma Q_A(x)` and `delta_s=M-q_s>=0`.
For a probability measure `mu` on signed states put
`r_e(mu)=E_mu[sigma x_i x_j]`.

Bilinear minimax on the compact box and finite probability simplex gives

\[
 \begin{aligned}
 \min_{|B_e|\le C}\Phi(A+B/\sqrt N)
 &=\min_B\max_\mu\left(\mathbb E_\mu q_s+
                 N^{-1/2}\sum_e B_e r_e(\mu)\right)\\
 &=\max_\mu\left(\mathbb E_\mu q_s-
                 CN^{-1/2}\sum_e|r_e(\mu)|\right).
 \end{aligned}                                               \tag{4}
\]

All extrema are attained. By (3), there is one common measure `mu_C`
satisfying

\[
 \boxed{\mathbb E_{\mu_C}\delta_s+
       \frac C{\sqrt N}\sum_e|r_e(\mu_C)|\le E_C.}          \tag{5}
\]

The measure does not depend on a subsequently chosen test matrix. For
fixed positive `C` this gives mean energy deficit `O_C(N^(5/4))` and
off-diagonal correlation `l1=O_C(N^(7/4))`. More generally,

\[
 C_N=o(\sqrt N),\qquad C_N\sqrt N\longrightarrow\infty       \tag{6}
\]

give deficit `o(N^(3/2))` and `l1=o(N^2)`, directly from the exact
formula for `E_C`.

Let `Gamma=E_mu[sigma xx^T]`. Its diagonal equals `E_mu sigma`, not
necessarily zero, and

\[
 \|\Gamma\|_F^2
 \le N+2\sum_e|r_e|
 \le N+\frac{2\sqrt N}{C}E_C=o(N^2)                         \tag{7}
\]

under (6). For independent samples from this same measure and any real
diagonal matrix `D` with `||D||op<=1`,

\[
 \left|\mathbb E[\sigma\sigma'(x^TDx')^2]\right|
 =|\operatorname{tr}(D\Gamma D\Gamma)|\le\|\Gamma\|_F^2.
                                                               \tag{8}
\]

This is a signed correlation matrix, not an unsigned positive-semidefinite
covariance bound. To obtain literal support in a near-maximizing cap, set
`eta_N=E_C/N^(3/2)` and condition on
`delta_s<=sqrt(eta_N)N^(3/2)`. Markov's inequality removes mass at most
`sqrt(eta_N)`; each bounded correlation changes by at most twice that
mass. The conditional measure therefore still has `l1=o(N^2)` under
(6). Nothing here requires exact maximizers to have this property.

## 3. Near-flat comparison for symmetric pressure

For a fixed `beta>0` define

\[
 F_\beta(H)=\log\mathbb E_x\cosh(\beta Q_H(x)),\qquad
 \mathcal P_N(\beta)=\min_A F_\beta(A).
\]

The pressure is convex and nonnegative, has `F_beta(0)=0`, and satisfies
`|F_beta(H)-F_beta(K)|<=beta Phi(H-K)`. The Section 1 rounding gives

\[
 \mathcal P_N(\beta)\le F_\beta(A')
 \le F_\beta(H/(1+\epsilon))+\beta T_C
 \le F_\beta(H)+\beta T_C.
\]

The last inequality is convexity toward zero and nonnegativity. Hence

\[
 \boxed{\mathcal P_N(\beta)-\beta T_C
 \le\min_{A,\ |B_e|\le C}F_\beta(A+B/\sqrt N)
 \le\mathcal P_N(\beta).}                                  \tag{9}
\]

For `beta=c/sqrt(N)`, fixed `c>0`, and `C=o(sqrt(N))`, the error is
`o(N)`. Again no norm cap on the underlying signing is needed.

## 4. Actual pressure-minimizing Gibbs measures are diffuse

Give each edge a positive strength `u_e` and let

\[
 H_A(x)=\sum_e u_e A_e x_ix_j,\qquad
 F(A)=\log\mathbb E_x\cosh H_A(x).
\]

Use the actual augmented Gibbs measure on `(sigma,x)` proportional to
`exp(sigma H_A(x))`, and write
`Gamma_e=<sigma x_i x_j>`, `r_e=A_e Gamma_e`.
Assume only that flipping any single sign of `A` cannot lower `F`.
Every global pressure minimizer satisfies this assumption.

The exact edge-flip identity gives

\[
 e^{F(A^e)-F(A)}=\cosh(2u_e)-r_e\sinh(2u_e)\ge1,
 \qquad r_e\le\tanh u_e.                                  \tag{10}
\]

For every signing, radial differentiation also gives

\[
 \sum_eu_e r_e
 =\frac{\mathbb E_x[H_A(x)\sinh H_A(x)]}
        {\mathbb E_x\cosh H_A(x)}\ge0.
\]

Combining this with (10) proves

\[
 \boxed{\sum_eu_e|\Gamma_e|
 \le2\sum_eu_e\tanh u_e\le2\sum_eu_e^2.}                  \tag{11}
\]

In particular, if `u_e=beta`,

\[
 \sum_e|\Gamma_e|\le2k\tanh\beta,\qquad
 \|\Gamma\|_F^2\le N+4k\tanh\beta.                        \tag{12}
\]

The matrix diagonal is `<sigma>`. At `beta=c/sqrt(N)`, (12) is
`O_c(N^(3/2))=o(N^2)`. More generally the last small-oh conclusion
requires `beta_N -> 0`. Formula (8) applies to this actual Gibbs measure;
no change of measure or zero-temperature identification is made.

### 4.1 Rowwise signed operator bound

At equal strengths `beta`, fix vertex `i` and write
`h_i=sum_(j!=i) A_ij x_j`, with `H_rest` containing the other edges.
Scaling only the incident strengths by `s>=0` and integrating `x_i`
gives the partition function
`E_rest[cosh(H_rest) cosh(s beta h_i)]`. Its derivative is nonnegative.
At `s=1` this implies
`sum_(j!=i) A_ij Gamma_ij>=0`, without an optimality assumption.
Combining with (10), at an edge-local pressure minimum,

\[
 \sum_{j\ne i}|\Gamma_{ij}|\le2(N-1)\tanh\beta,\qquad
 \|\Gamma\|_{\rm op}\le1+2(N-1)\tanh\beta.
\]

The second inequality follows from the maximum absolute row sum of the
real symmetric matrix, including its diagonal of magnitude at most one.
At fixed critical `c` this is `O_c(sqrt(N))`. It remains a bound on
the signed matrix, not on either one-sided phase covariance.

### 4.2 Full-row replacement and the cavity normalizer

For arbitrary positive strengths, now assume global pressure optimality
(or, more narrowly, optimality under replacement of this entire row).
Let `nu_i` be the probability measure on the other spins proportional to
`cosh(H_rest)`, and put `g_i=sum_(j!=i) u_ij A_ij x_j`. Factorization
after integrating `x_i` shows that the incident row minimizes
`E_(nu_i) cosh(g_i)` over all incident sign choices. Averaging independent
fair replacement signs gives, exactly,

\[
 \mathbb E_{\nu_i}\cosh(g_i)
 \le\prod_{j\ne i}\cosh(u_{ij}).
\]

This is `O_c(1)` when all `u_ij<=C_c/sqrt(N)`. However, the actual
full-Gibbs marginal on the remaining spins is `nu_i` tilted by
`cosh(g_i)`. The displayed normalizer bound is not a bound on
`E_(nu_i)[g_i^4 cosh(g_i)]`, nor on an unsigned star fluctuation in
that tilted measure. These stronger local facts do not supply the
defect comparison below.

### 4.3 Slowly cooled pressure minima are asymptotic norm minima

For uniform strengths, the finite signed state space has `2^(N+1)`
elements, so
`beta Phi(A)-(N+1)log2<=F_beta(A)<=beta Phi(A)`. If `A_N` globally
minimizes this pressure, comparison with a norm minimizer proves

\[
 0\le\Phi(A_N)-m_N\le\frac{(N+1)\log2}{\beta}.
\]

For its own Gibbs measure, relative entropy from the uniform signed-state
measure gives `beta <sigma Q_A>=F_beta(A)+D_KL>=F_beta(A)`. Therefore

\[
 0\le\Phi(A_N)-\langle\sigma Q_{A_N}\rangle
 \le\frac{(N+1)\log2}{\beta}.
\]

Taking `beta_N=c_N/sqrt(N)` with `c_N->infinity` and
`c_N=o(sqrt(N))` yields asymptotic norm optimality, a vanishing relative
mean deficit, and signed operator norm `o(N)` from Section 4.1. This
constructs near-optimal sources; it does not assert that every exact
norm minimizer is a pressure minimizer or control unsigned fluctuations.

## 5. Uniform control on the balanced two-block path

Let `N=2n`, with internal and cross groups `I,C`, and strengths

\[
 u_I(t)=c\sqrt{(2-t)/N},\qquad u_C(t)=c\sqrt{t/N},\quad
 0\le t\le1.
\]

There are `K_I=n(n-1)` internal edges and `K_C=n^2` cross edges. For
a fixed signing write `I_A=Q_(A_L)(x)+Q_(A_R)(y)` and `C_A=x^TBy`.
Pairing `y` with `-y` proves

\[
 \mathbb E\cosh(uI_A+vC_A)
 =\mathbb E[\cosh(uI_A)\cosh(vC_A)].                        \tag{13}
\]

Both partial derivatives are nonnegative: their numerators contain,
respectively, `I_A sinh(uI_A) cosh(vC_A)` and
`C_A sinh(vC_A) cosh(uI_A)`. Thus `sum_(e in g) r_e>=0` separately
for the internal group as a whole and for the cross group. This does
not assert separate positivity for the two individual diagonal blocks.
At an edge-local pressure minimum, (10) therefore gives

\[
 \sum_{e\in g}|\Gamma_e|\le2K_g\tanh u_g,
 \qquad
 \sum_e|\Gamma_e|
 \le\frac{cN^{3/2}}2(\sqrt{2-t}+\sqrt t)\le cN^{3/2}.       \tag{14}
\]

At `t=0` the cross correlations are zero by block-flip symmetry, so
the bound holds uniformly including that endpoint. In particular,

\[
 \left|\frac12\sum_e(u_e^2)'r_e^2\right|
 \le\frac{c^2}{2N}\sum_e|\Gamma_e|
 \le\frac{c^3\sqrt N}{2}.                                  \tag{15}
\]

## 6. Exact optimized-path identity, not a cross-order bound

Set `F_A(t)=log E cosh(u_I(t)I_A+u_C(t)C_A)` and
`f_N(t)=min_A F_A(t)`, over all complete signings. Each branch is
continuous on `[0,1]`, differentiable inside, and obeys
`|F_A'(t)|<=K_I|u_I'(t)|+K_C|u_C'(t)|`. This common bound is
integrable, including its `t^(-1/2)` endpoint behavior. The branches
and their finite minimum are consequently absolutely continuous.

Choose the lexicographically first active signing; it is Borel measurable.
Where `f_N` is differentiable, every active branch has the same derivative,
because `F_A-f_N` has a minimum zero there. Thus differentiating explicit
weights with the selected optimizer held fixed is valid almost everywhere.

Define

\[
 h(u)=\tanh u-u\operatorname{sech}^2u,
 \qquad 0\le h(u)\le\tfrac23u^3.
\]

The bound follows by integrating
`h'(u)=2u sech^2(u) tanh(u)` between zero and `u`. At the selected
optimizer define

\[
 d_e=u_e(1-r_e^2)+h(u_e)-r_e\ge0.                           \tag{16}
\]

As a function of `r in [-1,tanh(u)]` the expression is concave; its
endpoint values are `1+h(u)` and zero. This proves the inequality.
Define the signed group imbalance

\[
 \mathcal D_N(t)=
 \sum_{e\in C}\frac{c\,d_e}{2\sqrt N\sqrt t}
 -\sum_{e\in I}\frac{c\,d_e}{2\sqrt N\sqrt{2-t}}.           \tag{17}
\]

It is measurable and integrable. The almost-everywhere identity is

\[
 f_N'(t)=\frac12\sum_e(u_e^2)'(1-r_e^2)
              +\sum_eu_e'h(u_e)-\mathcal D_N(t).             \tag{18}
\]

Here `sum_e (u_e^2)'/2=c^2/4` exactly, since `K_C-K_I=N/2`.
The integral of the squared-correlation term is bounded by (15).
Also `sum_e u_e^2=c^2(N/2-1+t/2)<=c^2N/2` gives

\[
 \sum_e|u_e'|h(u_e)
 \le\frac{c^2}{3N}\sum_eu_e^2\le\frac{c^4}{6}.
\]

Consequently

\[
 \boxed{f_N(1)-f_N(0)=\frac{c^2}{4}
       -\int_0^1\mathcal D_N(t)\,dt+\mathcal E_N(c),\qquad
 |\mathcal E_N(c)|\le\frac{c^3\sqrt N}{2}+\frac{c^4}{6}.}    \tag{19}
\]

The error is `o(N)` for fixed `c`. Nonnegativity of each `d_e` does
not give a sign or lower-order bound for the imbalance integral.

## 7. Endpoint scope and provenance

Put `P_k(c)=mathcal P_k(c/sqrt(k))`. Then `f_N(1)=P_N(c)`, but

\[
 f_N(0)=\min_{A_L,A_R}\log\mathbb E_{x,y}
 \cosh\left(\frac c{\sqrt n}[Q_{A_L}(x)+Q_{A_R}(y)]\right)
                                                               \tag{20}
\]

is not automatically `2P_n(c)`. There is one shared sign, not independent
signs for the two blocks. More precisely, set

\[
 Z_\pm(A)=\mathbb E\exp(\pm cQ_A/\sqrt n),\quad
 a_A=\tfrac12(\log Z_+(A)+\log Z_-(A)),\quad
 b_A=\tfrac12(\log Z_+(A)-\log Z_-(A)).
\]

The block pressure equals `a_A+a_B+log cosh(b_A+b_B)`. Therefore

\[
 f_N(0)=2\min_A a_A\le2P_n(c).                             \tag{21}
\]

The lower bound follows from `log cosh>=0`; equality is attained by
an `a_A`-minimizer paired with its negative. The inequality follows from
`P_n(c)=min_A[a_A+log cosh(b_A)]`.

For the direction relevant to dyadic subadditivity, an estimate
`integral_0^1 D_N(t) dt >= -o(N)` would already give, by (19)--(21),
`P_(2n)(c)<=2P_n(c)+o(n)`. Endpoint equality or equivalence of the two
optimized pressures is not an extra gate for that implication. The
required signed integral estimate has not been proved. Moreover, a
bare dyadic inequality with an unspecified small-oh error is not itself
a proof of all-orders convergence. This optional pressure path imposes
no necessary conditions on other solutions of the original problem.
Convergence and its possible value remain open.

The complete sparse-rounding/minimax proof and complete weighted-pressure
path proof were read independently before integration. This note does
not reuse the earlier generic dense-rounding assertion: its changed
premise is vanishing per-edge variance in a shrinking neighborhood of
the sign domain. It also distinguishes an actual global norm minimizer,
an actual pressure minimizer, and an arbitrary low-norm signing.
Any accompanying order-four computation checks finite formulas only;
the all-orders statements above rest on the proofs, not that regression.
