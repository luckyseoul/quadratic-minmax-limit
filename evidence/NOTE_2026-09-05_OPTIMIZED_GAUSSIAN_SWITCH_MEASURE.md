# Optimized Gaussian smoothing: the switch measure and a trace bound

Fully and independently reviewed analytic proof, 2026-09-05. No computation is
used. The optimizer depends on the Gaussian field throughout this proof.
The result does not identify a deterministic profile derivative with a
Gaussian heat derivative and does not prove a cross-order inequality.

## 1. Setup and the negative switch measure

Let `N>=2`, let `E` be the unordered edges of `K_N`, and put
`k=binom(N,2)`. Fix a real edge profile `u=(u_e)` and `epsilon>=0`.
For every complete signing `A`, define, for `g in R^k`,

\[
 F_A(g)=\log\mathbb E_{\sigma,x}
 \exp\left(\sigma\sum_{e=\{i,j\}}u_e(A_e+\epsilon g_e)x_ix_j\right),
 \qquad \Psi(g)=\min_A F_A(g),                                \tag{1}
\]

where `(sigma,x)` is uniform on `{+1,-1}^{N+1}`. Each branch is real
analytic. At a point where the active branch is differentiable as a
minimum, write

\[
 \Gamma_e=\langle\sigma x_ix_j\rangle,
 \qquad \tau_e=\sigma x_ix_j
\]

for its actual Gibbs correlations. Then

\[
 \partial_e\Psi=\epsilon u_e\Gamma_e\quad\text{almost everywhere}.
                                                                    \tag{2}
\]

In the sense of matrix-valued distributional derivatives,

\[
 \boxed{D^2\Psi=S(g)\,dg-\mathsf M,\qquad
 S_{ef}(g)=\epsilon^2u_eu_f\operatorname{Cov}(\tau_e,\tau_f),}
                                                                    \tag{3}
\]

where `M` is a positive-semidefinite matrix-valued Radon measure. In
particular,

\[
 S_{ee}=\epsilon^2u_e^2(1-\Gamma_e^2).                       \tag{4}
\]

The minus sign in (3) is important. Along any line, a finite minimum of
analytic branches is smooth between crossings. At a crossing of active
branches its one-sided derivative can only jump downward. Its weak
second derivative therefore equals the selected smooth second derivative
minus the nonnegative jump masses. Identical branches along the line
have identical derivatives and cause no ambiguity. The multivariate
version gives the positive-semidefinite measure in (3).

For completeness, the smooth approximation below both justifies this
description and supplies the global Gaussian integrability needed later.
Let

\[
 \Psi_\kappa(g)=-\kappa\log\sum_A e^{-F_A(g)/\kappa},
 \qquad w_A(g)=\frac{e^{-F_A(g)/\kappa}}{\sum_B e^{-F_B(g)/\kappa}}.
\]

Then `|Psi_kappa-Psi|<=kappa log(2^k)` and

\[
 \begin{aligned}
 \partial_e\Psi_\kappa&=\epsilon u_e\overline\Gamma_e,
       &&\overline\Gamma_e=\sum_Aw_A\Gamma_{A,e},\\
 D^2\Psi_\kappa&=S_\kappa-M_\kappa,\\
 (S_\kappa)_{ef}&=\epsilon^2u_eu_f
                       \sum_Aw_A\operatorname{Cov}_A(\tau_e,\tau_f),\\
 M_\kappa&=\frac1\kappa
             \operatorname{Cov}_{w}(\nabla F_A)\succeq0.       \tag{5}
 \end{aligned}
\]

Distinct analytic branch functions tie only on a Lebesgue-null set; any
globally identical branches can first be merged. Thus `S_kappa` tends
almost everywhere to `S`, and the distributional limit of (5) gives
(3). Equivalently, a common global Hessian upper bound makes the finite
minimum semiconcave, so its singular Hessian part is negative
semidefinite. No independent-optimizer hypothesis appears.

## 2. Gaussian integration by parts and the explicit trace bound

Let `gamma_k` be the standard Gaussian density on `R^k`, and let
`G` have that density. Define

\[
 B_N(u,\epsilon):=
 \epsilon^2\|u\|_2^2+
 \epsilon\sqrt{2(N+1)\log2}\,\|u\|_2.                       \tag{6}
\]

Then

\[
 \boxed{
 0\le\int\operatorname{tr}(d\mathsf M(g))\,\gamma_k(g)
 \le B_N(u,\epsilon).
 }                                                                  \tag{7}
\]

To prove it, Gaussian integration by parts in (5) gives

\[
 \begin{aligned}
 \mathbb E\operatorname{tr}M_\kappa(G)
 &=\epsilon^2\sum_eu_e^2
             \mathbb E\left[1-\sum_Aw_A\Gamma_{A,e}^2\right]
       -\epsilon\,\mathbb E\sum_eu_eG_e\overline\Gamma_e\\
 &\le\epsilon^2\|u\|_2^2+
       \epsilon\,\mathbb E
       \max_{\sigma,x}\sigma\sum_eu_eG_ex_ix_j.               \tag{8}
 \end{aligned}
\]

The last inequality holds pointwise even with adaptive `A` and adaptive
weights `w_A`: the vector `overline(Gamma)` is the expectation of
`(sigma x_i x_j)_e` under a probability mixture of actual Gibbs measures.

Each Gaussian variable indexed by a signed state in (8) is centered and
has variance `||u||_2^2`. There are `2^(N+1)` signed states, possibly
with repetitions. For any `lambda>0`, the exponential-moment bound gives

\[
 \mathbb E\max_{\sigma,x}\sigma\sum_eu_eG_ex_ix_j
 \le\frac{(N+1)\log2}{\lambda}
       +\frac\lambda2\|u\|_2^2.
\]

Minimizing in `lambda` gives exactly the second term of (6).
Thus (8) is bounded uniformly in `kappa` by `B_N(u,epsilon)`.

Passing to the limit proves (7), and also the exact weak identity

\[
 \int\operatorname{tr}(d\mathsf M)\,\gamma_k
 =\epsilon^2\sum_eu_e^2\mathbb E(1-\Gamma_e^2)
      -\epsilon\mathbb E\sum_eu_eG_e\Gamma_e.                \tag{9}
\]

One can justify the noncompact integration by first inserting smooth
cutoffs into the Gaussian test function. Uniformly in `kappa`,
`||grad Psi_kappa||_2<=epsilon ||u||_2` and
`tr S_kappa<=epsilon^2 ||u||_2^2`; the Gaussian tails and their first
derivatives therefore make the cutoff errors tend to zero uniformly.
This also rules out loss of Gaussian-weighted switch mass at infinity.

In a critical profile with `||u||_2^2=O_c(N)`, equation (7) becomes

\[
 \int\operatorname{tr}(d\mathsf M)\,\gamma_k
 =O_c((\epsilon^2+\epsilon)N).                              \tag{10}
\]

Consequently the entire standard-Gaussian-weighted switch mass is `o(N)`
when `epsilon->0`. This uses no signed or unsigned covariance bound
beyond the elementary fact that a Gibbs correlation vector is an
expectation of signed-state vectors.

For any symmetric covariance-derivative matrix `W` with
`||W||op<=K`, positive semidefiniteness of the switch measure also gives
`|integral tr(W dM) gamma_k|<=K B_N(u,epsilon)`. Uniform boundedness
of this multiplying matrix is required to retain the stated smallness.

## 3. Exact one-edge cusp

For `N=2`, one edge of strength `u=beta>0`, and `epsilon>0`,

\[
 \Psi(g)=\min_{a\in\{-1,1\}}\log\cosh[\beta(a+\epsilon g)]
        =\log\cosh[\beta(1-\epsilon|g|)].                    \tag{11}
\]

Its derivative at zero jumps from `+beta epsilon tanh(beta)` on the
left to `-beta epsilon tanh(beta)` on the right. Hence

\[
 D^2\Psi=
 (\beta\epsilon)^2
 \operatorname{sech}^2[\beta(1-\epsilon|g|)]\,dg
 -2\beta\epsilon\tanh\beta\,\delta_0.                     \tag{12}
\]

The switch measure is exactly
`2 beta epsilon tanh(beta) delta_0`, with negative sign in the Hessian.
Its Gaussian-weighted mass is
`2 beta epsilon tanh(beta)/sqrt(2 pi)`.
Thus simply differentiating the selected smooth Gibbs branch and
discarding optimizer switches would already be wrong at order two.

## 4. What a genuine heat derivative does and does not give

Use `z` as the Gaussian argument in (1), and suppose independent
coordinate variances `v_e(t)>0` are differentiable. Wherever the heat
differentiation is justified (for example on a compact time interval
with variances bounded away from zero), its exact form is

\[
 \begin{aligned}
 \frac d{dt}\mathbb E_{G_{v(t)}}\Psi(G_{v(t)})
 &=\frac12\int\sum_e v'_e(t)\,\partial_{ee}\Psi\,
                                      \gamma_{v(t)}\\
 &=\frac{\epsilon^2}{2}\sum_e v'_e(t)u_e^2
                          \mathbb E(1-\Gamma_e^2)
   -\frac12\int\sum_ev'_e(t)\,d\mathsf M_{ee}\,
                                      \gamma_{v(t)}.          \tag{13}
 \end{aligned}
\]

If all `v'_e>=0`, the switch term has the favorable nonpositive sign.
In fact

\[
 \mathbb E_{G_v}\Psi(G_v)-\frac{\epsilon^2}{2}\sum_eu_e^2v_e
                                                                    \tag{14}
\]

is nonincreasing under coordinatewise variance increase. This follows
from (13): the derivative after subtracting the displayed quadratic
term is the sum of a nonpositive squared-correlation term and a
nonpositive switch term.

For mixed-sign variance derivatives the switch term has no automatic
sign. It can be bounded in absolute value by a suitable weighted trace
bound. For example, standardizing `G_v` and using the proof of Section 2
with perturbation profile `(u_e sqrt(v_e))_e` gives

\[
 \int\sum_e v_e\,d\mathsf M_{ee}\,\gamma_v
 \le\epsilon^2\sum_eu_e^2v_e
   +\epsilon\sqrt{2(N+1)\log2}\,
                        \sqrt{\sum_eu_e^2v_e}.               \tag{15}
\]

If `|v'_e|<=K v_e`, the absolute switch contribution in (13) is at most
`K/2` times the right side of (15). This assertion does not automatically
cover singular endpoints where that relative derivative bound fails.

Most importantly, varying the deterministic profile in
`sum_e u_e(t)(A_e+epsilon G_e)x_i x_j` has a first derivative containing
`sum_e u'_e A_e Gamma_e`. That is not the heat derivative (13).
Algebraically solving an integration-by-parts identity for a correlation
or for a normalized deterministic term can introduce `1/epsilon` or
`1/epsilon^2` factors. Such division can destroy the smallness in (10).
In particular, passing from the standard Gaussian argument `G` to the
physical additive perturbation `z=epsilon G` divides the corresponding
Gaussian-weighted Hessian switch bound by `epsilon^2`. Equation (10)
must not be quoted unchanged in those physical coordinates.

Thus this result retains and controls the actual optimizer-switch
measure. It does not justify equating a deterministic critical
two-block interpolation with monotone Gaussian heat, and does not
establish the remaining deterministic defect inequality or convergence
of the normalized complete-signing optimum.

## 5. A global-envelope correlation refinement and the actual derivative

At a positive strength `u_e`, coordinatewise Gaussian integration by parts
in (3), with the same cutoff justification as in Section 2, gives

\[
 \epsilon u_e\mathbb E[G_e\Gamma_e]
 =\epsilon^2u_e^2\mathbb E(1-\Gamma_e^2)-m_e,
 \qquad m_e:=\int d\mathsf M_{ee}\,\gamma_k\ge0.            \tag{16}
\]

Consider uniform strengths or one entire internal/cross group of the
balanced profile in `NOTE_2026-09-05_ADAPTIVE_PERTURBATION_CORRELATIONS.md`.
Physical group radial positivity and single-sign optimality give the
pointwise inequality
`L_g<=2K_g tanh(u_g)+epsilon sum_(e in g) G_e Gamma_e`.
Divide (16) by `u_e`, sum within the group, and retain the nonpositive
switch contribution. This proves

\[
 \boxed{\mathbb E L_g
 \le2K_g\tanh u_g+\epsilon^2u_gK_g.}                       \tag{17}
\]

At zero cross strength the cross correlations vanish separately, so
the same conclusion holds without division by zero. In particular, for
uniform strength `beta`,

\[
 \mathbb E\|\Gamma\|_F^2
 \le N+4k\tanh\beta+2\epsilon^2\beta k.                    \tag{18}
\]

The condition `(1+epsilon_N^2)beta_N->0` suffices for signed diffuseness
in mean and probability. This calculation uses the actual GLOBAL
Gaussian envelope, not an arbitrary edge-local selection. The alternative
Boolean-energy estimate in Section 5 of the adaptive-perturbation note
has weaker optimality assumptions and different noise-scale dependence.
Neither estimate asserts a corresponding expected maximum-row bound.

For clarity, apply (16) to the derivative when the deterministic profile
itself varies. Let each `u_e(t)>0` be continuously differentiable on
compact interior intervals, as holds for the balanced path, and set

\[
 f_\epsilon(t)=\mathbb E_G\min_A
 \log\mathbb E_{\sigma,x}
 \exp\!\left(\sum_eu_e(t)(A_e+\epsilon G_e)\sigma x_ix_j\right).
\]

Finite-branch absolute continuity and envelope differentiation, followed
by (16) at each profile, give almost everywhere

\[
 \boxed{
 f_\epsilon'(t)
 =\sum_eu_e'\mathbb E[A_e\Gamma_e]
   +\epsilon^2\sum_eu_eu_e'\mathbb E(1-\Gamma_e^2)
   -\sum_e\frac{u_e'}{u_e}m_e(t).}                         \tag{19}
\]

For fixed `N`, interchange with the Gaussian expectation is justified
by the integrable bound `sum_e |u_e'|(1+epsilon |G_e|)`.
Here the switch measure is always defined in standard Gaussian argument
coordinates, as in (3). Signed diffuseness controls the squared-correlation
part at bounded noise; it does not eliminate the first term of (19).
The last term has mixed weights on the balanced path. A monotone variance
path gives a favorable switch sign, but leaves the deterministic mean
term whenever the mean profile is changed. Covariance transport at fixed
mean cannot supply zero cross means at one endpoint and nonzero complete
signing cross means at the other. This remains a cross-order gap.

## 6. Sign-flip generators recover the already isolated defect

This section is noiseless. For an edge of strength `u>0` and actual
correlation `r=A_e Gamma_e`, its exact log pressure increase on flipping
the sign is

\[
 \phi_u(r)=\log[\cosh(2u)-r\sinh(2u)].                     \tag{20}
\]

At an edge-local minimum it is nonnegative. Recall the defect in Section 6
of `NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md`:
`d_u(r)=u(1-r^2)+h(u)-r`, where
`h(u)=tanh u-u sech^2 u` and `0<=h(u)<=2u^3/3`.
For every `u>0` and `r in [-1,1]`,

\[
 \left|\frac{\phi_u(r)}{2u}-d_u(r)\right|
 \le\frac43u^2|r|+\frac43u^3.                             \tag{21}
\]

To verify this, put `b_r(s)=log(cosh s+r sinh s)`.
Its first three derivatives at zero are `r`, `1-r^2`, and
`-2r(1-r^2)`. Writing `m=b_r'`, one has
`b_r''''=-2(1-m^2)(1-3m^2)` and `|b_r''''|<=2`.
Taylor expansion at `s=-2u` gives
`phi=-2ur+2u^2(1-r^2)+(8/3)u^3r(1-r^2)+R`,
with `|R|<=4u^4/3`. Divide by `2u` and subtract `d_u(r)` to obtain (21).
The formulas also hold at `r=+1,-1` by direct evaluation or continuity.

At the selected edge-local minimum on the noiseless balanced profile,
define, for `0<t<=1`,

\[
 \mathcal J_N(t)=\sum_e\frac{u_e'}{2u_e}\phi_{u_e}(r_e)
 =\frac1{4t}\sum_{e\in C}\phi_{u_C}(r_e)
  -\frac1{4(2-t)}\sum_{e\in I}\phi_{u_I}(r_e).             \tag{22}
\]

Since `|u_e'|u_e=c^2/(2N)`, `sum_e u_e|r_e|<=2 sum_e u_e^2`, and
`sum_e u_e^2<=c^2N/2`, (21) proves

\[
 |\mathcal J_N(t)-\mathcal D_N(t)|
 \le\frac43\frac{c^2}{2N}
       \left(\sum_eu_e|r_e|+\sum_eu_e^2\right)
 \le c^4.                                                 \tag{23}
\]

The already proved integrability of `D_N` and (23) also proves
integrability of `J_N`; their integrals differ by at most `c^4`.
Thus the needed defect estimate can equivalently, up to this bounded
error at fixed `c`, be stated for a signed imbalance of actual edge-flip
pressure gaps. Nonnegativity of each gap still does not sign that imbalance.

For a fixed branch, sign-flip rate `lambda` and amplitude drift `u'`
have pressure generator `u'r+lambda phi_u(r)`. Making the physical
coefficient `uA` a martingale requires `lambda=u'/(2u)`, which is a
nonnegative rate only when `u'>=0`. Under optimization over all signings,
any multiplicative sign-noise vector is absorbed exactly by relabeling
the admissible signing set. The optimized jump contribution is therefore
zero; reoptimization cancels the fixed-branch jump gain. Equation (23)
identifies the resulting signed reset imbalance with the existing defect,
not a new error and not a heat-based proof of its sign. Decreasing internal
amplitudes cannot be represented by negative stochastic jump rates.

The refinements and exact identities in Sections 5--6 were independently
reviewed. No new finite computation, all-orders comparison, or convergence
claim accompanies them.
