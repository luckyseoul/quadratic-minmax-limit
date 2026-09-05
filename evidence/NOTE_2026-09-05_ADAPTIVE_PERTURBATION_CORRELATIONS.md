# Adaptive perturbations preserve signed correlation bounds

2026-09-05. Analytic extension, fully and independently reviewed; no computation is used. The statements
below concern actual edge-local pressure minima with the perturbation held
fixed during an edge flip. They do not prove a cross-order pressure
comparison or control optimizer-switching terms in a heat equation.

## 1. Uniform strengths: pointwise bounds for every perturbation

Let `N>=2`, index `e` by the unordered edges of `K_N`, and put
`k=binom(N,2)`. Let `E=(E_e)` be a real edge perturbation, with
`||E||_2^2=sum_e E_e^2`. Fix `beta>0` and set

\[
 H_{A,E}(x)=\beta\sum_{e=\{i,j\}}(A_e+E_e)x_ix_j,
 \qquad F_\beta(A+E)=\log\mathbb E_x\cosh H_{A,E}(x).
\]

The variable `A_e` is a sign, while the physical coefficient `A_e+E_e`
can be arbitrary. Use the augmented Gibbs measure on `(sigma,x)`
proportional to `exp(sigma H_(A,E)(x))`, and define

\[
 \Gamma_{ij}=\langle\sigma x_ix_j\rangle,\qquad
 r_e=A_e\Gamma_e,\qquad
 L=\sum_e|\Gamma_e|.
\]

Here `Gamma_ii=<sigma>`. Assume that replacing any one `A_e` by `-A_e`,
with `E` fixed, cannot decrease `F_beta(A+E)`. In particular, every global
minimizer over signings has this property, even when it was selected after
observing `E`.

The physical coefficient changes by `-2A_e`, not by `-2(A_e+E_e)`.
Consequently the exact edge-flip identity is

\[
 e^{F_\beta(A^e+E)-F_\beta(A+E)}
 =\cosh(2\beta)-r_e\sinh(2\beta)\ge1,
 \qquad r_e\le\tanh\beta.                                  \tag{1}
\]

Radially scaling the entire physical Hamiltonian gives, for every signing,

\[
 \sum_e(A_e+E_e)\Gamma_e
 =\frac{\mathbb E_x[Q_{A+E}(x)\sinh(\beta Q_{A+E}(x))]}
        {\mathbb E_x\cosh(\beta Q_{A+E}(x))}\ge0.
                                                               \tag{2}
\]

Thus `sum_e r_e>=-sum_e E_e Gamma_e`. Since `|r_e|=|Gamma_e|`,

\[
 \begin{aligned}
 L&=2\sum_e(r_e)_+-\sum_e r_e\\
  &\le2k\tanh\beta+\sum_e E_e\Gamma_e\\
  &\le2k\tanh\beta+\|E\|_2\sqrt L.
 \end{aligned}                                                \tag{3}
\]

The last step uses Cauchy--Schwarz and `sum_e Gamma_e^2<=L`.
Applying `ab<=(a^2+b^2)/2` to the last term proves

\[
 \boxed{L\le4k\tanh\beta+\|E\|_2^2,\qquad
 \|\Gamma\|_F^2\le N+8k\tanh\beta+2\|E\|_2^2.}             \tag{4}
\]

The second inequality includes the diagonal and both copies of each
off-diagonal entry. Hence, along any deterministic sequence with
`beta_N->0` and `||E_N||_2=o(N)`, the actual selected Gibbs matrix has
`||Gamma||_F^2=o(N^2)`. No entrywise bound on `E`, probabilistic
independence of `A` from `E`, or norm-optimality premise is needed.

## 2. Uniform strengths: rowwise operator bound

Fix vertex `i` and put
`h_i=sum_(j!=i)(A_ij+E_ij)x_j`. Scale all physical coefficients incident
to `i` by `s>=0`, leaving the remaining Hamiltonian `H_rest` unchanged.
Integration over `x_i` gives

\[
 Z_i(s)=\mathbb E_{x_{-i}}
       [\cosh(H_{\rm rest})\cosh(s\beta h_i)].
\]

Its derivative is nonnegative. At `s=1` this implies
`sum_(j!=i)(A_ij+E_ij)Gamma_ij>=0`. Apply the proof of (3) to this row,
using (1) on each incident edge. With `L_i=sum_(j!=i)|Gamma_ij|`,

\[
 \boxed{L_i\le4(N-1)\tanh\beta+\sum_{j\ne i}E_{ij}^2,\qquad
 \|\Gamma\|_{\rm op}
 \le1+4(N-1)\tanh\beta+\max_i\sum_{j\ne i}E_{ij}^2.}        \tag{5}
\]

The operator bound is the maximum absolute row-sum bound for a real
symmetric matrix; the diagonal contributes at most one. It gives
`||Gamma||op=o(N)` when `beta_N->0` and the maximum perturbation row
squared norm is `o(N)`. The weaker total condition in (4) alone does
not imply this row condition. These remain signed-matrix estimates.

## 3. Balanced two-block weights: pointwise group bounds

Let `N=2n`, with fixed left and right vertex blocks. Let `I,C` denote
the internal and cross edge groups, with
`K_I=n(n-1)` and `K_C=n^2`. For fixed `c>0` and `0<=t<=1`, put

\[
 u_I(t)=c\sqrt{(2-t)/N},\qquad u_C(t)=c\sqrt{t/N}.
\]

Keep the perturbation `E` fixed, and write the physical forms as

\[
 I_{A,E}=Q_{A_L+E_L}(x)+Q_{A_R+E_R}(y),\qquad
 C_{A,E}=x^T(B+E_C)y.
\]

The Hamiltonian is `u_I I_(A,E)+u_C C_(A,E)`. At each `t`, choose any
edge-local pressure minimum over `A`, with this `E` and these weights
held fixed. No consistency of the selected signing between different
values of `t` is required for the bounds below.

Pairing `y` with `-y` proves

\[
 \mathbb E\cosh(u I_{A,E}+v C_{A,E})
 =\mathbb E[\cosh(u I_{A,E})\cosh(v C_{A,E})].
                                                               \tag{6}
\]

Both partial derivatives are nonnegative for `u,v>=0`. Therefore
`sum_(e in g)(A_e+E_e)Gamma_e>=0` separately for the entire internal
group and for the cross group. This does not assert separate positivity
for the two individual diagonal blocks.

For positive `u_g`, the edge-flip identity gives `r_e<=tanh u_g`.
Repeating (3) within each group yields

\[
 \boxed{L_g:=\sum_{e\in g}|\Gamma_e|
 \le4K_g\tanh u_g+\sum_{e\in g}E_e^2.}                     \tag{7}
\]

At `t=0`, all cross correlations are zero by block-flip symmetry, so
(7) remains valid there; the zero-strength case is not deduced by
dividing the edge-flip identity by `sinh(2u_C)`.

Using `tanh u<=u` and `sqrt(2-t)+sqrt(t)<=2` gives, uniformly in `t`,

\[
 L\le2cN^{3/2}+\|E\|_2^2,\qquad
 \|\Gamma\|_F^2\le N+4cN^{3/2}+2\|E\|_2^2.                 \tag{8}
\]

In particular, the squared-correlation term in the optimized-path
calculation satisfies

\[
 \boxed{\left|\frac12\sum_e(u_e^2)'r_e^2\right|
 \le\frac{c^2}{2N}L
 \le c^3\sqrt N+\frac{c^2}{2N}\|E\|_2^2.}                 \tag{9}
\]

For fixed `c` and `||E||_2=o(N)`, this is `o(N)`, uniformly including
the endpoints. This is only one term of a pressure comparison. For a
fixed active branch and `0<t<1`, the actual weight derivative also contains

\[
 \sum_e u_e' E_e\Gamma_e,
 \qquad
 \frac{dF_{A,E}}{dt}=\sum_e u_e'(r_e+E_e\Gamma_e).           \tag{10}
\]

Equations (7)--(9) do not assign a favorable sign to (10), to the
internal-versus-cross defect imbalance, or to an optimizer-switching
contribution in a Gaussian calculation.

## 4. Gaussian consequences with adaptive sign choices

For each `N`, let the unordered entries `G_e` be independent standard
real Gaussians, and put `E_e=epsilon_N G_e`, with deterministic
`epsilon_N>=0`. Symmetrize with zero diagonal. The signing may be any
measurable edge-local minimizer selected after observing all of `G`.
The deterministic inequalities apply pointwise, so no independence
between that signing, its Gibbs correlations, and `G` is used below.

First, `E||E||_2^2=epsilon_N^2 k`. In the uniform-strength model,

\[
 \mathbb E\|\Gamma\|_F^2
 \le N+8k\tanh\beta_N+2\epsilon_N^2k.                       \tag{11}
\]

If `beta_N->0` and `epsilon_N->0`, then
`||Gamma||_F^2/N^2->0` in mean and hence in probability. Also
`||E||_2/N->0` in mean square, so the perturbation hypothesis of (4)
holds in probability. No almost-sure assertion across different `N`
is needed or implied here.

For the operator norm, put `d=N-1` and
`Y_i=sum_(j!=i)G_ij^2`. Each `Y_i` has Gaussian-square exponential
moment `E exp(Y_i/4)=2^(d/2)`. The rows need not be independent.
The log-sum-exp bound and Jensen's inequality give

\[
 \mathbb E\max_iY_i
 \le4\log N+2d\log2.
\]

Alternatively, a union bound gives, for every `0<delta<1`,

\[
 \Pr\{\max_iY_i>2d\log2+4\log(N/\delta)\}\le\delta.
\]

Thus (5) implies the explicit expectation estimate

\[
 \mathbb E\|\Gamma\|_{\rm op}
 \le1+4d\tanh\beta_N+
       \epsilon_N^2(2d\log2+4\log N).                       \tag{12}
\]

Under the same assumptions `beta_N->0`, `epsilon_N->0`, it follows
that `||Gamma||op/N->0` in mean and in probability. At fixed critical
`beta=c/sqrt(N)`, (12) is `O_c(sqrt(N)+epsilon_N^2 N)`; retaining an
`O(sqrt(N))` bound requires, for example, `epsilon_N=O(N^(-1/4))`.
Merely assuming `epsilon_N->0` guarantees `o(N)`, not that sharper rate.

For the balanced profile at fixed `c`, (7) gives
`E L_g<=4K_g tanh(u_g)+epsilon_N^2 K_g`. The same deterministic random
upper bounds apply simultaneously at every `t` and to every chosen
edge-local minimizer at that `t`. In particular,

\[
 \mathbb E\|\Gamma(t)\|_F^2
 \le N+4cN^{3/2}+2\epsilon_N^2k,
\]

and (9) has expected upper bound

\[
 c^3\sqrt N+\frac{c^2\epsilon_N^2k}{2N}.                    \tag{13}
\]

If `epsilon_N->0`, these yield respectively signed diffuseness and an
`o(N)` squared-correlation error, uniformly in the profile parameter.
For a growing `c_N`, retain the explicit bounds: sufficient conditions
for (13) to be `o(N)` are `c_N^3/sqrt(N)->0` and
`c_N^2 epsilon_N^2->0`. Those stronger conditions should not be confused
with the weaker `c_N=o(sqrt(N))` condition sufficient for diffuseness.

These corollaries repair signed correlation control for actual
noise-adaptive minima. They do not establish the still-open deterministic
profile-amplitude comparison, nor any bound on unsigned phase covariance.

## 5. Energy-norm refinement: bounded Gaussian noise is also diffuse

The estimate before Cauchy--Schwarz in (3) gives a different pointwise
bound. For any real edge array `B`, put
`Phi(B)=max_x |sum_e B_e x_i x_j|`. The actual Gibbs expectation satisfies
`sum_e E_e Gamma_e=<sigma Q_E><=Phi(E)`, even when the signing was
selected after observing `E`. Therefore every edge-local minimum obeys

\[
 L\le2k\tanh\beta+\Phi(E),\qquad
 \|\Gamma\|_F^2\le N+4k\tanh\beta+2\Phi(E).                \tag{14}
\]

This can be combined with (4) by taking the smaller upper bound. In
particular, the deterministic hypothesis `Phi(E_N)=o(N^2)` and
`beta_N->0` suffices for signed Frobenius diffuseness; a small Euclidean
norm relative to `N` is not required for this alternative.

The same argument applies separately to the two entire groups in (6).
Embed `E_g` in the complete edge array by setting all other edges to zero.
Then, uniformly on the balanced path, including zero cross strength,

\[
 L_g\le2K_g\tanh u_g+\Phi(E_g),\qquad
 L\le cN^{3/2}+\Phi(E_I)+\Phi(E_C).                         \tag{15}
\]

For `E=epsilon G` with independent standard Gaussian edges, each signed
state's value of `Q_(E_g)` is centered Gaussian of variance
`epsilon^2 K_g`. The exponential-moment bound for the maximum of the
`2^(N+1)` signed-state values proves

\[
 \mathbb E\Phi(E_g)
 \le\epsilon\sqrt{2(N+1)\log2}\sqrt{K_g}.                  \tag{16}
\]

No independence of these state values is needed. For uniform strengths,
(14) consequently gives

\[
 \mathbb E\|\Gamma\|_F^2
 \le N+4k\tanh\beta+
        2\epsilon\sqrt{2(N+1)\log2}\sqrt k.                \tag{17}
\]

Thus `beta_N->0` and `epsilon_N=o(sqrt(N))` imply signed diffuseness
in mean and probability for every measurable edge-local selection.
In particular, fixed bounded noise is permitted. This strengthens the
vanishing-noise sufficient conditions in Section 4, without changing
the separately stated row-operator estimate there.

For the balanced groups, `sqrt(K_I)+sqrt(K_C)<=N`. Equations (15)--(16)
and `|(u_e^2)'|=c^2/N` therefore prove, uniformly in `0<=t<=1`,

\[
 \mathbb E\left|\frac12\sum_e(u_e^2)'r_e^2\right|
 \le\frac{c^3\sqrt N}{2}
      +\frac{c^2\epsilon}{2}\sqrt{2(N+1)\log2}.             \tag{18}
\]

For fixed `c` and bounded `epsilon`, this is `O_(c,epsilon)(sqrt(N))`,
hence `o(N)`. For growing parameters, a sufficient condition is
`(c_N^3+c_N^2 epsilon_N)/sqrt(N)->0`. All statements in this section
use only edge-local optimality, not a global Gaussian envelope. The
separate global-envelope expectation refinement is proved in Section 5
of `NOTE_2026-09-05_OPTIMIZED_GAUSSIAN_SWITCH_MEASURE.md`.

The argument bounds the perturbation's signed Gibbs energy by its
Boolean norm. It does not bound either unsigned phase covariance,
improve the maximum-row estimate at bounded noise, or remove the
deterministic derivative term (10). This full refinement was also
independently reviewed; no finite computation is used.
