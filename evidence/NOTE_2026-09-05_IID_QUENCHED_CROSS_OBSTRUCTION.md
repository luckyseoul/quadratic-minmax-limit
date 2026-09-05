# Quenched iid cross blocks and the information cost of selection

2026-09-05. **All-orders, method-scoped theorem; convergence OPEN.**

This note concerns the average of the actual log pressure, not an
annealed partition function or a quadratic moment-generating-function
proxy. At sufficiently large fixed temperature and a fixed full balanced
step, independent cross signs have a linear excess over the optimized
paired endpoint. The conclusion is uniform over the internal host,
including hosts selected after the cross block is seen.

The proof also excludes proposal distributions whose relative entropy
from independent fair cross signs is `o(n)`, and `exp(o(n))` samples
with that iid matrix marginal. It does not exclude arbitrary dependent
rounding, a selected rare cross block, the actual Gram--Schmidt law, or
the earlier correlated Gaussian-sign law.

No finite signing census, Gaussian simulation, numerical value of the
SK constant, or numerical integration is used.

## 1. Precise objects and conclusion

Retain the notation of the finite-step rounding note. For a complete
symmetric zero-diagonal signing `A` of order `n`, put

\[
 Q_A(x)=\sum_{i<j}A_{ij}x_ix_j,\qquad
 a_A(s)=\tfrac12[\log\mathbb E e^{sQ_A}+\log\mathbb E e^{-sQ_A}],
 \qquad R_n(s)=\min_A a_A(s).
\]

Fix `c>0` and `0<t<=1` independently of `n`, and define

\[
 \beta=c/\sqrt n,\qquad
 \eta=\beta\sqrt{1-t/2},\qquad
 \gamma=\beta\sqrt{t/2},\qquad m=n^2.
\]

For `B in {−1,1}^{n by n}`, including independent choices on its
diagonal, define the actual paired log-cosh pressure

\[
 F_{A,B}(t)=\log\mathbb E_{x,y}
 \cosh[\eta(Q_A(x)-Q_A(y))+\gamma x^TBy].                 \tag{1}
\]

Let `P_0` be the law of all `n^2` independent fair signs in `B`.
The constant used below is explicit:

\[
 K_0=\frac4{3\sqrt\pi}>\frac1{\sqrt2}.                   \tag{2}
\]

**Theorem.** Uniformly over every complete host `A`,

\[
 \mathbb E_{P_0}F_{A,B}(t)
 \ge [c\sqrt{2t}\,K_0-2\log2]n-o_{c,t}(n).              \tag{3}
\]

In fact (3) remains valid with `F_{A,B}` replaced by `min_A F_{A,B}`.
Since the established all-orders construction bound implies

\[
 2R_n(\beta)\le cn+o_c(n),
\]

putting

\[
 \Delta(c,t)=c(\sqrt{2t}\,K_0-1)-2\log2,                 \tag{4}
\]

gives a positive linear excess whenever `Delta(c,t)>0`. Such fixed
parameters exist at `t=1`, by (2). In particular, this refutes the
candidate iid-quenched estimate

\[
 \mathbb E_{P_0}F_{A,B}(1)\le2a_A(\beta)+o(n)
\]

for actual `a_A(beta)` minimizers at all sufficiently large fixed `c`.

More generally, for any law `D_n` on cross signs,

\[
 \mathbb E_{D_n}\min_A F_{A,B}(t)
 \ge [c\sqrt{2t}\,K_0-2\log2]n
 -c\sqrt{tn\,D(D_n\Vert P_0)}-o_{c,t}(n).                \tag{5}
\]

Thus relative entropy `D(D_n||P_0)=o(n)` does not remove the excess.
The proposal law may depend arbitrarily on the order, parameters,
and a specified host. The bound is uniform in those choices.

## 2. Primary input and normalization

The only spin-glass input is the rigorous zero-temperature Parisi
formula of Auffinger and Chen, *Parisi formula for the ground state
energy in the mixed p-spin model*, Annals of Probability 45 (2017),
4617--4631. The exact version inspected here is
[arXiv:1606.05335v2](https://arxiv.org/pdf/1606.05335v2).
Its PDF page 1 sets the SK covariance function to `xi(q)=q^2/2`.
Theorem 1, equation (5), is on PDF page 3. The control representation
is Theorem 2 on PDF pages 5--6 and Corollary 2, equation (12), on
PDF page 7; the following paragraph extends it to all admissible
integrable parameters.

In the normalization

\[
 H_N(z)=\frac1{\sqrt N}\sum_{i<j}g_{ij}z_i z_j,
 \qquad g_{ij}\text{ iid }N(0,1),
\]

write

\[
 P_{\rm SK}=\lim_{N\to\infty}\frac1N\mathbb E\max_z H_N(z).
\]

The source uses covariance `N xi(q)=Nq^2/2`; adding a common
independent `N(0,1/2)` variable to `H_N` gives that covariance exactly.
Such a common centered shift changes neither the expected maximum
nor the expected log partition function. Thus the normalizations agree.
The almost-sure limit in the source also gives the expected limit:
Gaussian union-bound tails for `max_z |H_N(z)|/N` give uniform
integrability.

For every nonnegative, nondecreasing, right-continuous function
`a:[0,1)->[0,infinity)` with `integral a<infinity`, let `Phi_a`
be the zero-temperature Parisi PDE solution. The imported identities are

\[
 P_{\rm SK}=\inf_a\left[\Phi_a(0,0)
                  -\tfrac12\int_0^1t a(t)\,dt\right],    \tag{6}
\]

\[
 \Phi_a(0,0)=\sup_{|u_t|\le1}
 \mathbb E\left[\left|W_1+\int_0^1a(t)u_t\,dt\right|
             -\tfrac12\int_0^1a(t)u_t^2\,dt\right],       \tag{7}
\]

where the supremum uses progressively measurable controls for standard
Brownian motion `W`. Equivalently,
`partial_t Phi=−(partial_xx Phi+a(partial_x Phi)^2)/2`, with terminal
condition `Phi(1,x)=|x|`. Equations (6)--(7), rather than an approximate
decimal value of `P_SK`, are the external input.

## 3. An explicit lower bound for the Parisi constant

Use the Brownian motion `W` of (7) and set

\[
 q(t)=2t-t^2,\qquad
 M_t=\int_0^t\sqrt{2-2s}\,dW_s,\qquad
 S=\operatorname{sgn}M_1,\qquad
 u_t=\mathbb E(S\mid\mathcal F_t^W).                     \tag{8}
\]

The martingale `M` is Gaussian with `E M_t^2=q(t)` and
`E M_1^2=1`. For `t<1`, the independent remaining increment has
variance `1-q(t)=(1-t)^2`, so, writing `Phi` for the standard normal
distribution function,

\[
 u_t=2\Phi\!\left(\frac{M_t}{1-t}\right)-1,\qquad
 |u_t|\le1,\qquad \mathbb E(Su_t)=\mathbb E u_t^2.        \tag{9}
\]

This control is adapted and continuous on `[0,1)`. Almost surely
`M_1` is nonzero, so `u_t` tends to `S` as `t` increases to one.
Consequently, on setting `u_1=S`, it has almost surely continuous
bounded paths and is an admissible progressive control. Its
martingale is uniformly integrable.

To compute its second moment, condition on `M_t` and take two
independent standard normal residuals. The two terminal normal
variables each have variance one and correlation `q(t)`. The normal
sign-correlation identity therefore gives

\[
 \mathbb E u_t^2=\frac2\pi\arcsin q(t)\ge t.             \tag{10}
\]

The sign identity follows directly from rotational symmetry of a
standard planar Gaussian: for correlation `r`, the sign disagreement
region has angular probability `arccos(r)/pi`, so the product
expectation is `1-2 arccos(r)/pi=2 arcsin(r)/pi`. To verify the
inequality in (10), set `s=1-t`. Concavity of sine on `[0,pi/4]`
gives `sin(pi s/4)>=s/sqrt(2)`. Thus

\[
 1-\cos(\pi s/2)=2\sin^2(\pi s/4)\ge s^2,
 \qquad q(t)\ge\sin(\pi t/2).                           \tag{11}
\]

Since `|z|>=Sz` pointwise, insertion of (8) into (7) gives

\[
 \begin{aligned}
 \Phi_a(0,0)-\tfrac12\int_0^1t a(t)\,dt
 &\ge\mathbb E(SW_1)
      +\tfrac12\int_0^1a(t)(\mathbb E u_t^2-t)\,dt\\
 &\ge\mathbb E(SW_1).
 \end{aligned}                                          \tag{12}
\]

All exchanges are justified by boundedness of `S,u` and integrability
of `a`; the last inequality uses `a>=0`. The jointly Gaussian pair
`(M_1,W_1)` satisfies

\[
 \operatorname{Cov}(M_1,W_1)
   =\int_0^1\sqrt{2-2s}\,ds=\frac{2\sqrt2}{3}.
\]

Conditioning `W_1` on `M_1` and using `E|M_1|=sqrt(2/pi)` now yields

\[
 \mathbb E(SW_1)=\sqrt{\frac2\pi}\frac{2\sqrt2}{3}
                =\frac4{3\sqrt\pi}=K_0.                \tag{13}
\]

Taking the infimum in (6) proves

\[
 P_{\rm SK}\ge K_0=\frac4{3\sqrt\pi}>\frac1{\sqrt2}.     \tag{14}
\]

The strict comparison is equivalent to `pi<32/9`, which follows
already from `pi<22/7`. No approximate ground-state constant is used.

## 4. Pointwise reduction to a bipartite random pressure

Block reversal preserves the internal energy difference and reverses
the cross energy. Hence, for each fixed host and cross block,

\[
 \begin{aligned}
 e^{F_{A,B}(t)}
 &=\mathbb E_{x,y}\cosh[\eta(Q_A(x)-Q_A(y))]
                       \cosh(\gamma x^TBy)\\
 &\ge\mathbb E_{x,y}\cosh(\gamma x^TBy)
  =\mathbb E_{x,y}e^{\gamma x^TBy}.
 \end{aligned}                                          \tag{15}
\]

Write the logarithm of the last expression as `L_n(B)`. This lower
bound contains no host. It therefore also bounds `min_A F_{A,B}`,
even if the host is chosen after observing `B`.

First replace the entries of `B` by independent standard Gaussians
`G_ij`. Index the Gaussian process `X_z=x^T G y` by
`z=(x,y)` in `{−1,1}^{2n}`. For two indices write
`q_x=n^{-1}x dot x'`, `q_y=n^{-1}y dot y'`; then

\[
 \mathbb E X_zX_{z'}=n^2q_xq_y.
\]

Independently form

\[
 Y_z=\frac1{\sqrt2}\sum_{i<j\le2n}g_{ij}z_iz_j
                               +\sqrt{n/2}\,g_0.
\]

The common Gaussian term gives exactly

\[
 \mathbb E Y_zY_{z'}=n^2\left(\frac{q_x+q_y}{2}\right)^2,
 \quad \mathbb E Y_z^2=\mathbb E X_z^2=n^2,
\]

and the covariance difference is
`E YY'−E XX'=n²(q_x−q_y)²/4>=0`. Gaussian integration by parts
along `sqrt(s)X+sqrt(1−s)Y` therefore gives

\[
 \mathbb E\log\mathbb E_z e^{\gamma X_z}
 \ge\mathbb E\log\mathbb E_z e^{\gamma Y_z}.              \tag{16}
\]

More explicitly, the derivative of the expected log partition function
is `(gamma²/2)` times the Gibbs average of the nonnegative covariance
difference above; the diagonal covariance difference is zero.
The common `g_0` adds a centered scalar to the log partition function
and drops out of its expectation. With `N=2n`, the remaining SK
inverse temperature is `gamma sqrt(N/2)=c sqrt(t/2)`. The elementary
maximum-term lower bound and (14) now imply

\[
 \begin{aligned}
 \mathbb E L_n(G)
 &\ge c\sqrt{t/2}\,\mathbb E\max_z H_{2n}(z)-2n\log2\\
 &\ge[c\sqrt{2t}\,K_0-2\log2]n-o_{c,t}(n).
 \end{aligned}                                          \tag{17}
\]

This uses neither a claimed formula for the bipartite ground state
nor a reversal through an outer minimum.

## 5. Direct Bernoulli-to-Gaussian comparison at fixed temperature

Consider `L_n(B)` as a smooth function of one real entry `b=B_ij`,
with all other entries held fixed. The associated observable
`chi=x_i y_j` is always `+-1`. Its first three derivatives give

\[
 \partial_b^3L_n
 =\gamma^3\mathbb E_{\rm Gibbs}(\chi-\langle\chi\rangle)^3
 =-2\gamma^3\langle\chi\rangle
                         (1-\langle\chi\rangle^2),
\]

so `|partial_b^3 L_n|<=2 gamma^3` uniformly. Taylor expansion through
degree two about zero has matching expectations for a fair sign and a
standard Gaussian, because both have mean zero and variance one.
The sum of the remainder bounds in one replacement is at most

\[
 \frac{\gamma^3}{3}(1+\mathbb E|g|^3)
 =\frac{\gamma^3}{3}(1+2\sqrt{2/\pi}).
\]

Replace the `n^2` entries successively, keeping all un-replaced
variables independent. This proves

\[
 |\mathbb E_{P_0}L_n(B)-\mathbb E L_n(G)|
 \le\frac{1+2\sqrt{2/\pi}}3 n^2\gamma^3
 =O_{c,t}(\sqrt n).                                     \tag{18}
\]

The comparison is at fixed temperature; no ground-state universality
theorem or exchange of limits is required. Combining (15), (17),
and (18) proves (3), including its uniform host minimum version.

Finally `a_A(beta)<=beta Phi(A)` and the existing complete-signing
upper bound `m_n<=(1/2+o(1))n^(3/2)` give
`2R_n(beta)<=cn+o_c(n)`, proving the excess in (4).

## 6. Relative entropy and selected outcomes

Changing one sign entry of `B` changes `L_n(B)` by at most `2 gamma`.
The independent-coordinate martingale proof of the bounded-differences
inequality, using Hoeffding's lemma at each increment, yields

\[
 \log\mathbb E_{P_0}
   e^{\lambda(L_n-\mathbb E_{P_0}L_n)}
 \le\frac{m\gamma^2\lambda^2}{2}
 \quad(\lambda\in\mathbb R).                             \tag{19}
\]

For any probability law `D` on cross signs, the entropy variational
inequality applied to `−lambda(L_n−E_P0 L_n)` gives

\[
 \mathbb E_D L_n\ge\mathbb E_{P_0}L_n
      -\frac{D(D\Vert P_0)}\lambda
      -\frac{m\gamma^2\lambda}{2},\qquad\lambda>0.
\]

Optimizing in `lambda`, including the zero-entropy limiting case,
and noting `2m gamma²=c²tn`, proves

\[
 \mathbb E_D L_n\ge\mathbb E_{P_0}L_n
                       -c\sqrt{tn D(D\Vert P_0)}.        \tag{20}
\]

Equation (5) follows from the host-free pointwise bound (15).

In particular, if `Delta(c,t)>0` and a proposal law satisfies
`E_D min_A F_{A,B}(t)<=2R_n(beta)+o(n)`, equations (4), (20) force

\[
 D(D\Vert P_0)\ge
 \left[\frac{\Delta(c,t)^2}{c^2t}+o(1)\right]n.
\]

Thus success in mean requires an extensive amount of information
relative to the iid signs. The same necessary condition holds with
any fixed or cross-block-dependent choice of internal host in place
of the minimum.

There is also a direct rare-outcome consequence. Suppose
`Delta(c,t)>0`, fix any `0<kappa<Delta(c,t)`, and let `epsilon_n->0`.
For all sufficiently large `n`, the event

\[
 E_n=\{B:\min_A F_{A,B}(t)
                      \le2R_n(\beta)+\epsilon_n n\}
\]

is contained in `{L_n<=E_P0 L_n−kappa n}`. The lower tail from (19)
therefore gives

\[
 P_0(E_n)\le\exp\left[-\frac{\kappa^2}{c^2t}n\right].    \tag{21}
\]

Consequently `exp(o(n))` samples cannot encounter `E_n` with probability
bounded away from zero if each matrix sample has marginal law `P_0`.
No independence between different samples is needed for this union
bound. In particular, choosing the best among that many iid proposals
does not supply the desired pressure comparison.

Likewise, binary relative entropy gives

\[
 D(E_n)\log\frac1{P_0(E_n)}
 \le D(D\Vert P_0)+\log2.
\]

Thus `D(D||P_0)=o(n)` implies `D(E_n)->0`. An unrestricted selected
cross block is not covered: a point mass has relative entropy
`n² log2`, far outside this hypothesis.

## 7. A conditional consequence for Gaussian-sign rounding

Let `Sigma` be a positive-definite `m by m` covariance matrix with
unit diagonal, and let `D_Sigma` be the law of the componentwise signs
of `G~N(0,Sigma)`, reshaped into a cross block. Componentwise signs of
`N(0,I_m)` have law `P_0`. Data processing for relative entropy and
the exact Gaussian density calculation give

\[
 D(D_\Sigma\Vert P_0)
 \le\tfrac12(\operatorname{tr}\Sigma-m-\log\det\Sigma)
 =-\tfrac12\log\det\Sigma.
\]

Consequently any such family with `-log det Sigma=o(n)` is excluded
both in mean and in probability by Section 6. This is a conditional
test, not a claim that the earlier Gaussian-rounding covariance
always satisfies that hypothesis.

One sufficient condition is explicit. Let
`Sigma=I_m-rho T`, where `T` is symmetric, has zero diagonal and
`||T||_op<=1`, and `0<=rho<1`. For its eigenvalues `lambda_i`, the
trace vanishes. The convergent power series for the logarithm gives

\[
 \begin{aligned}
 -\log\det(I_m-\rho T)
 &=\sum_i\sum_{k\ge2}\frac{\rho^k\lambda_i^k}{k}\\
 &\le\frac{\rho^2}{2(1-\rho)}\sum_i\lambda_i^2
  =\frac{\rho^2}{2(1-\rho)}\|T\|_F^2.
 \end{aligned}
\]

It follows that

\[
 D(D_\Sigma\Vert P_0)
 \le\frac{\rho^2}{4(1-\rho)}\|T\|_F^2.
\]

For fixed `rho<1`, a squared Frobenius norm `o(n)` therefore makes
the rounding family insufficient at the parameters of (4). More
generally the same conclusion follows whenever the displayed upper
bound is `o(n)`. No claim is made here for singular covariances,
uncontrolled determinant loss, or arbitrary dependent sign laws.

## 8. Exact planted-channel formulation and what remains open

The preceding obstruction tests a stronger iid-quenched candidate,
not the general selected-outcome problem. There is an exact identity
that retains the actual paired Gibbs measure and distinguishes them.

Let `nu_eta` be the symmetric mixture of the two opposite-temperature
product Gibbs laws from the finite-step note. Define a probability law
on cross sign matrices by

\[
 Q_{\nu_\eta}(B)=P_0(B)
     \frac{\mathbb E_{\nu_\eta}e^{\gamma x^TBy}}
          {(\cosh\gamma)^m}.                             \tag{22}
\]

Normalization follows because for every fixed spin pair, averaging
the independent signs gives `(cosh gamma)^m`. This is the output law
of the binary rank-one planted channel mixed over the full prior
`nu_eta`. All its probabilities are positive.

The exact finite-step identity
`F=2a_A(eta)+log E_nu exp(gamma x^TBy)` now gives, for every law `D`,

\[
 \boxed{\quad
 \mathbb E_D F_{A,B}(t)
 =2a_A(\eta)+m\log\cosh\gamma
   +D(D\Vert P_0)-D(D\Vert Q_{\nu_\eta}).
 \quad}                                                  \tag{23}
\]

For `D=P_0`, the correction is the reverse relative entropy
`−D(P_0||Q_nu)`, not the forward entropy of the planted law.
The candidate reverse-entropy lower bound needed to pay the entire
radial deficit at actual `a_A(beta)` minimizers is therefore false
at the parameters of (4), by the theorem above.

For genuinely dependent selection, (23) still applies, but the needed
difference of relative entropies has not been bounded in the required
direction. The theorem shows that `o(n)` information relative to iid
signs is insufficient at the stated fixed parameters; it does not
show that `Theta(n)` or larger information cannot succeed. Neither
the correlated Gaussian-sign law nor an arbitrary Gram--Schmidt law
has been shown to satisfy the excluded entropy hypothesis.

The general integral finite-step comparison, and hence the original
MO convergence problem, remain open.

## Provenance

The complete arXiv v2 primary paper, including both the Parisi theorem
and its bounded-control representation, was read before this argument
was recorded. Its downloaded PDF SHA-256 is
`19abfa99c606191e3c33d0c90492b79f5aa1f3a3a7ac32d44ff43aec15bf6978`.
No new mathematical job or finite calculation was run. The explicit
constant and its strict inequality are derived analytically above.
