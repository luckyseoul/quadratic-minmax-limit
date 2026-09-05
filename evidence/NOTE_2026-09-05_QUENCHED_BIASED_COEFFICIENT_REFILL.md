# Quenched Gaussian replacement for actual biased coefficient refills

2026-09-05. **Explicit finite-refill comparison with a subextensive
remainder.** This replaces independent interaction COEFFICIENTS, not the
Boolean spin prior. The entire log-partition posterior is retained.
It does not evaluate the resulting Gaussian expectation or establish
an order-transport bound.

## 1. Fixed profile and the two refill distributions

Let N=2n, c,lambda>0, 0<=t<=1, and let

\[
 M=\frac{\sqrt{2-t}A_I+\sqrt t A_C}{\sqrt N},\qquad
 F(M)=\log\mathbb E_{\sigma,x}
                  e^{\sigma c x^TMx/2},\qquad
 V(M)=\operatorname{tr}M^4.
\]

Both sigma and the Boolean spins have uniform reference laws, so this
is exactly the symmetric log-cosh pressure. The row squared norm is
`d=1-(2-t)/N<=1`, and every edge magnitude is at most `sqrt(2/N)`.

Choose any unordered edge set E. For each e in E choose a deterministic
bias rho_e in [-1,1]. All signs and Gaussians below are independent.
The signed refill M^R keeps the other edges fixed and replaces

\[
 m_e\longmapsto X_e=m_eR_e,\qquad
 \Pr(R_e=1)=\frac{1+\rho_e}{2}.
\]

The Gaussian refill M^G instead replaces these coefficients by

\[
 Y_e=\rho_e m_e+|m_e|\sqrt{1-\rho_e^2}\,G_e,
 \qquad G_e\sim N(0,1).
\]

Thus X_e and Y_e have the same mean and second moment. Both matrices
are symmetric and zero diagonal. M^R is an admissible complete signing
at the original profile weights; M^G need not be a signing. Singular
one-coordinate cases rho_e=+/-1 and zero profile weights are included.

## 2. Uniform quenched pressure replacement

For a fixed real choice of every other coefficient, differentiating
the FULL log-partition function with respect to m_e gives

\[
 \partial_e F=c\langle\tau_e\rangle,\qquad
 \partial_e^3 F=c^3\langle(\tau_e-\langle\tau_e\rangle)^3\rangle,
 \quad \tau_e=\sigma x_i x_j\in\{-1,1\}.
\]

If r is the current actual posterior mean of tau_e, its third centered
moment is `-2r(1-r^2)`. Consequently

\[
 |\partial_e^3F|\le2c^3                                      \tag{1}
\]

uniformly over ALL real coefficient values, including intermediate
Gaussian ones. The posterior is differentiated in (1), not frozen.

Taylor expansion at the common mean, with matching first two moments,
therefore gives a one-coordinate replacement error at most

\[
 \frac{c^3}{3}\left(
   \mathbb E|X_e-\rho_em_e|^3+
   \mathbb E|Y_e-\rho_em_e|^3\right).
\]

The absolute third moments are exactly

\[
 \mathbb E|X_e-\rho_em_e|^3=|m_e|^3(1-\rho_e^4),
\]
\[
 \mathbb E|Y_e-\rho_em_e|^3
 =2\sqrt{2/\pi}\,|m_e|^3(1-\rho_e^2)^{3/2}
 \le2\sqrt{2/\pi}\,|m_e|^3(1-\rho_e^4).
\]

Condition on all other independent coefficients at each replacement
and telescope. With `C_3=(1+2sqrt(2/pi))/3`, this proves

\[
 \boxed{\quad
 |\mathbb E F(M^R)-\mathbb E F(M^G)|
 \le\varepsilon_E:=C_3c^3
             \sum_{e\in E}|m_e|^3(1-\rho_e^4).
 \quad}                                                     \tag{2}
\]

This compares expectations of LOG pressures. It does not move the
coefficient expectation outside the logarithm, integrate a conditional
moment-generating upper bound, or replace a posterior covariance by
its value in the original host.

## 3. Exact quartic correction

The quartic expectations have an exact difference:

\[
 \boxed{\quad
 \mathbb E V(M^G)-\mathbb E V(M^R)
 =4\sum_{e\in E}(1-\rho_e^4)m_e^4.
 \quad}                                                     \tag{3}
\]

Indeed, expand the trace over closed walks of length four. A walk using
only one edge contributes its fourth power, twice for each unordered
edge. Every other nonzero walk either uses two distinct edges twice
each or is a simple four-cycle. Their expectations agree under the two
refills, since independence and the matching first and second moments
are enough. There is no nonzero closed four-walk with an edge of
multiplicity three and a distinct edge of multiplicity one in a
zero-diagonal matrix. Finally,

\[
 \mathbb E Y_e^4=(3-2\rho_e^4)m_e^4,
 \qquad \mathbb E X_e^4=m_e^4,
\]

which gives (3), including arbitrary unchanged coefficients elsewhere.

## 4. An actual-optimizer finite-noise inequality

Now assume A is an ACTUAL global minimizer of `F(M)+lambda V(M)`
over complete signings at the fixed c,lambda,t. Every outcome M^R is
admissible. Apply global optimality, then (2)--(3), to obtain

\[
 \boxed{\quad
 F(M)+\lambda V(M)
 \le\mathbb E\bigl[F(M^G)+\lambda V(M^G)\bigr]
       -4\lambda\sum_{e\in E}(1-\rho_e^4)m_e^4
       +\varepsilon_E.
 \quad}                                                     \tag{4}
\]

The expectation on the right is the actual finite Gaussian-refill
log-pressure. It remains an expectation of the full posterior problem;
no Gaussian annealing or prior quadratic proxy has replaced it.
The force-kernel or an upper local-field moment bound is not required
for this particular remainder.

For the edges incident to a single vertex, (2) gives
`epsilon_E<=C_3 c^3 d sqrt(2/N)=O_c(N^(-1/2))`.
More generally, if E is the set of edges incident to S, then

\[
 \sum_{e\in E}m_e^2\le d|S|,\qquad
 \varepsilon_E\le C_3c^3d\sqrt{2/N}\,|S|,
 \qquad
 4\lambda\sum_{e\in E}(1-\rho_e^4)m_e^4
                          \le\frac{8\lambda d}{N}|S|.    \tag{5}
\]

For all edges together the sharper row count gives

\[
 \varepsilon_E\le C_3c^3d\sqrt{N/2}=O_c(\sqrt N),\qquad
 4\lambda\sum_{e\in E}(1-\rho_e^4)m_e^4\le4\lambda d.
                                                               \tag{6}
\]

Thus even a macroscopic independent coefficient refill is compared at
subextensive quenched error. This does NOT mean that the refill changes
the pressure by o(N); only the difference between its two distributions
has that bound.

## Scope of the remaining comparison

Equation (4) is a necessary finite-noise variational inequality for
every actual optimizer and every specified collection of independent
coefficient biases. It retains the exact spectral correction and the
entire quenched Gaussian expectation. The biases may be chosen after
the deterministic optimizing host is fixed, but cannot depend on the
independent refill variables within this argument.

To obtain order transport one still needs a favorable comparison of
the Gaussian expectation in (4) with the appropriate smaller-order
endpoint. Neither (2) nor (4) supplies it. In particular, an independently
refilled Gaussian law is not the correlated, pair-dependent boundary
law of full-strength canonical sign rounding. No such law substitution
or conclusion about its covariance integral is made here.
