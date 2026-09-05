# Full-strength canonical failure on leading half-product near-minimizers

2026-09-05. **An actual near-minimizer consequence of the new exact-minimum
operator theorem.** This upgrades the construction-cap application, not
the exact-minimizer conclusion. No new numerical computation is used.

Fix c>0 and 0<r<1 such that

\[
 g_r=c(\sqrt2 K_0-1)-2\log2-\frac{c^2}{2\pi}\arcsin r>0,
 \qquad K_0=\frac4{3\sqrt\pi}.
\]

Such fixed parameters exist as established in the spectral-deficit note.
Along an unbounded sequence of orders N there exist complete signings A_N
with, at beta_N=c/sqrt(N),

\[
 0\le a_{A_N}(\beta_N)-R_N(\beta_N)=o(N),                     \tag{1}
\]

whose actual canonical full-strength Gaussian-sign cross law obeys

\[
 \mathbb E\min_{A'}F_{A',B}(1)
       -2a_{A_N}(\beta_N)\ge[g_r-o(1)]N.                     \tag{2}
\]

The probability of reaching that host's own paired endpoint plus o(N)
tends to zero. The covariance-generating host is fixed before B is drawn.
There is no assertion about exact minima or about norm near-minimality.

## 1. Reuse of one twin module above an exact minimizing base

For each n choose an exact half-product minimizer B_n at c/sqrt(n).
The new exact-minimum operator theorem gives

\[
 D_n:=\|B_n\|_{\rm op}+\sqrt n=o(n^{3/4}).
\]

Set

\[
 \ell=2\left\lceil\frac{\sqrt{D_n n^{3/4}}}{2}\right\rceil,
 \qquad m=2\ell,\qquad N=n+m.
\]

Then `ell/D_n` tends to infinity, ell=o(n^(3/4)), and N/n tends to one.
Attach exactly ONE twin module from Sections 1--2 of
`evidence/NOTE_2026-09-05_NEAR_MINIMIZER_OPPOSITE_PHASE_COUNTERFAMILY.md`:
two communities of ell vertices, internal pattern `[[1,1],[1,-1]]`
with zero diagonal, and old-to-twin edges `R_(pi)(1,-1)`. There are
no intermodule edges to choose. The reviewed elementary net argument
supplies a filler realization with operator norm at most 16 sqrt(N).

The identical uniform energy estimate from that proof reads

\[
 |Q_{A_N}(z,x)-Q_{B_n}(z)|
 \le E_N:=\ell^2+16N\sqrt{2\ell}+16\ell\sqrt N
                                      =o(N^{3/2}).            \tag{3}
\]

All filler signs are ordinary complete-host signs. No thermal event is
required, and the old signing remains unchanged.

## 2. The pressure is genuinely leading-order minimal at the final scale

For an exact pressure minimizer C_n at any tau with tau sqrt(n) in a
fixed compact subset of (0,infinity), random-edge comparison and Jensen
give `a_(C_n)(tau)<=binom(n,2) log cosh(tau)=O(n)`. The elementary
extreme-configuration bound gives

\[
 \tau\Phi(C_n)\le2a_{C_n}(\tau)+n\log2,
 \qquad \Phi(C_n)=O_c(n^{3/2}).                               \tag{4}
\]

Each a_C is Phi(C)-Lipschitz in temperature. Apply (4) to minimizers
at c/sqrt(n) and at beta_N. Comparison of their two minimizing values
then gives

\[
 0\le a_{B_n}(\beta_N)-R_n(\beta_N)
 \le O_c(n^{3/2})\left|\frac c{\sqrt n}-\frac c{\sqrt N}\right|
 =O_c(m)=o(N).                                               \tag{5}
\]

Uniform spin averaging on deleted vertices gives R_n(beta)<=R_N(beta)
at the same beta. Equation (3) gives
`|a_(A_N)(beta_N)-a_(B_n)(beta_N)|<=beta_N E_N=o(N)`.
Together with (5) these prove (1). They also give
Phi(A_N)=O_c(N^(3/2)) by (3)--(4). The established construction cap
`2R_N(c/sqrt(N))<=cN+o(N)`, also used in the spectral-deficit note,
therefore implies

\[
 2a_{A_N}(\beta_N)\le cN+o(N).                               \tag{6}
\]

## 3. Only two directions have a fixed negative normalized eigenvalue

The invariant community-constant matrix is exactly
`[[ell-1,ell],[ell,-ell+1]]`, with simple eigenvalues +M and -M,
where M=sqrt(2ell^2-2ell+1). All remaining eigenvalues have magnitude
at most `||B_n||_op+1+16sqrt(N)=o(ell)`. Thus the two extremes of the
full host are exactly +M,-M for large n.

Use the ACTUAL opposite-temperature Gibbs covariances to form the
canonical centering alpha and tensor H. Their phase energy means have
opposite signs, so `|alpha|<=Phi(A_N)/N=O_c(sqrt(N))=o(M)`.
The exact normalization is mu=M^2 because the extreme eigenvalues
have equal magnitude. For T=H/mu its eigenvalues are

\[
 \frac{xy-\alpha(x+y)}{M^2},
\]

over ordered pairs of eigenvalues x,y of A_N. Those with at least one
nonextreme eigenvalue tend uniformly to zero. The two same-sign extreme
pairs have eigenvalues `1 +/- 2alpha/M`, tending to one. The two mixed
extreme pairs have eigenvalue exactly -1. Hence, for every sufficiently
large n,

\[
 V_r=\operatorname{tr}[-rI-T]_+=2(1-r).                       \tag{7}
\]

This exact full-strength deficit uses the actual centering, not a
surrogate phase law. In particular it is o(N).

## 4. Full-strength quenched consequence and scope

Apply `evidence/NOTE_2026-09-05_FULL_STRENGTH_SPECTRAL_DEFICIT.md` to (7).
Its host-free floor is `[c+g_r-o(1)]N`, and (6) proves (2), even if
the internal host A' is selected after observing B. The spectral repair
coupling has expected pressure loss O_c(sqrt(N)), so its probability
transfer bounds success by O_c(N^(-1/2)) plus an exponentially small
repaired-law term at any fixed positive linear margin below the floor.

This is a leading half-product near-minimizer obstruction. The sparse
module is not asserted to preserve exact optimality; the exact-minimum
operator theorem is used only for the old base. No passage to norm
near-minimality, no exclusion of all selected outcomes, and no original
MO convergence or cross-order transport theorem is claimed.
