# Independent complete review: universal midpoint and Gaussian shell upper

2026-09-05. Reviewer: exact-proof agent.

The final sources were read completely:

- `/tmp/original_mo_universal_spectral_midpoint_covariance.md`, 424 lines,
  SHA256 `1fc6f5bbb69038b6ac4ed845d26e0724a0ceb0b5a9d96d01b4554a8e37e6f968`.
- `/tmp/original_mo_gaussian_energy_shell_upper.md`, 275 lines,
  SHA256 `8bd3507b722d13077cdb47e8eaa47024b8e95144900226ae4e38272795c5c728`.

**Both PASS. No mathematical corrections required.**

## Universal spectral midpoint

1. The trace of `(aI-A)(A+bI)` gives `ab>=n-1`. The interval
   `[-b,a]` is exactly the feasible interval for the stated covariance
   formula: its same-sign corner numerators impose the two bounds when
   mu is positive; when mu is negative their reversals are incompatible.
2. The positive tensor factorization has the stated coefficients. The
   actual extreme eigenvectors attain both same-sign corner values, and
   bilinearity excludes a larger interior value. Therefore the covariance
   operator norm is exactly `max(u_alpha,v_alpha)`.
3. Differentiation gives the strict opposite monotonicities in (9), so
   their intersection `alpha=(a-b)/2` is the unique operator-norm
   minimizer. Its covariance norm is `(a+b)^2/(a^2+b^2)<=2`, with no
   bound assumed on the source operator norm.
4. The exact disjoint tensor supports of the sign covariance are
   preserved in the arcsine remainder. The estimates
   `L^2<=2mu`, `alpha^2<=mu/2`, and `mu>=n-1` imply the displayed
   absolute `O(1/n)` operator error, uniformly in rho.
5. The Gaussian covariance interpolation uses the current posterior
   Hessian and its trace bound. The free-alpha derivative likewise
   retains the actual posterior and has the correct factors a^2,b^2;
   operator-optimal midpoint and pressure-optimal alpha are not equated.
6. The full epsilon-dependent sign-to-Gaussian theorem applies with
   the absolute covariance bound 2, uniformly in the arbitrary fixed
   internal prior. The auxiliary choices `c=n^(1/22)` and
   `epsilon=n^(-1/11)` give the raw expected-maximum error
   `O(n^(16/11))`, with an absolute constant for every source A.
7. The paired block identity is the ORIGINAL maximum-absolute Boolean
   norm. Every rounded block is admissible, giving the direction
   `m_(2n)<=E Phi(Gaussian pair)+error` even for exact original
   minimizers. The final source correctly leaves the Gaussian dyadic
   target open and explicitly states that unspecified little-o dyadic
   error alone would not establish all-orders convergence.

## Actual Gaussian energy-shell upper

1. Contracting the exact covariance with a fixed y in shell k gives
   `T_k=v_k I+t_k A` with the two displayed coefficients. The positive
   tensor decomposition proves `T_k>=b_0 nI`, including all endpoint
   choices of alpha.
2. On fixed source-energy shells all three positive-sector squared
   radii are constant. The exact increment difference in (9) factors
   into nonnegative products by PSD Cauchy--Schwarz. Adding a common
   Gaussian equalizes point variances without changing expected maxima;
   the soft-maximum interpolation has the asserted upper direction.
3. The two resulting Gaussian field maxima separate exactly. Combining
   both phase signs and every pair of shells uses Gaussian Lipschitz
   concentration, not independence between shell maxima. The remainder
   is exactly `n sqrt(2 Lambda log(2 d_A^2))`; the complete-sign energy
   parity count gives `d_A<=binom(n,2)+1`.
4. The universal midpoint has `Lambda<=1+2/pi`, so this remainder is
   uniformly `O(n sqrt(log n))` without a source operator cap.
5. The constrained-width calculation has the correct mean energy
   `n(n-1)arcsin(t/v)/pi`. Centering A spectrally gives the Hamming
   inequality `|Q_A(x)-Q_A(s)|<=2R sqrt(D(n-D))`; its convex inverse
   and Jensen imply `E D>=n delta` without a sign-independence premise.
6. The half-normal truncation loss uses only the univariate Gaussian
   marginals. The optimizing quantile gives the exact factor
   `2 exp(-q_delta^2/2)-1`. The elementary density bound and
   `delta>=r^2/4` then yield both explicit polynomial deficits in (15).

The shell inequality is a genuine expected-max upper bound. Its current
right side has not been bounded by the source-optimal constant. These
reviews do not assert convergence, a selected signing estimate, or a
replacement of the coupled posterior by the source Gibbs law.
