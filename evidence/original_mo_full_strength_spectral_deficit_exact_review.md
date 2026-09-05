# Independent review of the full-strength spectral-deficit criterion

2026-09-05. Reviewer: optimized_profile_exact.

Reviewed the complete 190-line proof
`/tmp/original_mo_full_strength_spectral_deficit.md`, SHA-256
`c37123564ec9bba8c8f16048a3ce0d1a40348990cd82d87070f223ce7aa51ad6`.

**PASS.** No correction required.

- Adding independent Gaussian covariance
  `D=[(1-r)I-Sigma]_+` gives the claimed coordinate correlation and
  exact angular mismatch `arctan(sqrt(D_ee))/pi`, including singular Sigma.
- The sign-pressure Lipschitz bound and Cauchy--Schwarz give precisely
  `(2 gamma/pi) sqrt(n^2 V)`.
- A centered biased sign of mean m has third absolute moment
  `1-m^4<=1`; the Gaussian replacement therefore has uniform
  `O(gamma^3 n^2)` error. Convexity and global evenness justify removing
  the conditional mean only after Gaussian replacement.
- The variance loss is `gamma^2 sum m_e^2/2`. The displayed derivative
  of `arcsin((r+s)/(1+s))` is correct and at most one, giving the stated
  additional `gamma^2 V/pi` term.
- Every n factor in the final pressure floor is correct. Under a
  successful mean comparison, the subsequence argument proves the exact
  necessary lower limit `V/n >= pi^2 g_r(c,t)^2/(2 c^2 t)`.
- The actual covariance operator bound supports exponential concentration
  for the repaired law. The expected coupling loss transfers this only
  to an o(1) success probability for the original law when V=o(n).
  It does not exclude exponentially many original-law proposals.
- The sufficient spectral hypothesis is not proved for exact global
  minimizers. No convergence or unrestricted-selection claim follows.

No signing census, simulation, or mathematical computation was run.
