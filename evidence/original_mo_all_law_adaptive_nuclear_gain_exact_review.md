# Contributing full-source review: all-law adaptive nuclear gain

2026-09-06. Reviewer/contributor: optimized_profile_exact.

## Frozen source and reads

I directly read the ENTIRE final 553-line source
`/tmp/original_mo_all_law_adaptive_nuclear_gain.md`, SHA256
`0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51`.
No source correction is required.

All named prerequisites were fully read in this analytic work sequence,
and all hashes were directly rechecked for this review:

- Near-flat source612: `7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298`.
- Nuclear budget147: `ee8ad5ff3dbf9aa9e251c4190e98ee1671c9a2140c759ba6f768f8c9c03ef13d`.
- Complete-cross Gaussian proof411: `b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
- Original phase262: `7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
- Pi-enclosure source444: `106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

The 147-line nuclear budget was newly refreshed and fully read during
this review; the other complete reads also underlie the preceding full
near-flat review. The scalar209 fixed-update lemma is not a logical
premise of this new adaptive argument. No pi certificate was rerun.

## Contribution disclosure

This is NOT an independent whole-new-source review. I derived the
actual diagonal-normalized frame control and all-law higher-chaos mean,
supplied the finite positive-part bound now in (5.4), and checked the
adaptive identity and scalar proposals before the source was written.
I also contributed to the older Gaussianization/phase/trace chain.
The adaptive clipping, radial bound, convex averaging, and full source
were authored by the proof worker. Their written implementation is
checked here, with my contributing role retained. The separate new
cross-to-nuclear transfer is my work but is not a premise of source553.

## Substantive full-source checks

1. Actual h_i=|M|_ii lie in [q/C,sqrt(q)], so both unequal-diagonal
   spectral phases are genuine correlations with norm at most 2C^2/q.
   The exact oriented phase average in (4.1) has the correct nuclear
   baseline kappa q/(2ell), enabling an added gain on the SAME phases.

2. The exact diagonal-congruence identity in (5.3), the complete-entry
   O(sqrt(n)) trace error, and M^2>=t|M|-t^2I/4 prove (5.4) uniformly
   over all odd chaos orders. Summing positive coefficients is valid.
   The conditional variances obey w_i>=v_i, mean w_i at least
   (1-kappa)ell^2-o(1), row sigma_i^2<=2C^2, and mean c_i^2<=kappa C^2.

3. The reused contractions establish the needed joint two-variable
   Gaussian limit for arbitrary local covariance, including singular
   limits. Only the LINEARLY growing clipped gain is passed through
   that limit. Bounded second moments supply its uniform integrability;
   no unsupported uniform integrability of F_i^2 is asserted.

4. Conditional independent masks with epsilon_i=min(r_i/(2C),1)
   give the exact adaptive quadratic expansion. Zero diagonal removes
   all mask-variance terms. The spectral penalty is the actual C,
   and the resulting gain is precisely g_C(r), with admissible masks.

5. The Rayleigh radial identity (8.1) is correct. Since the per-angle
   coefficient squared is at most sigma_i^2<=2C^2, the clipping loss
   is at most exp(-1) times the Gaussian negative-field second moment,
   for either covariance sign and degenerate pairs alike.

6. The angular integral, negative-covariance comparison, and continuous
   perspective Psi(a,w) are correct. Its stated Hessian is PSD, it
   decreases in a and increases in w, so Jensen yields (9.5) without
   local homogeneity or first-chaos alignment. Averaging the two phase
   lower bounds adds this common gain to the exact nuclear baseline.

7. Integrating the w derivative before applying the arctangent lower
   proves the positive envelope (10.1). The exponential-series bound,
   product bound on kappa(1-kappa), denominator comparison
   11296/5625<81/40, and exact 3/1250 gain are all correct.
   The ell<=3/4 branch is covered separately by the nuclear bound.
   Compact subsequences and C=5/3+epsilon followed by epsilon->0
   correctly handle both limsup hypotheses, without moment convergence.

Verdict: full-source contributing review PASS, no correction. The
uniform all-law gain and the entire region (1.4) are proved as stated;
the latter has lower 2/5+7/55000. This removes spectral-shape and
diagonal-homogeneity premises, not the actual operator cap. It does
not settle the complementary nuclear region or global original MO limit.

No mathematical computation, solver, checker, numerical integration,
optimization, or signing search was run on any host. Only complete
reads, hashes, and this /tmp review write were performed. No canonical
repository edit, publication, or backup was performed by this reviewer.
