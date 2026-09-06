# Independent full-source review: all-law second-moment-tail gain

2026-09-06. Role: INDEPENDENT reviewer of the NEW 440-line extension.
I authored its older 553-line all-law prerequisite and other older
imports, and independently reviewed the 325-line cap-free prerequisite.
I did not derive or contribute the NEW extension before this review.
That prerequisite overlap is disclosed; it is not new-source authorship.

## Complete frozen-source read

    /tmp/original_mo_all_law_second_moment_tail_gain.md
    440 lines
    SHA256 c1288ded09bcfef09b6ae9e123b0816e159c10f1fa3745991177814b9bef77e2

I directly read the ENTIRE frozen source in the untruncated ranges
1--220 and 221--440. Verdict: PASS, with no correction required.
Every consequential new link and all stated corollaries were checked.
The source was not edited; this is the only new review artifact.

## New construction and gain checks

1. Exact complete-row second moments give the full nuclear tail bound,
   ell_b=ell-o(1), and an eventual positive lower on ell. The common
   good-coordinate set controls BOTH the missing diagonal energy and
   the signed truncated diagonal. Its complement has o(n) coordinates.

2. Each phase is normalized by its OWN positive diagonal. The printed
   denominator in B_n is positive eventually. Comparison with the common
   h_i normalization is uniform in operator norm; the independent bad
   identity blocks cancel. No source coordinates or variables are removed.

3. The Hermite series and its odd-power trace bound are uniform. Averaging
   the SAME two phase baselines gives tr(M J M_b J). Frobenius pairing
   permits M_b->M at o(n) cost. Literal complete off-diagonal squares and
   Jensen over the common good-coordinate pairs give kappa/(2ell)-o(1)
   with the correct factors. Unequal actual phase diagonals are retained.

4. K_s R_s K_s-P_s has rank at most twice the bad-coordinate count and
   operator norm at most 4C. Its nuclear error is therefore o(n).
   Even Schur powers retain the whole diagonal sum tr P_s. All error
   pairings are uniform in the Hermite index, including the full M_o
   nuclear pairing. Positivity allows L_b->|M|, and the optimized scalar
   quadratic yields the full ell^2-o(1) higher-chaos bound before summing.

5. First-chaos Cauchy--Schwarz gives w_i>=v_i. Actual row lengths bound
   every local variance, while M^2=M_b^2+M_o^2 bounds the averaged c_i^2.
   Neither first-chaos alignment nor a source operator cap is inserted.

6. The imported joint marginal CLT uses ACTUAL complete rows, not rows of
   M_b. Fixed U>C makes both printed finite variance bounds literal.
   Radial clipping, the jointly convex Psi perspective, and its derivative
   bound justify the common gain despite nonconverging local parameters.

7. The bulk clipping identity tolerates its nonzero diagonal. Nuclear
   and second-moment errors control the actual input phase and fields.
   For fixed epsilon, odd smoothed means have bounded covariance by
   Gaussian Poincare, paying the FULL original outlier output energy.
   Independent Boolean rounding is used only for full zero-diagonal M.
   The successive n, epsilon, and U-down-to-C limits are valid; compact
   ell ranges give the last uniform-continuity step. No uniform error
   depending only on C,n is claimed without a spectral-tail rate.

## Corollaries and global quantifiers

8. Cutoffs strictly above 5/3, followed by endpoint continuity, give the
   stated all-law region and 2/5+7/55000 lower. The limiting ell range is
   [3/5,4/5]. The distinction between endpoint atoms and strict cutoffs
   is correct. The all-positive example correctly disproves rank-alone
   second-moment control without additional source-norm information.

9. The stronger rank corollary uses an alpha cap ONLY on a hypothetical
   violating subsequence. The existing norm-only cubic bootstrap has the
   stated normalization; spectral Holder then turns every o(n)-rank tail
   into a vanishing second-moment tail on that subsequence. This proves
   the asserted corollary without an unannounced global alpha cap.

10. The exact nuclear baseline forces ell->1 on any hypothetical
    alpha->1/pi sequence. Completeness then gives mean (|M|-I)^2->0.
    The displayed scalar tail estimate works at EVERY fixed C'>1;
    applying the new theorem and then C' down to 1 produces a strictly
    positive contradiction. This proves an unspecified eventual uniform
    epsilon0>0 above 1/pi. It does NOT identify epsilon0 with
    F_1(1)-1/pi, give F_1(1) as an unconditional floor, or prove convergence.

## Import verification and boundaries

The previously fully read 553- and 325-line prerequisites were rehashed;
the relevant complete CLT, radial, Psi, scalar-region, and smoothing
sections were refreshed in this review sequence. Their hashes remain
0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51 and
0dfa5f62baaa57850a661bbc98d33d32440c783cccb11eaf5446feffbd81f7d4.
I directly refreshed the ENTIRE 176-line phase-moment bootstrap, hash
3736db69d904b5a63ade46b32f6fddcc0505019f45ef483110c3ee67b24c8915,
and the complete nuclear spectral-budget proof underlying saturation.

No mathematical job, canonical edit, frozen-source edit, or backup ran.
The new strict universal lower is proved; convergence and the remaining
original all-order transport question remain OPEN.
