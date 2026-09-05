# Independent complete review: expected original paired norm

2026-09-05. Reviewer: exact-proof agent.

Reviewed all 265 lines of
`/tmp/original_mo_expected_paired_norm_gaussian_equivalence.md`, SHA256
`bff778718c0f357598c035edba4598f2ed67b1c49359c668958afe1c39207df3`.

**PASS. No corrections required.**

Analytical checks:

1. The canonical covariance has operator bound `1+8K^2` for every
   `n>=2`, uniformly in the actual opposite temperatures used to
   determine alpha. Those temperatures are correctly fixed separately
   from the later auxiliary zero-temperature parameter c.
2. The internal energy is absorbed into a deterministic prior, which
   may vary with c. The reviewed sign-to-matched-Gaussian theorem is
   uniform in that prior, so its full epsilon-dependent estimate may
   be applied with the growing choices in Section 4.
3. The exact tensor covariance linearization gives the displayed
   `64(1-2/pi)K^2/n^2+16(1-2/pi)K^4/n` operator error. Gaussian
   interpolation then contributes at most `D_K c^2` to pressure,
   or `D_K c/n` to the normalized maximum estimate.
4. Both maximum-term defects in (8) lie in the SAME interval of
   length `(2n+1)log 2`. Their difference costs only that length.
   Division by `beta n^(3/2)=cn` yields every term in (9).
5. For `c=n^(1/22)`, `epsilon=n^(-1/11)`, the four main normalized
   terms all equal `n^(-1/22)`. The other two are smaller. Hence the
   raw error is `n^(3/2-1/22)=n^(16/11)`.
6. The constants in the full comparison have growth at most
   `D(1+K^4)`: the Gaussian covariance cap is `1+8K^2`, and only its
   square appears in the smooth OU term. The exact tensor error has
   the same maximum fourth power of K. This justifies using a slowly
   growing cap, not only a fixed-K asymptotic statement.
7. The cited same-order regularization applies to ORIGINAL norm
   minimizers and provides both inequalities in (10). For
   `r_n=n^(1/99)`, its objective loss is `n^(-1/198)` and the
   expected-norm comparison is also of that order, because
   `4/99-1/22=-1/198`.
8. The paired block quadratic form is exactly
   `Q_A(x)-Q_A(y)+x^TDy`; thus its maximum absolute value is precisely
   the original Boolean norm, with no width or half-product
   substitution.

All statements correctly concern expectations of maxima. Neither a
coupled pointwise estimate nor a selected signing guarantee is inferred.
The Gaussian expected norm still requires the separate upper comparison.
