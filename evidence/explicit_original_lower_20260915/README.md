# Explicit original lower bound: independent numerical corroboration

Evidence for `../NOTE_2026-09-15_EXPLICIT_ORIGINAL_LOWER_BOUND.md`
(SHA-256 `b62d94c49487652f83c1b0b304acda4b7ac1a652db441f605ef101dca9a7f27d`,
including the note's 2026-09-17 addendum).

**Scope.** The all-orders argument is analytic. Nothing here proves any
all-orders statement; these are finite corroborations of the note's
estimates on exhaustive small orders and on adversarial families, plus an
exact-rational certificate for the scalar contradiction. No independent
human review and no formalization is claimed. The review record is
`SELF_REVIEW.md`.

## What is checked

For each signing `A` the scripts reconstruct every object of the note:
`M=A/sqrt(n)`, `S=sign(M)` (kernel `+1`), `|M|`, `ell=tr|M|/n`,
`d=q+1-2*ell`, the rounding frame `R` of Section 2, the exact
sign-Gaussian covariance `C_X=(2/pi)*arcsin(R)`, `v_i=(SC_XS)_ii` and
`E F_i^2=(MC_XM)_ii`, and then test

  * (4) `R` is a correlation matrix with `0 <= R <= 4I`;
  * (5) `||R-(I+S)||_F <= 4 sqrt(n d)`;
  * (6) `e >= kappa ell/2 - 2 kappa sqrt(d) - 2(1-kappa)/sqrt(n)`;
  * (7) `||C_X-(I+kappa S)||_F <= (4 kappa + 5(1-kappa) a_n) sqrt(n d) + (1-kappa)/sqrt(n)`;
  * the diagonal-variance estimate `(1/n) sum |v_i-1| <= (5 + 5(1-kappa)(n^{-1/2}+n^{-1})) sqrt(d) + (1-kappa)/n`;
  * the field-variance estimate `(1/n) sum |sqrt(E F_i^2) - sqrt(v_i)| <= 2 sqrt(d)`;
  * (8) `f >= sqrt(kappa)(1-b_n) - e_4(n)` (Monte Carlo for `f`);
  * (9) `|Q_M(z)| <= n alpha` for the blend `z=(1-p)X+pY`, `p=1/10`
    (exhaustive orders 2..5, 256 sampled states each; `Phi(A)` exact);
  * identities `||M-S||_F^2 = n d`, `sum_i S_ii^2 <= n d`.

Float64 note: entries of `R` can equal `+/-1` exactly, so the entrywise
`arcsin` evaluation loses about `1e-8` (derivative blow-up at the endpoint);
the PSD/diagonal checks use a `1e-6` tolerance, far below any genuine
deviation, and the observed maximal slack is `1.4e-8`.

## Results (soulkiller, 2026-09-17)

Exhaustive, deterministic estimates, zero violations at every order:

| order n | signings | c5 max ratio | c7 max ratio | cv max ratio | cfv max ratio | update9 |
| --- | --- | --- | --- | --- | --- | --- |
| 5 | 1,024 | 0.4116 | 0.2876 | 0.1471 | 0.2873 | max 1.0, 0 viol. |
| 6 | 32,768 | 0.4105 | 0.2658 | 0.1254 | 0.5000 | - |
| 7 | 2,097,152 | 0.4255 | 0.2471 | 0.1231 | 0.2743 | - |

(`ratio` = lhs/rhs for the upper-bound estimates; `c6` is inactive at these
orders because its right side is negative there. `update9` ratio 1.0 means a
corner state attains the Boolean maximum exactly, as the bound requires.)

Families (65 signings; Paley conference orders `q+1` for prime
`q = 1 mod 4`, `q <= 157`; conference with 1/3/10 forced flips; uniform
random orders 8..256; all-ones off-diagonal): zero violations. Max ratios:
`c5 0.3384`, `c6 0.8042` (active on conferences), `c7 0.1940`,
`cv 0.0693`, `cfv 0.5000`, `update9 1.0`.

Conference regime (near saturation): `R = I + S` identically, so `c5`
vanishes (float64: `<= 2.2e-15`); `e` approaches `kappa*ell/2` from above (0.3176 at `n=158`); the
Monte Carlo `f` approaches `sqrt(kappa)` (0.7934 +/- 0.0001 at `n=62`);
`b_n` decays like `7 sqrt(d)`, and the empirical slack `e_4` is zero at
every tested order. `d` decays like `1/(4 n^2)` as expected.

Exact rational certificates (`sharpen_constant.py`; same proof as the note,
tighter enclosures, sharper `d <= 6.3 eps` bookkeeping):

  * claimed: `liminf alpha_n > kappa/2 + eps` certified for
    `eps = 1e-6, 2e-6, 3e-6, 3.5e-6, 4e-6` (margins `+0.4295`, `+0.2699`,
    `+0.1473`, `+0.0938`, `+0.0440` in the normalized `101`-unit chain);
  * controls `eps = 4.5e-6, 5e-6` fail as expected (beyond the ceiling);
  * ceiling of this method: `eps = 305/68378763 = 4.460449...e-6`, certified;
  * the note's own published margin `88577/250000` is reproduced exactly
    from its stated constants (`note_exact_reproduction`).

Exploratory (descriptive, not used in any claim): hill-climb search for
signings maximizing the `(5)` ratio finds `0.4171` at `n=10` and `0.3944`
at `n=14`, so estimate (5) carries real slack.

## Cross-machine replay

`verify_estimates.py` and `sharpen_constant.py` were replayed in fresh
staging directories on two other mesh nodes with the staged hashes verified
before execution:

  * orin (aarch64, 6 cores), `/tmp/qml-explicit-lower-3ZC8qz`,
    exit status 0 (3.5s exhaustive order 6 / 0.5s order 5 / 1.3s families);
  * jellyfin (x86_64, 16 cores), `/tmp/qml-explicit-lower-rcFX6v`,
    exit status 0 (1.9s / 0.3s / 0.9s).

Replay receipts are `replay_orin_*.json`, `replay_jellyfin_*.json`.
Agreement: the exact-rational certificates are bit-identical on all three
hosts; the floating-point aggregates agree to `1e-9` (last-bit
eigendecomposition differences between architectures).

## Artifacts and SHA-256

    verify_estimates.py      e6db5e83c4feeffdd7f9f787217ade6e53a9f4325d617786f9c13351ac6cb635
    sharpen_constant.py      6ac88f483bee1e1ce10f4fccbfec77c788f6189859a9881d5bc564b27df6955b
    check_constants.py       8a7d23b8620762f12a6da2bc80755269f541e9bda7dd1f06a99dca18f9fe4c29  (author's)
    result_estimates_ex5.json   d470abe5c3bbbbb53de5a2e2423d28f0832169e3be822adedc343dc56f1a4ff5
    result_estimates_ex6.json   c01337dc57b0a6501f80919c3dbc887a1e5899ba4320a457f093cc11fb7b858f
    result_estimates_ex7.json   8f67492a9d16df74e6f5be67101141d18c598c67ba9756ea13c2600fd180d7f4
    result_families.json        0e07fc68d1c26797c0e6b62f36d78fc79d818214486a39323a9d57101fa4f56b
    result_hillclimb_n10.json   54e3138d7822a05947f949cfa7aca5b175cdda14098619773707a73a35df04f6
    result_hillclimb_n14.json   17a661449d6cfbfbdb3c0149dc972aff8d98cb0f3bee9598223c75f2a5bca8b3
    result_sharpen.json         e1c1c9c5cd5e60a478a862a5768b433f27e9c6bf5ce2afcd9b057b1fd20f0fc6

Toolchain: Python 3.14.4, numpy 2.4.4 (soulkiller); numpy 2.5.2 (orin),
2.3.5 (jellyfin).

## Reproduce

    python3 verify_estimates.py --exhaustive 7 --workers 84 --out result_estimates_ex7.json
    python3 verify_estimates.py --exhaustive 5 --update-check --workers 16 --out result_estimates_ex5.json
    python3 verify_estimates.py --families --mc-samples 200000 --update-check --out result_families.json
    python3 verify_estimates.py --hillclimb 14 --restarts 16 --steps 30 --out result_hillclimb_n14.json
    python3 sharpen_constant.py --out result_sharpen.json
    python3 check_constants.py            # the author's scalar checker

`verify_estimates.py` is self-contained (numpy + stdlib) so it can be
replayed on any node; the Paley builder inside it is validated against
Gram checks (`C^2 = q I`, error 0) and against the repo
builder `scripts/paley_structural_lift_family.py` at `q=13` (exact match).

## Limits

  * finite checks cannot prove the all-orders estimates; they corroborate
    them and pin their finite-order slack;
  * the Monte Carlo layer estimates `E|F_i|` with standard errors reported
    in `result_families.json`; no MC-based claim is used in the certificate;
  * the sharpening is arithmetic only (same proof, tighter enclosures); it
    changes no analytic estimate and does not raise the method's ceiling
    beyond `4.46e-6`.
