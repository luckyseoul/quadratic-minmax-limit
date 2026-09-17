# Whole-proof self-review: NOTE_2026-09-15_EXPLICIT_ORIGINAL_LOWER_BOUND.md

Reviewer: 2026-09-17 operating session (machine review; this is not
independent human review). Reviewed artifact SHA-256 (as committed,
including the 2026-09-17 addendum):

    NOTE_2026-09-15_EXPLICIT_ORIGINAL_LOWER_BOUND.md
    b62d94c49487652f83c1b0b304acda4b7ac1a652db441f605ef101dca9a7f27d
    (recomputed and pinned in result.json)

(At first review, before this session appended its own addendum, the note
hashed `04b1314f30f0de32ef717faf98f6c65b1e7276fe1613971285c869fbe08d7a51`.)

Scope: this review re-derives each display in the note and checks the two
imported inputs against their source notes. It is a reading of the same
argument by the same lineage of machine sessions, not an independent
human referee report. The finite computations cited below cannot prove any
all-orders statement.

## 1. Imports

- (2) `alpha >= kappa q/(2 ell)`: taken from
  `NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md`, display (4.1). Checked:
  the two-phase average, the arcsine monotonicity step, the AM-GM step
  `1/sqrt(h_i h_j) >= 2/(h_i+h_j)`, the Cauchy-Schwarz step
  `sum_{i<j} 1/(h_i+h_j) >= k^2/((n-1) n ell)` with `k=n(n-1)/2`, and the
  trace identity `ell = (1/n) sum_i h_i = tr|M|/n`. All correct as displayed.
- (3) `|E|F| - sqrt(kappa) sqrt(E F^2)| <= e_4(n)`: taken from the Lemma in
  `NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md`, Section 6,
  display (6.1), with the uniformity over `|d_j| <= 1/sqrt(n)`,
  `sum d_j^2 <= 1`, `||R||op <= C`. The note applies it with `C = 4` to the
  rows of `M`, whose coefficients are `M_ij in {+-1/sqrt n, 0}`: the two
  coefficient conditions hold exactly (`max_j |d_j| = n^{-1/2}`,
  `sum_j d_j^2 = (M^2)_ii = q <= 1`). Uniformity in the source (sup over the
  data at fixed C, error -> 0) is what the argument needs; the source states
  it. The full proof of the source lemma was not re-verified here.
- Rational enclosure `7/11 < kappa < 16/25`: rechecked (`7/11 < 2/pi` and
  `2/pi < 16/25`). The sharper enclosures used in the sharpening below come
  from the frozen enclosure `31415926/10^7 < pi < 31415927/10^7` recorded in
  `NOTE_2026-09-05_SOURCE_CROSS_NUCLEAR_TRACE_BOUNDARY.md`.

## 2. Display-by-display verdicts

- Section 2, frame construction and (5): re-derived by hand.
  * `d = q + 1 - 2 ell = ||M-S||_F^2/n` (Frobenius identity via `MS = |M|`,
    `tr S^2 = n`).
  * `sum_i S_ii^2 <= n d`: write `S_ii = sum_k eps_k w_k(i)` with
    `eps_k = sign(lambda_k)`, `w_k(i) = u_k(i)^2`; `sum_k w_k(i) = 1` and
    `(M)_ii = sum_k lambda_k w_k(i) = 0` give
    `S_ii = sum_k (eps_k - |lambda_k|) w_k(i)`; Jensen over `w(i)` then
    `sum_i S_ii^2 <= sum_k (1-|lambda_k|)^2 = n d`. Correct.
  * The normalization inequality `|(1+t)^{-1/2}-1| <= 2(sqrt2 - 1)|t|` on
    `|t| <= 1/2`: the stated rationalization
    `|t|/(sqrt(1+t)(1+sqrt(1+t)))` with the denominator minimized at
    `t = -1/2` gives exactly `2(sqrt2 - 1)`. Correct.
  * Bookkeeping of the three blocks (good `[I,I]` block, removed cross
    blocks, `-S[B,B]`): total squared Frobenius cost `<= 16 sum_i S_ii^2`
    using `2|B| <= 8 sum_{i in B} S_ii^2` (from `|S_ii| > 1/2` on B) and
    row norm one for the cross blocks. Hence (5) with constant 4. Correct.
  * Numerical corroboration: (5) holds with zero or small ratios on all
    signings of order `n <= 7` (2,097,152 signings), all tested Paley
    conference families and random/perturbed families; for exact conference
    signings `R = I + S` identically, so (5) is an identity there (numeric
    residual at most `2.2e-15`).
- Section 3, (6): re-derived. `tr(MR) = tr(M(I+S)) + tr(M(R-(I+S))) = n ell
  + O(||M||_F ||R-(I+S)||_F) = n ell - 4 n sqrt(q d)` worst case; the
  higher-chaos tail obeys `|tr(M R^{circ k})| <= n^{-1/2} tr R^2 <= 4 sqrt n`
  uniformly in `k` (using `|R_ij| <= 1`, `R <= 4I`, `tr R = n`). Dividing by
  `2n` gives (6). Correct. (Note the estimate is inactive for large `d`; it is
  the near-saturation bound.)
- Section 4, (7) and its diagonal/field consequences: re-derived.
  * `C_X - (I + kappa S) = kappa (R - I - S) + H`, `H = C_tail - (1-kappa) I`
    with zero diagonal; `|H_ij| <= (1-kappa)|R_ij|^3`;
    `r^3 <= (1+h+h^2)|r-h| + h^3` for `h = n^{-1/2}` (immediate for `r <= h`,
    else factor `r^3-h^3`); `|||R_off| - |M|_entry||_F <= ||R_off - M||_F <=
    5 sqrt(n d)`. Gives (7). Correct.
  * Diagonal part: conjugation by `S` fixes `I + kappa S` and preserves the
    Frobenius norm; `(1/n) sum |(S C_X S)_ii - 1 - kappa S_ii| <= ||SC_XS -
    (I+kappa S)||_F/sqrt(n)`; `(1/n) sum |S_ii| <= sqrt(d)`. The arithmetic
    `5 kappa + 5(1-kappa) a_n = 5 + 5(1-kappa)(n^{-1/2}+n^{-1})` is exact.
    Correct.
  * Field variance: `(1/n) sum |sqrt((M C_X M)_ii) - sqrt((S C_X S)_ii)| <=
    n^{-1/2} ||C_X^{1/2}(M-S)||_F <= 2 sqrt d` via `C_X <= 4I`. Correct.
- Section 5, (9): re-derived. `z = (1-p)X + pY` with `Y = sign(MX)`;
  rounding identity `E Q_M(z) <= n alpha` for every `z in [-1,1]^n`
  (independent coordinate rounding; multilinearity); `X^T M Y = ||MX||_1`
  by symmetry (`sum_j Y_j (MX)_j`); `Q_M(Y) >= -n alpha`. Gives
  `(1+p^2) alpha >= (1-p)^2 e + p(1-p) f`. Correct. Numerical check of
  `|Q_M(z)| <= n alpha` on exhaustive orders `n <= 5` (1024 signings,
  256 sampled states each): zero violations, extremal ratio exactly 1.
- Section 6, contradiction and (11): re-derived. `d <= 4 eps/(kappa + 2 eps)`
  from `ell >= kappa q/(2 alpha)`, `q <= 1`; `4/kappa < 44/7 < 7`; the final
  rational arithmetic reproduces `88577/250000` exactly (see
  `sharpen_constant.py`, `note_exact_reproduction`).

## 3. Findings

1. No mathematical error was found in the note's chain as displayed, under
   the two imported inputs.
2. The stated constant is conservative: with the tighter rational enclosure
   of `kappa` (from the frozen `pi` enclosure) and the sharper `d`-bound
   `d <= 6.3 eps`, the same argument certifies `liminf alpha > kappa/2 + eps`
   for every `eps <= 4e-6`, with exact rational margins; the same-method
   ceiling is `4.4604e-6` (see `SHARPENING` section of README.md and
   `result_sharpen.json`). No analytic estimate is changed.
3. The note's error/slack budgets are dominated by the `sqrt(d)` terms; the
   method's ceiling scales like `(income/coeff)^2/6.3`, so no improvement of
   this constant by rational bookkeeping alone can exceed `~4.5e-6` without
   sharpening (5), (7), or (8).

## 4. Boundaries

- Not verified here: the full proof of the source Gaussianization lemma; the
  analytic estimates (5)-(8) as all-orders statements (they were re-derived
  by hand and checked numerically on finite families only); any formal
  machine-checked proof.
- The finite checks cover: all signings of orders 2..7 (deterministic
  estimates); 65 family signings (conference, perturbed conference, random,
  all-ones) with deterministic checks plus Monte Carlo for (8); exhaustive
  update-bound checks at orders 4..5.
