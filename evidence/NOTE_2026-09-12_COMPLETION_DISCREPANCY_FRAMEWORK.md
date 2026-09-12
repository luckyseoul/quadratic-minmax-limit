# The completion-discrepancy framework for the doubling estimate

**Prior-art notice (added 2026-09-12 after the user's audit request).** This
note's framework largely duplicates pre-existing material. Verified against
the repository: `NOTE_2026-09-01_RG2_EQUAL_ENDPOINT_PALEY_SHIELD.md` already
contains, in strictly stronger form, (i) the cross-term floor — its eq. (3):
`max|x^T R y| = max_x ||Rx||_1 >= (sqrt(2/pi)+o(1)) n^{3/2}` for skew `R`;
(ii) the independent-budget obstruction (its §2); (iii) the exact residual
reduction to a band-style condition (its eq. (14)); plus the equal-endpoint
doubling frame, the balanced near-conference skew construction with a
**proved** Dini-summable error, and three geometric shields. The two-block
frame used here (`[[A,B],[B^T,-A]]`, `D = Q(x)-Q(y)`) is a weaker parallel of
that machinery; Propositions 1–2 below are special/weaker cases of results
already on record. The incremental content of this note is limited to: the
corrected order-8 completion-floor data (breadth/depth negative), the
partition observation on the known `Phi=30` objects, and the session
correction record. **Work from the RG2 note and its successors, not from
this one.**

**Status:** working framework note (2026-09-12). Contains two elementary
propositions (proved here), an exact reformulation, empirical ingredients from
this session, three named open lemmas, and the algorithm the reformulation
suggests. **No convergence claim; no proof of the doubling estimate.** Draft
for review.

## 0. Where this sits

CORE §7 reduces convergence of `alpha_n` to a Dini-summable doubling and
tripling estimate: with `H(n) = m_n^(2/3)`, it suffices that for some
`eta >= 0` with `E(N) = sum_j eta*(2^j N) -> 0`,

    H(2n) <= 2 H(n) + 2 n eta(n),   H(3n) <= 3 H(n) + 3 n eta(n)

for all large `n`. Equivalently (for `eta(n) = O(n^-delta)`):
`m_{2n} <= 2 sqrt(2) m_n (1 + O(n^-delta))`. This note maps that estimate,
via the two-block completion family, onto a finite-dimensional discrepancy
problem, proves the fragments available today, and names exactly what is
missing.

## 1. Exact reformulation

For an order-`n` signing `A` and `B in {+-1}^{n x n}` put

    K(A,B) = [[A, B], [B^T, -A]].

For states `(x, y) in {+-1}^n x {+-1}^n` (x projective),

    Q_{K(A,B)}(x,y) = Q_A(x) - Q_A(y) + x^T B y,

so with `D(x,y) := Q_A(x) - Q_A(y) in [-2 Phi(A), 2 Phi(A)]`:

    Phi(K(A,B)) = max_{(x,y)} | D(x,y) + x^T B y |.

Therefore the doubling estimate **via completions** is exactly:

    (**)  exists A with Phi(A) ~ m_n and B in {+-1}^{n x n} with
          max_{(x,y)} |D(x,y) + x^T B y| <= 2 sqrt(2) m_n (1 + o(1)).

Since `2 sqrt(2) = 1.414... > 2`, the target exceeds the trivial
no-cross bound `max|D| <= 2 m_n`: there is room for the cross term of size
up to `(sqrt(2)-1) 2 m_n ~ 0.828 m_n` on the pairs where `D` is extreme.

## 2. Two provable fragments

**Proposition 1 (cross-term floor).** For every `B in {+-1}^{n x n}`,

    max_{x,y in {+-1}^n} x^T B y = max_y ||B y||_1 >= E_y ||B y||_1
        = n E|S_n| = sqrt(2/pi) n^{3/2} (1 - O(1/n)),

where `S_n` is a Rademacher sum. In particular no completion can keep the
cross term below `sqrt(2/pi) n^{3/2} (1-o(1)) ~ 0.798 n^{3/2}` in range.

*Proof.* `max_x x^T B y = ||By||_1` (choose `x = sign(By)`); maximum is at
least the average over `y`; each row `(By)_i` is a Rademacher sum
`sum_j B_ij y_j`; `E|S_n| = n 2^{1-n} binom(n-1, floor((n-1)/2)) ~
sqrt(2n/pi)`. ∎

**Consequence (the suppression factor).** With `Phi(A) = m_n ~ c n^{3/2}`
`(c ~ 0.45)` the target forces, on the pairs with `D` near its maximum
`2 m_n`, the bound

    x^T B y <= 2 m_n (sqrt(2)-1) + slack ~ 0.37 n^{3/2},

i.e. the cross term must be **suppressed by a factor ~ 2 below its
universal floor** — but only on the high-`D` pairs, not globally. This is
why uniformly random `B` fails (its cross range `~ n sqrt(4 n ln 2) ~
1.67 n^{3/2}`, and hundreds of thousands of high-`D` pairs exceed the
bound; the measured midpoint-law overshoot `R ~ 1.3` is exactly this
excess), while structured (conference) `B` succeeds at `n = 7, 9`.

**Proposition 2 (redundancy reduction).** With `C_max = max |x^T B y|`, all
constraints `|D + x^TBy| <= T` with `D(x,y) <= T - C_max` are automatic.
Hence the binding constraints live only where

    D(x,y) >= T - C_max >= 2 sqrt(2) m_n - 0.798 n^{3/2}
              = (1.273 c^{-1} - 0.798) n^{3/2}  (c = m_n / n^{3/2}),

a strict high-`D` band. *Proof.* `|D + cross| <= |D| + C_max <= T`. ∎

So the whole doubling question, through completions, is carried by a thin
band of near-extreme pairs — a **tensor-structured discrepancy problem**:
find `B` with `x^T B y <= 2m_n(sqrt(2)-1) + slack` on the band
`X_+^delta x X_-^delta`, where `X_+^delta, X_-^delta` are the near-extremal
level sets of `Q_A`.

## 3. Empirical ingredients (this session)

- Completion exactness: `min_B = m_{2n}` at `n = 2..7`; `c`-ratios
  (min/target `2 sqrt(2) m_n`): **0.82 (n=7), 0.97 (n=9)** — algebraic
  (character-sum) blocks; search floors `1.13-1.36` elsewhere.
- **Source variance at n = 8 (WITHDRAWN — see correction below):** the
  first measurement (12 sources) used a half-space evaluator and reported
  completion floors {30 (2 sources), 32 (8), 34 (2)}. Those numbers are
  withdrawn; a corrected sweep is in progress.

**Correction (2026-09-12, same day).** The `disc2`/`deep_b` free-block state
construction fixed the second block's first coordinate (`y_0 = +1`),
evaluating a half-space (`2^{2n-2}` of the `2^{2n-1}` projective states).
Half-space maxima are **lower bounds** on the true values, which gives a
precise directional audit of everything measured with the defect:

- **Survive (non-improvement certificates):** the exhaustive diagonal-lift
  family closures at orders 16/18/20 and every "no `D`-opt improvement"
  statement. Their logic needs only `half <= true` for every tested object
  plus one independently certified true anchor (the archived full
  enumerations: 32 / 39 / 40 / 108 / 121), all of which hold and match.
- **Withdrawn (upper-bound / record claims):** the {30, 32, 34} source table
  and the claim "completion-type `Phi = 30` objects exist". The "`Phi = 28`"
  candidates never landed: independent full enumeration (correct space)
  evaluated them at 42 before any commit — the verification pipeline
  worked as intended. Corrected measurements (fixed evaluator, all values
  re-verified by full enumeration): **12-source deep sweep floors: 32 (6
  sources) / 34 (6 sources), min 32; 2,000-source quick scan: zero
  completions <= 30 found.** The completion route at `n = 8` therefore
  bottoms out at 32 > 28.28 across breadth and depth: it does not settle
  the m_16 straddle, which remains open for non-completion constructions.
- The reformulation and the two propositions are unaffected (they concern
  the exact state space).
- At `n = 16`, a Hadamard-aligned diagonal completion gives `80 <=
  2 sqrt(2) * 32` (c = 0.88): a below-target completion exists there.
- Hence: completions reach the target at `n = 7, 9, 16` and have not yet at
  `n = 8` (best 30 vs 28.28). The **m_16 straddle** — `m_16 in [27,30]`
  against `2 sqrt(2) m_8 = 28.28`, the first open doubling instance — is a
  ≤ 28-completion question for order-8 sources.

## 4. The open lemmas (tolls on the road)

- **T1 (tensor discrepancy).** For minimizer-grade `A`, there exists `B`
  with `Phi(K(A,B)) <= 2 sqrt(2) m_n (1 + o(1))`. Equivalently: a signing
  `B` near-orthogonal (with the sign and size above) to the tensor family
  `{x (x) y : (x,y) in X_+^delta x X_-^delta}` on the high-`D` band.
  *Not proved. The special tensor structure (not generic vectors) is
  essential: generic discrepancy bounds (`6 sigma sqrt(log M)`) give range
  `~ n^{3/2} sqrt(log M)` which needs `|X_+-|` polylogarithmic; the observed
  extremal multiplicities (e.g. 2,240 at order 32) are far larger.*
- **T2 (all-n uniformity).** The estimate is required for all large `n`
  (or a ratio-dense variant; see the program's alternative sufficient
  statements, e.g. the Paley-tail deficit route). Completion constructions
  are presently known only at sporadic scales.
- **T3 (barrier at 8?).** Whether every order-8 source has completion
  floor `>= 29` is open; the source-scan in flight (2,000 sources) probes
  it. A single `<= 28` completion settles the n=8 doubling instance in the
  positive; a theorem that none exists would isolate the completion route's
  first failure while leaving `m_16` itself open.

## 5. Algorithm the reformulation suggests: active-set cutting-max-plane

`disc2.py`: full exact pass `->` extract the binding pairs (top-|value|)
`->` re-optimize `B` against the active constraints (linear in `B`)
`->` repeat. Reaches `25 <= target` at `n = 7` in two iterations (the exact
optimum 21 is not reached — the method meets targets, it does not
certify optimality) and reaches the `n = 8` search floor (32) in
approximately one iteration, versus thousands of ILS rounds previously.

## 6. Repro

`~/scratch/disc2.py` (active set), `attack8_batch.py` (source-scan),
`run_discrepancy.py` (grid variant); this note's propositions are
elementary and self-contained.
