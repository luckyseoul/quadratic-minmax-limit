# Doubling debt across the proved window

**Status:** checked arithmetic, not a proof of the miss. The limit stays open.
No branch of this repo records these window-endpoint figures. Code search on `main` finds neither `0.336493` nor `0.519124`. The eleven other branches (`archive/*`, `codex/leftover-moment-attack`, `maxplus-p11-enumeration`, `navier-stokes-techniques`, `prop15586-maxplus-gram-reduction`, `research/paired-disagreement-floor-20260924`, `residual/p13-u6-common-moments`, `strategy/2026-08-28-unattempted-directions`) are earlier attacks; the paired-disagreement tip only touches `HANDOFF.md`. This note does not repeat Propositions 1–2 of `evidence/NOTE_2026-09-12_COMPLETION_DISCREPANCY_FRAMEWORK.md`. It replaces that note's sample ratio `c ~ 0.45` by the proved window.

## Setup

Let `F(n) = m_n` and `c_n = F(n) / n^{3/2}`. The field-plus-spin theorem and the Paley cap give, for every large `n`,

    gamma_* <= c_n <= 1/2,
    gamma_* = 2 phi(t_*) (2 Phi(t_*) - 1) = 0.336493364432... ,

where `t_* > 0` solves `2 phi(t) = t (2 Phi(t) - 1)`. For equal blocks the bridge satisfies `R(n,n) >= (sqrt(2/pi) - o(1)) n^{3/2}`, and `sqrt(2/pi) = 0.797884560803...`. The exact block identity is

    Phi(K) = max ( |Q_A(x) + sigma Q_B(y)| + |x^T W y| ),

so the sum of the separate maxima is only an upper estimate:

    Phi(K) <= 2 F(n) + R(n,n).

The doubling target is `2^{3/2} F(n) = 2 sqrt(2) c_n n^{3/2}`.

## Two different gaps

**Sum-of-maxima gap.** The estimate `2 F + R` exceeds the target by

    ( sqrt(2/pi) - 2 (sqrt(2) - 1) c_n ) n^{3/2}.

At the floor `c = gamma_*` this is `0.519124330410... n^{3/2}`: the joint maximum has to fall to `2 sqrt(2) gamma_* n^{3/2} = 0.951746959257... n^{3/2}`, and `2 gamma_* + sqrt(2/pi) = 1.470871289667... n^{3/2}`. At the conference end `c = 1/2` the same gap is still `0.383670998430... n^{3/2}`.

**Extreme-pair budget.** On a pair with `|Q_A(x) - Q_A(y)| = 2 F(n)`, the identity forces the cross term itself to be at most `2 (sqrt(2) - 1) F(n)`. Inside the window that budget is at most `(sqrt(2) - 1) n^{3/2} = 0.414213562373... n^{3/2}`, so the debt against the bridge floor is at least

    sqrt(2/pi) - (sqrt(2) - 1) = 0.383670998430... n^{3/2}

even at the conference end. The September 12 note's sample `~ 0.37 n^{3/2}` was this budget at the guessed ratio `c ~ 0.45`, not either endpoint. Its displayed band identity `(1.273 c^{-1} - 0.798) n^{3/2}` is not the expansion of `2 sqrt(2) c - sqrt(2/pi)`; the coefficient in `c` is `2 sqrt(2) c - sqrt(2/pi)`, and `1.273` is just `2 sqrt(2) * 0.45`.

## Why the bridge alone does not break `H`

Put `H(n) = F(n)^{2/3}`. At the floor,

    H(n) + H(n) = 0.967566890422... n,
    R^{2/3} = 0.860254013828... n,
    (2 F)^{2/3} = 0.767958349853... n,
    (2 F + R)^{2/3} = 1.293351131286... n.

Concavity lifts `H` by `0.092295663976... n` when the bridge replaces the two blocks, which still sits under `H+H`. Charging the bridge on top of the blocks sits `0.325784240864... n` over `H+H`. That overshoot is exactly the sum-of-maxima gap above; it is not a non-summable `eta` until the joint maximum is shown to miss `2 F + R` by `0.519... n^{3/2}` at the floor, and then by a further Dini-summable tail. A best-response fixed point has no such miss. This note does not produce one.
