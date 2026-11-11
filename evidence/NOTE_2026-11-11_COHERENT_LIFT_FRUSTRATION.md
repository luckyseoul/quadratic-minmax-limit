# The order-36 plateau at 108 is a matching-frustration cost

2026-11-11. Follow-up to
`NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md` (§5, prediction 2) and to
`NOTE_2026-11-11_C42_FALSIFIES_ROM_LAW.md`. The original MO limit is OPEN.
Nothing here bounds any `m_n`. The new material is an exact decomposition of
one shipped lift value, and the reason its natural relaxation cannot be
attained by any signing of that family.

## Result

Let `A` be a symmetric signing of order `m` and let

    K0 = [[A, A], [A, -A]],   order 2m

be the coherent lift. Its `2m` cross-block diagonal entries are the zero
diagonal of `A`, so `K0` is not a signing: each of those entries must be
replaced by `+1` or `-1`. A completion is `B = A + diag(d)` with
`d in {+-1}^m`, and these are *all* completions of `K0`, because the
cross-block diagonal is the only off-main-diagonal zero of `K0`.

For the Paley sources `q = 5, 13, 17`:

    q    order   Phi(K0)   Phi(B=A+I)   Phi(B=A-I)   relax attained
    5     12       18         18           18            yes
    13    28       70         70           70            yes
    17    36      102        108          108            NO

Every value is exact over the full projective cube: `x` runs over
`2^(m-1)` representatives (the form is invariant under
`(x, y) -> (-x, -y)`) and `y` over all `2^m`, so `2^(2m-1)` states —
`2^35` at order 36.

At order 36, `Phi(K0) = 102` is exactly the number prediction 2 of the
source note expects (`S = 108`, Haar `0.9579 * S = 103.4`, "expected 102").
The law's value is therefore realized, but by the *relaxation*: the
cross-block zeros are not signs, and no signing of this family reaches it.

## Mechanism: the extremal matching patterns are not co-orientable

At a state `(x, y)`,

    Q_K(x, y) = Q_K0(x, y) + d . z,     z = x o y in {+-1}^m,

so a completion shifts each state by a signed sum over the `m` matching
entries. The relax-extremal states and their matching patterns:

    q    distinct patterns at relax max   minimax over extremal set   floor
    5     10                               16  (<= 18)                  18
    13    78                               68  (<= 70)                  70
    17   204                              106  (> 102)                 108

At `q = 5` and `q = 13` one diagonal drives every extremal pattern's
contribution negative, so the completion loses nothing and the lift sits at
its relaxation value. At `q = 17` no diagonal does: the 204 patterns are
mutually incompatible, and the extremal set alone forces every completion to
106. The floor over all `2^18` diagonals is 108. That family was already
exhausted by the September campaign's wave-6 sweep: 32 order-18 sources x
262144 diagonals = 8,388,608 completions, global minimum **108**, 11 sources
attaining it, per-source minima in `{108, 110, 112, 114, 116}`
(`evidence/original_mo_broad_campaign_20260906/v100_campaign.tar.gz`,
`stage/diagonal_pool_wave6_results/`). The two spot values `B = A +- I`
above agree.

So 108 is not a search failure. It is the price of the sign alphabet on 18
entries: the relaxation that reaches the law's 102 is frustrated, and the
frustration costs 6.

## What this says about prediction 2

"Some signing of order 36 has `Phi <= 104`" is neither refuted nor proved
here. What changes is the reading of the plateau:

- The source note's remark that single-edge moves "cannot leave a crystalline
  maximiser set" now has a cause. The escapes are blocked by a matching
  obstruction, not by the search.
- If 104 is reachable at order 36, it must come from a signing *outside* the
  coherent lift family: every member of that family is `>= 108`, and its
  relaxation floor of 102 is unreachable by construction.
- The law's table has nine exact matches and four one-step-low rows. The
  single large gap is order 36 (`108` against about `103.4`), and it is
  exactly the order at which this coherent lift is frustrated. That gap is
  structural, not noise.

The wider reading is that the law's central number is right about the relaxed
problem. The sign alphabet, not the random-orthogonal prediction, supplies the
last six units at order 36.

## Reproduction

    python3 scripts/rom_law/lift_frustration.py 5 13 17

CuPy 14.1.0 on the Tesla V100; the order-36 pass is two complete `2^35`
enumerations (about 3 s each) plus an exhaustive `2^18` search for a common
orientation. Receipt:
`evidence/coherent_lift_frustration_20261111/result.json`.

Prior values this note relies on and does not recompute: `Phi(C18) = 33`,
`Phi(C14) = 21` (`scripts/paley_structural_lift_family.py`), the order-28
lift value 70, and the wave-6 diagonal exhaustion above. The 108 at order 36
is independently reproduced here in two ways (the `B = A +- I` passes and the
`lift_opt.py` anchor).
