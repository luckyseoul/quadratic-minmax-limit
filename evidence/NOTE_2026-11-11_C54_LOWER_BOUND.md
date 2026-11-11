# Φ(C54) ≥ 189: the law's order-54 prediction is refuted

2026-11-11. Follow-up to
`NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md` (§5, prediction 5). The
original MO limit is OPEN. Nothing here bounds any `m_n`; the result is one
explicit signing at an order that had no value on file.

## Result

Let `C54` be the Paley conference matrix of order 54 (`q = 53`), checked
entrywise: zero diagonal, `+-1` off the diagonal, `C54 C54^T = 53 I` exactly.
There is a Boolean `x` with

    (1/2) x^T C54 x = 189,

verified twice in exact integer arithmetic (int64 matmul and a pure-Python
double loop over all 2916 ordered pairs). Witness bitmask of the negative
coordinates: `10774371055909017`.

Hence

    189 <= Phi(C54) <= 196.

The upper end is the flat-spectrum ceiling `54 sqrt(53) / 2 = 196.5630` with
integrality. `x^T C54 x` is a sum of `54 * 53 / 2 = 1431` terms of `+-1`, so
every admissible value is odd; with `Phi <= 196` the bracket is

    Phi(C54) in {189, 191, 193, 195}.

Prediction 5 of the source note is `187` or `185` (`S = 196.56`,
`Haar * S = 188.5`, "the admissible integer below that, or one step lower").
**189 already exceeds both alternatives**, and it also exceeds the Haar mean
`188.5` itself, so the prediction fails at the level of its central estimate
and not only at its parity step.

## Search

`scripts/rom_law/c54_search.py`: batched 1-opt ascent on the halved form,
`2` signs x `16` rounds x `131072` random restarts, plus `12` kick rounds of
`131072` three-coordinate perturbations of the incumbent. The value 189 was
found in the first round of 131072 restarts, confirmed in round 0 of the
`-C54` pass, and never improved across the remaining ~7.2M local optima
(4.19M random restarts plus 3.15M kicked ones). The
search produces witnesses only: 189 is a lower bound, and no exact maximum is
claimed. Closing the remaining `{191, 193, 195}` needs either an exact
`2^53`-state enumeration (about 29 core-days at the measured `gs2s` rate) or a
branch-and-bound with a tighter bound than the spectral ceiling; neither was
run.

Numpy on the CPU was used rather than the V100: CuPy 14.1.0 on this host
fails NVRTC compilation for 2-D `max`/`argmax` reductions (the documented
V100 header problem, now also at `cuda_fp4.hpp`), and the workload is small
enough that BLAS on the host is not the bottleneck.

## Where this sits

`alpha = Phi / n^(3/2) = 189 / 396.817 = 0.47629`, against `0.47393` for C42
(129), `0.47800` for the order-50 plus-I lift (169) and `0.49504` for C50
itself (175). The new point continues the flat ~0.476 trend at an order not
previously on file and stays below the order-50 lift.

## Reproduction

    python3 scripts/rom_law/c54_search.py --rounds 16 --kick-rounds 12 \
        --out evidence/c54_20261111/witness.json

Receipt: `evidence/c54_20261111/witness.json` (witness state, exact integer
value, ceiling, and the prediction it decides).
