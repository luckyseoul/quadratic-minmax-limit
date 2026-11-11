# Order 50: both coherent lifts give 169, so that family cannot reach 167

2026-11-11. Follow-up to
`NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md` (§5, prediction 4). The
original MO limit is OPEN. Nothing here bounds `m_50` below 169.

## Result

The order-50 coherent lift of the September campaign's `B25`

    K = [[B, B + d I], [B + d I, -B]],   d = +-1,

has been evaluated exactly on both signs over the full `2^49` projective cube:

    d = +1 :  Phi = 169   (8311.102 s, evidence/k50_nuka_gpu_20260919/)
    d = -1 :  Phi = 169   (8245.568 s, this note)

Same source matrix (md5 `4877cb1cb3ed29b09e91d1dcab1372fd`), same tiling
(`Px=256`, `Py=1048576`), same 9070 XT / hipBLAS backend, both runs reaching
`1/1` quadrant completion. So `Phi(C50) = 175`, the plus-I lift is `169`, and
the minus-I lift is also `169`.

## What this decides

Prediction 4 is "some signing of order 50 has `Phi <= 167`". Both signs of the
only order-50 lift construction on file give `169`, so **this family is
exhausted without meeting the prediction**. The gap from the prediction to the
family is 2, and the gap from Paley `C50` to the family is 6.

That is the same shape as order 36 in
`NOTE_2026-11-11_COHERENT_LIFT_FRUSTRATION.md`: the coherent lift is a
fixed-point structure and the law's target is below its floor. If 167 is
reachable at order 50 it must come from outside this family, which no longer
has untested signs.

`alpha = 169 / 50^(3/2) = 0.47800`, unchanged from the plus-I receipt.

## Reproduction

    hipcc --offload-arch=gfx1201 -O3 scripts/k50_plusi_hip_phi.hip.cpp \
        -L/opt/rocm/core-10.0/lib -Wl,-rpath,/opt/rocm/core-10.0/lib \
        -lhipblas -o k50v
    ./k50v B25.txt 256 1048576 -1

Receipt: `evidence/k50_lift_20261111/minus_i_result.json`. The source patch
(4th argument `dsign`) is the one described in the HANDOFF entry for
2026-11-11; the September build is unchanged.
