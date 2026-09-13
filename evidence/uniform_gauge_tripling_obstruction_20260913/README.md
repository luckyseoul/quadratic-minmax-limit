# Uniform gauge tripling obstruction: exact corroboration

The all-orders proof is
`../NOTE_2026-09-13_UNIFORM_GAUGE_TRIPLING_OBSTRUCTION.md`. It excludes the
recently proposed uniform reverse-KL tripling criterion for all sufficiently
large orders, using the existing CORE upper bound. It does not exclude
one good gauge or prove convergence of the original minimum.

One fixed `4+4` input tests the proof's randomized Boolean strategy. The
fixture, including both internal matrices, the cross matrix, and the seed
set, is recorded in `result.json`. This is not a source-signing census.

NUKA reported 16 hardware threads and zero load averages before dispatch.
One serial standard-library process ran in a fresh directory; its staged
source hash was verified before execution:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /usr/bin/python3 /tmp/qml-gauge-average-floor.96Po3OmA/check_strategy.py
```

Exit status was zero. The recorded computational time was about 0.016
seconds. All 64 conditional contexts had exact quadratic self-interaction
cancellation. In 24 contexts that interaction was nonzero before the block
sign was averaged. All 832 strategy outcomes satisfied the actual gauge
diamond lower bound, and the exact mean score was `43/4`, matching the
analytic formula. Independent fair coins at zero fields are included.

Rational arithmetic also checked the polynomial identity giving the
elementary bound `pi<22/7`, the square-root bounds, and the strict constant
gap `3795009/354580000>1/100`. No floating-point constant establishes the
gap. The asymptotic inference uses the analytic proof and CORE, not finite
enumeration. Root reviewed that proof; this check is not an independent
human review.

SHA-256:

- `check_strategy.py`: `b3cde74b361e881958f5239f3aa4c899c759d58fec326228fd61ffcfcc14202e`
- `result.json`: `850c0cfe49ed64614a7b586dc4b8f4d9939f1d774286737527a2fd3638dc84ab`

The earlier conditional derivation is retained and marked as an impossible
hypothesis, rather than erased. This focused correction requires no fresh
major-milestone backup. The original convergence problem remains OPEN.
