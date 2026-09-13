# Positive-temperature gauge certificate: scope and exact check

The analytic error criterion is in Section 4.1 of
`../NOTE_2026-09-02_UPSTREAM_RELATIVE_GAUGE_BRIDGE.md`.

It combines the existing pressure/norm sandwich with the gauge partition
normalization. At `beta=[log(n+1)]^2/sqrt(N)`, a one-sided `O(N)` error
in the full log density gives a Dini-summable normalized norm error.
This is a conditional certificate: the leading model inequality and the
required uniform remainder bound for the actual blocks remain unproved.
It corrects a claimed universal need for integer-level occupancy accuracy.
No novelty is claimed for the underlying soft maximum inequality.

Root derived and reviewed the analytic argument. The computation below is
finite corroboration of its normalization, not an independent human review
or a proof of the all-orders hypotheses.

NUKA was reachable, reported 16 hardware threads and zero load averages,
and ran one serial standard-library process in a fresh directory:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /usr/bin/python3 /tmp/qml-gauge-temperature.yskm7IW8/check_partition.py
```

The staged source hash was verified before execution. One fixed `3+3`
source signing and its 32 projective gauges were checked. Integer Laurent
coefficients prove that the averaged partition polynomial is the product
of the three source partition polynomials for this fixture, at every
positive temperature. Omitting the `tau` average fails that same identity.
Exact rational evaluations at `exp(beta)=2` corroborate the norm sandwich
and its `2^N` normalization for all 32 gauges. Exit status was zero; the
recorded computational time was about 0.007 seconds.

This is a new check of the temperature normalization. It does not repeat
the previous `H4/H6/H8` collision search, evaluate a source-signing census,
or check the leading model inequalities (11).

SHA-256:

- `check_partition.py`: `0e6b7c021286110d47e031ff48669d822a23ef1872c8b167ab05926a29f1eb73`
- `result.json`: `70a35b1ff14cb26105ff571cc56d557b0e0fbf82b0a87ac3c8bd964615e98d50`

No fresh large-drive backup was required for this supporting proof-scope
correction. The original convergence problem remains OPEN.
