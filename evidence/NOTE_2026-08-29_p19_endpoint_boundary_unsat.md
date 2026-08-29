# p=19 second all-finite endpoint exclusion

**Date:** 2026-08-29
**Proposition:** 15.699
**Status:** proved computationally; p=19 endpoint closed

After Proposition 15.698, three profiles remained:

```text
slack 24: phase zero {0:5,16:5}, phase one {2:9,6:1}
slack 28: phase zero {0:5,16:5}, phase one {2:9,10:1}
slack 32: phase zero {0:5,16:5}, phase one {2:9,14:1}
```

For each profile, the exact native-XOR model imposes a 16-point boundary,
all 380 affine line parities, `r=A*x`, the inverse `x=A^T*r`, and the exact
directional weight histogram. No edge-lift variable, quadratic floor, or
repair assumption is used. A phase-zero `b=0` line supplies the lossless
pair normalization to field elements zero and one.

Completed CryptoMiniSat runs returned `UNSATISFIABLE`:

```text
slack 24: nuka 101.79 s; soulkiller ECC 121.83 s
slack 28: soulkiller ECC 79.23 s
slack 32: jellyfin 92.92 s; soulkiller ECC 98.28 s
```

All raw JSON files and hashes are in
`evidence/p19_endpoint_boundary_unsat/`. Orin's separate slack-32 run
returned `UNKNOWN` and is not archived or used.

Nonsquare dilation flips the direction type and `c_H` together while
preserving the phase profile, so the exclusion covers both signs. Therefore
the p=19 second all-finite endpoint is closed. Residual (ii), later boundary
sizes, p=17 and p=23 endpoints, R1, Type I, and the limit remain open.
