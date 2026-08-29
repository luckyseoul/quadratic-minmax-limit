# p=19 binary affine-Radon reduction

**Date:** 2026-08-29
**Proposition:** 15.692
**Status:** proved; sharpens the fourteen-profile endpoint but does not close it

Let `A` be the binary line-point incidence matrix of `AG(2,p)`, with its
rows grouped into the `p+1` parallel classes. For odd `p`,

```text
A^T A = I + J  over F_2.
```

The diagonal is `p+1=0 mod 2`, while two distinct affine points determine
one line. Therefore `A^T A` is the identity on even-weight point vectors.
Moreover, each directional block of `Ax` has parity `wt(x)`. Thus `A`
maps the even point space into the direct sum of the even subspaces of the
direction blocks. Both dimensions are

```text
p^2-1 = (p+1)(p-1).
```

Consequently this restriction is an isomorphism and its inverse is

```text
x = A^T r.
```

For the fourteen `p=19` profiles left by Proposition 15.689, choose `b_d`
lines in direction `d`, where the twenty block weights have the prescribed
phase histograms. Every `b_d` is even, so there are no further linear
compatibility conditions. The exact remaining condition is the nonlinear
inverse-weight equation

```text
wt(A^T r) = 16.
```

The first two moments do not supply the missing obstruction. If `N(v)` is
the number of chosen lines through `v`, then distinct directional stripes
are pairwise independent. Each profile fixes `E[N]` and `E[N(N-1)]`, but
both moments are matched by an explicit nonnegative distribution supported
entirely on `{4,6,8}`. Hence a second-moment argument cannot force even one
point with odd `N(v)`. The exact mod-four inverse-weight congruence also
accepts all fourteen profiles.

This rules out linear Radon compatibility, mod-four weight, and pairwise
moment bounds as closure mechanisms. A successful finite attack must use
the exact weight-16 inverse condition or genuinely higher concurrency
information. The endpoint, residual (ii), R1, Type I, and the limit remain
open.

Reproduction:

- `src/e1_gmin_m4_prop15692.py`
- `evidence/e1_gmin_m4_prop15692.json`
- `tests/test_prop15692.py`
- `scripts/p19_second_boundary_profile_cryptominisat.py` (exact native-XOR
  model; a timeout is not treated as evidence)
