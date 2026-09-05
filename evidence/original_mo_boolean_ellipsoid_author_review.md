# Author receipt: Boolean ellipsoid shell upper

2026-09-05. The source is
`/tmp/original_mo_boolean_ellipsoid_shell_upper.md`, 322 lines, SHA-256
`ede1b62a26a636179d918ba84a48d122ab013c38175bdb9cd164bcfd8bfeb9aa`.

The author checked the complete square-completion argument, both
nonnegative remainder quantities, singular-covariance and zero-product
limits, diagonal-affine inverse inequality, full positive and negative
scalar branches, endpoint infima, two-trace simplification, and the actual
independent-cushion factor in the diagnostic profile. The root agent
reports full-read mathematical PASS, and reports that the independent
docs agent has also fully read and passed the mathematically identical
initial source with SHA-256
`fdbe8a48112a7490e64c43fe116a48b5e28ea9611209ee0548dcce0a6d16ff78`.
The final source changes only the two introductory lines describing
the penalty: the deterministic completion-square center is not called
the random constrained ellipsoid maximizer. Every displayed formula and
proof is unchanged. There are no outstanding requested source changes.

Scope: this is an exact Gaussian upper theorem and its explicit scalar
specialization. It does not prove the all-source, all-shell leading
comparison, a uniformly flat majorizer, or original convergence. An
indefinite formal reference matrix is not treated as a covariance.

## Elementary-scalar evaluation provenance

The printed decimal diagnostic in equation (21) was computed once on the
local workspace host, using the system `python` interpreter, through the
following exact command. This was a computation, not a proof step, and
was not performed on an approved remote compute host. It is recorded
honestly here rather than replaying the unchanged scalars for a different
provenance. Future computational checks must follow the repository's
offload rule.

```sh
python - <<'PY'
import math
k=2/math.pi
kg=math.pi/(2*math.log(1+math.sqrt(2)))
u=1/kg
v=-k*u
d=1-k-u
b=d+k*v
Delta=(1-u)*(2*k-1+u)*(1-v*v)
f=(k-d*v+math.sqrt(Delta))/2
eta=-b/(k+d*v+math.sqrt(Delta))
print('kappa',repr(k),'K_G',repr(kg),'u',repr(u),'v',repr(v))
print('B',repr(b),'Delta',repr(Delta),'eta',repr(eta),'factor_squared',repr(f),'factor',repr(math.sqrt(f)),'width_coefficient',repr(2*math.sqrt(f)))
PY
```

Its output was

```text
kappa 0.6366197723675814 K_G 1.7822139781913693 u 0.56109985233918 v -0.35720726027165234
B -0.42512482942894814 Delta 0.3194667143878699 eta 0.334096674476977 factor_squared 0.5656033764926796 factor 0.7520660718930748 width_coefficient 1.5041321437861497
```

The theorem depends on the exact displayed expressions and their analytic
proofs, not on these rounded decimal values. No source-matrix enumeration,
solver run, simulation, or statistical inference was used for this note.
