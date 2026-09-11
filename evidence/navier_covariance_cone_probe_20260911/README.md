# Navier-style covariance-cone scout

**Status: finite floating-point diagnostic; no all-orders conclusion.**

The useful source-side transfer from the Navier--Stokes construction is its
strictly positive covariance cone.  The active all-law source proof uses the
two correlation phases induced by `|M|+M` and `|M|-M`, where
`M=A/sqrt(n)`.  This scout tested the natural extension

```text
P + M >= 0,   P - M >= 0,
```

using a common PSD majorizer `P` selected by independently weighted trace
minimizations.  For each candidate, it compares the diagonal-only source
baseline

```text
sum_(i<j) 1/sqrt(P_ii P_jj)
```

against the `P=|M|` baseline.  This is the exact diagonal quantity that
appears before the arcsine lower bound in the existing two-phase argument;
the optimization and comparison are floating point, so this is a scout only.

Soulkiller ran 88 independent random complete signings at order 10, with 12
log-uniform positive weight rays per signing, using CLARABEL SDP solves.  No
weighted majorizer improved the `|M|` score:

```text
instances                  88
positive improvements       0
best relative change       -0.01581048839720356
worst relative change      -0.06416032905042823
concatenated receipt SHA-256 d457cce73eb2b0b637c7c5452f38d9c8081e139797663b73f0369197cc4985a6
```

Each receipt is `SEED.json`.  The scoped driver is
`scripts/navier_covariance_cone_probe.py`.  A negative result here does not
prove that `|M|` is optimal over the full cone: only 12 weighted rays and a
finite sample were tested.  It does rule out treating this obvious
weighted-majorizer variant as a promising accelerator without a new
all-orders cone theorem or a different source objective.

