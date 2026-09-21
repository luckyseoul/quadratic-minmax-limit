# Actual Gibbs conditional-score check

This supports the local algebra in
`../NOTE_2026-09-21_ACTUAL_GIBBS_COVARIANCE_FLOOR.md`.
The proof is author-reviewed, not independently human-reviewed or formally
verified. The original MO limit and cross-order pressure comparison are open.

NUKA was selected for one small serial SymPy check, not a parallel signing
search. Preflight on 2026-09-21 at 09:11:21 UTC confirmed SSH reachability,
16 logical CPUs, load averages 0.00/0.00/0.00, and SymPy 1.14.0. No other
node was needed. The mesh-dispatch procedure determined this allocation.

Fresh staging: `/tmp/mo-score-floor-20260921.Gl3z06Ul`.
The command in that directory was:

```sh
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 check_scores.py
```

Result: exit 0, twelve symbolic identities passed. `result.json` preserves
the JSON emitted by that run. The checker used four two-spin configurations
with symbolic positive Boltzmann weights, so no numerical tolerance or
finite matrix census entered the identities. Positivity of the parameters
allows arbitrary real fields and coupling through their exponentials.

SHA-256 of the tested checker (matched on controller and NUKA):
`3c8b2d5aa6664765e0263dac06219e1f344d0622fa42062e0a6240ff692ab4c7`.
SHA-256 of the final author-reviewed proof:
`b540a85b9a210cbb62226f58f239d5856713c3b3e61ee527fb2d7436c3aaeb38`.
SHA-256 of the retained result JSON:
`748fd9ed8a81c9ea61c12c6e3dcc2f2f8ff21db730a100f50c36dac21e78fc22`.
The proof was reviewed after the algebra run; an already-known scalar
row-response corollary was removed, without changing the tested checker.

What the run checks: conditional centering, spin/score covariances, the
exact two-spin score identity and its symmetry, the derivative and scalar
integrals used in the Lipschitz bound, and edge-flip threshold algebra.
The matrix Schur bound, Cauchy--Schwarz argument, optimality implications,
and cavity change-of-measure bound rest on the written analytic proof;
the checker does not formally verify them or prove convergence.

Deduplication checked the active notes and solution, preserved branch
heads, linked worktrees, and the indexed campaign's preserved frontier.
The old field-response theorem already supplies a scalar removed-row
response. The new claim is the simultaneous covariance matrix floor,
including its actual nonzero-field cavity version.

This is a supporting analytic lemma, not a major research or structural
milestone and not a destructive edit. No new large-drive backup was made;
prior verified checkpoints were preserved. No remote research job remains
from this foreground run.
