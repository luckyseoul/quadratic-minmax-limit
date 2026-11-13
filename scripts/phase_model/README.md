# Phase model of the Paley conference form

Scripts behind `evidence/NOTE_2026-11-13_PALEY_PHASE_MODEL.md`.

- `phase_lib.py` — split-torus coordinates, Jacobi sums `J_k`, square-wave
  states, the closed formula `rho*(k, delta)`, and exact integer evaluation of
  `Q` for the Paley conference matrix (Kronecker substitution, no floating
  point).
- `verify.py [q ...]` — Lemma 1, Lemma 2, Proposition 3 and Theorem R against
  exact integers. Output: `evidence/phase_model_20261113/identities.json`.
- `scan.py q ...` — best square-wave state per prime, exact `Q`, deterministic
  1-opt polish, Weyl sums and box counts for the phases at `k, 3k, 5k`.
  One JSON record per line; `evidence/phase_model_20261113/scan.json` is the
  list for the thirteen primes of the note.

Depends on `scripts/parity_ceiling/paley_lib.py`. Python 3 with numpy only.
