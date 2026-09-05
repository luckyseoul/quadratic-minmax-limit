# Quadratic minimax limit

Does the sequence

`alpha_n = n^(-3/2) min_{a_ij in {+1,-1}} max_{x in {+1,-1}^n} |sum_{i<j} a_ij x_i x_j|`

converge as `n` tends to infinity?

The original MO limit is OPEN. A proof of existence need not identify its
value. The proved bounds are `1/pi <= liminf alpha_n <= limsup alpha_n <= 1/2`.

Start with [CORE.md](CORE.md) for the route-neutral mathematics,
[STATUS.md](STATUS.md) for current claims, and [HANDOFF.md](HANDOFF.md) for
continuation. [ARTIFACTS.md](ARTIFACTS.md) records the 2026-09-05 reset and
the complete archived research states.

The long [solution.md](solution.md), proof notes under `evidence/`, and
theorem checkers under `src/` retain scoped mathematical results. They are
not a complete solution. Historical Paley, residual-(ii), and amplification
ledgers are optional research routes, not global acceptance criteria.

Global status is tracked in `src/original_mo_status.py`; the old E1 status
interface is retained only with separate route-local diagnostics.
