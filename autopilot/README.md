# Autonomous research controller

This directory defines a cruise-control research loop for the original MO
quadratic-minmax limit problem. It is deliberately route-neutral: the goal is
still only to prove convergence or nonconvergence. Finite searches are evidence,
not completion.

## Operating model

The controller owns a persistent hypothesis graph and launches short, isolated
worker jobs. Workers communicate only through JSON result artifacts. A result is
never promoted to a theorem merely because an agent says it is true.

Each task has:

- one precise implication to establish;
- a route and parent hypothesis;
- an experiment/proof budget;
- an independent destroyer task when it reports a positive result;
- a deterministic result schema and content hash.

The controller scores routes by information gained per compute-hour, not by how
many matrices were searched. Duplicate routes are suppressed by canonical task
hashes. Repeated failures lower a route's allocation. New positive evidence
spawns generalization and adversarial tasks automatically.

## Default portfolio

1. scaling/composition: attack inequalities for `H(n)=m_n^(2/3)`;
2. restriction/lifting: turn observed good finite matrices into unconditional
   all-orders comparisons;
3. optimizer structure: mine spectra, switching invariants, row/cut moments,
   active states, and principal restrictions;
4. spectral/rounding: search for unconditional source inequalities;
5. counterexample: actively seek separated subsequential limits;
6. proof synthesis: convert recurring numerical invariants into quantified lemmas;
7. adversarial review: try to break every promising claim;
8. exploration: reserve a small fraction for genuinely unrelated approaches.

The allocation is adaptive. A route that stops producing new implications is
starved automatically; a route producing independently verified progress gets
more slots.

## Cruise-control rules

- No human input is required for ordinary execution, retries, pivots, or
  archiving of routine artifacts.
- The controller may stop, restart, resize, or replace jobs.
- It may request approval only for a genuinely consequential external action
  (for example, publishing a claimed proof, destructive filesystem changes,
  or an unavailable credential).
- It must never edit or overwrite reviewed mathematical artifacts in place.
- A major milestone uses the existing `scripts/milestone_backup.sh`; routine
  runs do not create backup churn.
- Compute is checkpointed frequently so a host loss only loses the active task.
- The controller fails closed on malformed results, missing provenance, stale
  inputs, or contradictory theorem claims.

## Worker contract

A worker receives a task JSON plus read-only copies/references to `CORE.md`,
`STATUS.md`, `HANDOFF.md`, `ARTIFACTS.md`, and relevant evidence. It must emit
one JSON result containing:

`status`: `success`, `counterexample`, `partial`, `dead_end`, or `error`;
`claim`: one sentence;
`implication`: the exact mathematical statement established or tested;
`evidence`: artifact paths and hashes;
`scope`: finite, asymptotic, conditional, or all-orders;
`dependencies`: prior claims used;
`attack_targets`: assumptions a destroyer should test;
`next_tasks`: bounded follow-ups;
`compute`: wall time, host, and workload counts.

Free-form prose is retained as a human-readable sidecar, but JSON is the
machine control plane.

## Launch

The reference implementation is `scripts/research_autopilot.py`. It is an
orchestrator, not a model provider: set `AUTOPILOT_WORKER_CMD` to the existing
local runner/agent launcher. The command receives the task JSON path as its
first argument and must write the result JSON path supplied by the environment.

A node may therefore use an existing runner without changing the mathematics.
SSH hosts, local CPU/GPU slots, and per-host concurrency live in
`autopilot/config.json`.

The first run should be a dry run. After checking the generated task portfolio,
launch normally. No broad unchanged permutation census is seeded by default.
