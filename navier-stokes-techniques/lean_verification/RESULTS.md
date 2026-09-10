# Verification results (soulkiller, 2026-09-10)

Machine: 88 cores, 60 GB RAM, Tesla V100 (unused — this is a CPU workload).
Repo: `openai/NavierStokesAndEuler` @ `main`, Lean `v4.34.0-rc2`.

## Full build

`lake exe cache get && lake build` — **succeeded, 11,251 jobs.**

Lean's own kernel reports the top-level theorems' axiom dependencies during
the build:

```
Euler/Solution.lean:73:0: 'Euler.euler_breakdown_R3' depends on axioms:
  [propext, Classical.choice, Quot.sound]
Euler/Solution.lean:75:0: 'Euler.exists_compact_smooth_euler_singularity' depends on axioms:
  [propext, Classical.choice, Quot.sound]
NavierStokes/ComparatorSolution.lean:31:0: 'NavierStokes.Comparator.navier_stokes_breakdown_R3' depends on axioms:
  [propext, Classical.choice, Quot.sound]
NavierStokes/ComparatorSolution.lean:32:0: 'NavierStokes.Comparator.navier_stokes_breakdown_periodic' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

Only the three standard classical axioms. **No `sorryAx`**, no extra axioms.

## Comparator, Navier–Stokes challenge

```
lake exe comparator ComparatorChallenges/NavierStokes.json
```

```
Running nanoda kernel on solution
nanoda kernel accepts the solution
Running Lean default kernel on solution.
Lean default kernel accepts the solution
Your solution is okay!

real    8m39.374s
user    25m47.476s
sys     0m25.493s
COMPARATOR_NS_EXIT_0
```

Both theorems checked: `navier_stokes_breakdown_R3` (Clay alternative (C))
and `navier_stokes_breakdown_periodic` (alternative (D)), against
`permitted_axioms = [propext, Quot.sound, Classical.choice]`.

## On the `sorry` warnings

The build emits:

```
warning: ComparatorChallenges/NavierStokes.lean:273:8: declaration uses `sorry`
warning: ComparatorChallenges/NavierStokes.lean:280:8: declaration uses `sorry`
```

These are **not** in the proof. `ComparatorChallenges/NavierStokes.lean` is the
*challenge* file: it only **states** the two theorems (adapted from DeepMind's
formal-conjectures formalization of the Clay problem statement), with `sorry`
standing in for proofs. The actual proofs are in the separate `solution_module`,
`NavierStokes.ComparatorSolution`. Comparator's job is precisely to check that
the solution proves *those* statements — and since `sorryAx` is absent from
`permitted_axioms`, a vacuous solution would be rejected.

## Parallelism

Stock Comparator runs nanoda single-threaded (it omits `num_threads` from the
config it generates; nanoda defaults that field to `0`, taking the serial
branch). With the one-line patch documented in `README.md`:

* observed `nanoda_bin` at **nlwp = 89** (88 checker threads + main), peaking
  around **3488 % CPU**;
* nanoda's own phase completed in roughly 1–2 minutes.

Note the first ~10 s of nanoda's run is single-threaded — that is the serial
export-file *parse* phase; worker threads spawn only once parsing finishes.
`lean4export` is likewise serial and runs twice (challenge, then solution).
So the 8m39s wall time is dominated by the serial export/replay phases, not by
the kernel check; the patch shortens only the last part.

## Comparator, Euler challenge

```
lake exe comparator ComparatorChallenges/Euler.json
```

```
Running nanoda kernel on solution
nanoda kernel accepts the solution
Running Lean default kernel on solution.
Lean default kernel accepts the solution
Your solution is okay!

real    13m3.484s
user    38m13.031s
sys     0m22.826s
COMPARATOR_EULER_EXIT_0
```

Both theorems checked: `euler_breakdown_R3` and `exists_compact_smooth_euler_singularity`
(finite-time singularity for smooth, compactly supported, divergence-free
Euler initial data), against the same axiom allowlist. Same `sorry`-in-the-
challenge-file caveat as above applies (`ComparatorChallenges/Euler.lean:85,170`);
the proofs are in `Euler.Solution`.

Both challenges: **passed, independently, on both kernels.**
