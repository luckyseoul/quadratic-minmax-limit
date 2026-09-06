# Broad CPU/GPU campaign, 2026-09-06

The original MO convergence problem remains **OPEN**. This is finite
computational evidence, not a proof of convergence, nonconvergence, a limit
value, or optimality of the sampled matrices. The reviewed asymptotic bounds
remain `1/pi < liminf alpha_n <= limsup alpha_n <= 1/2`.

The campaign began from main commit
`6bcd7d46f307d0e55141244daa3e58044f32cf55` after Nick requested a broad
multi-GPU/many-core search. It used the Tesla V100 on soulkiller, the
RX 9070 XT on Nuka, and 100 CPU search workers allocated 76/12/12 across
soulkiller/Jellyfin/Nuka. Independent integer field-update verification ran
on Jellyfin. Orin was not used for a workload already supported by the V100.

## What the finite results say

The explicit winner bank gives these UPPER bounds on m_n, in matching
increasing-order columns. Each is the exact norm of a retained matrix,
not a proof that the bound is optimal.

| Orders | Upper bounds on m_n, in the same order |
| --- | --- |
| 11, 12, 13, 14, 15, 16, 17, 18 | 17, 18, 20, 21, 29, 32, 32, 33 |
| 19, 20, 21, 22, 23, 24, 25, 26 | 41, 42, 44, 49, 53, 56, 60, 61 |
| 27, 28, 29, 30, 31, 32 | 67, 70, 72, 75, 79, 80 |

In particular, one order-26 matrix has norm 61 and one order-32 matrix has
norm 80. Independent integer enumeration exhausted all 33,554,432 and
2,147,483,648 projective states respectively. The latter matrix has 2,240
states at each extremum +/-80. The norm 80 was independently recovered
through an order-16 source with norm 32 whose signed diagonal completion
is Hadamard. No new equivalence class or external record claim is made.

The searches varied whole signings, cross blocks, source classes, diagonals,
structured/principal seeds, spectral proposals, edge/row/clique/rectangle
mutations, and symmetric Hadamard trades. Same source norm did not imply
the same lifting behavior: different norm-32 order-16 sources gave optimized
diagonal lifts with norms 96, 88, and 80. Minimizing source peak multiplicity
did not automatically improve the lift. The final corrected order-19 source
with norm 41 gave a coherent order-38 diagonal-family minimum of 121,
improving the older norm-43 source's minimum 129 within that family.
This is not an order-38 record: CORE's conference construction already
gives a better general upper bound there.

For the coherent lifts, set `K0=[[A,A],[A,-A]]`. Its zero matching entries
are not admissible signs. Stored profiles evaluate its norm exactly, and
changing those n zeros to signs changes the norm by at most n. The
small-source audit is only a finite observation. The stronger bound
`Phi([[A,A+D],[A+D,-A]]) <= 2sqrt(2) Phi(A)+n` for every complete
source with `Phi(A)=O(n^(3/2))` MUST NOT be inferred: the existing
[coherent clique-flip counterfamily](../NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md)
already disproves any o(n^(3/2)) repair on the entire class
`Phi(A)=O(n^(3/2))`, via `Phi(K0)=2 L_cl(A)`. That counterfamily has a large
leading constant and does not settle the near-minimizer implication.

## Counts, corrections, and evidence limits

The campaign records at least **28,555,014 exhaustive candidate-score
evaluations**: 9,633,741 CPU, 18,509,509 V100, and 411,764 Nuka.
These are calls, not globally unique matrices or permutation/switching
classes. Per-wave uniqueness has only its explicitly recorded scope.
Validation replays are not new candidate variants. Nuka also performed
102,592 additional heuristic-only candidate records, separately from these
totals. Some exactly scored candidates also received heuristic pre-screens;
these record counts are not the number of all heuristic scoring passes.

- Larger-order local spin searches give LOWER bounds on the norm. All 16
  shortlisted order-64 trade candidates, including heuristic scores as low
  as 244, had exact norm 256 once Boolean eigensign witnesses were found.
- Canonical Sylvester matrices have checked witnesses of energy 676 at
  order 128 and 5,418 at order 512. These are lower bounds, not upper bounds
  or exhaustive maxima. They cannot give candidate upper bounds on m_n.
  The traded order-128 bank is likewise heuristic.
- The new Nuka kick kernel initially compared an unsigned absolute value
  with a negative sentinel. Its first witness check rejected the run.
  Accepted repaired and rejected versions are separately retained; the
  earlier exact and trade kernels were unaffected.
- Staging adaptive CPU v2 accidentally overwrote parts of the OLD v1 logs.
  Authentic snapshots and mixed raw logs are preserved; only matching,
  valid per-host records enter reconstructed v1 counts. Those counts are
  lower bounds. Initial, sidecar, and v2 logs were unaffected; all three
  v2 processes reached their natural 1,800-second caps.
- V100 source wave7 used an older seed bank at orders 19--22. Both input
  versions are retained; the corrected norm-41 order-19 pass is separate.
  Wave1's 98,688 source-search evaluations are reproducible from code and
  seed, but only its six selected sources were archived. Later source
  histories, cross candidates, and all diagonal norm vectors were retained.

Full counts, archive checksums, and extraction roots are pinned in
[manifest.json](manifest.json). All owned search jobs finished naturally.
The research commit, verified large-drive backup and exact backup coverage
are recorded in [publication.json](publication.json).

## Preserved packages and replays

The four archives retain campaign data rather than only the winners:
scripts, input snapshots, retained raw logs/candidates, exact profile
arrays, witnesses, scoped counts, integrity checks, and replay instructions.

- `cpu_campaign.tar.gz`: structured/block search, adaptive v1/v2, principal
  sidecars, winner bank, independent integer verification including the
  odd-order GPU winners, and order-64 eigensign certificates. See its
  `PROVENANCE.md`.
- `v100_campaign.tar.gz`: source/cross searches and amortized all-diagonal
  sweeps. See its `README.md` and `manifest.json`.
- `nuka_campaign.tar.gz`: exact/heuristic evolutionary searches, Hadamard
  trade scores, odd principal subsets, and strong Walsh witnesses. See its
  `INTEGRITY_REPORT.md` and self-contained replay launcher.
- `independent_and_hadamard.tar.gz`: integer field-update verifiers and
  original receipts, all three 2,048-matrix Hadamard construction banks,
  and explicitly unpublished analytic scratch notes. Those scratch notes
  are not reviewed theorem artifacts or a global proof registry entry.

For a compact independent norm replay, extract the last archive into a
fresh directory, enter its extraction root, then run on a compute host:

```sh
g++ -O3 -march=native -std=c++17 original_mo_broad_independent_verify32.cpp -o verify
./verify independent_winners/v100_n26_61.txt
./verify independent_winners/rx_n32_0.txt
```

The research snapshot's focused registry/documentation check passed 38 tests.
From a compatible repository environment, its portable replay is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q -n 14 tests/test_original_mo_status.py tests/test_main_chain_docs.py::test_main_chain_L_open_and_docs_ok
```

Replaying an adaptive trajectory is scheduling/time dependent. Replaying
an actual retained matrix's exact norm is not. Do not relaunch unchanged
broad searches merely for cleaner receipts. The unresolved target remains
an all-orders convergence or nonconvergence argument; none of these finite
experiments supplies it.
