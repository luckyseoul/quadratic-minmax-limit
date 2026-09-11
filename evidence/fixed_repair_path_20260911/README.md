# Integral correction and fixed repair paths

2026-09-11. Finite diagnostics on the archived order-11 source and aligned
order-21 target. The original convergence problem remains OPEN.

## Actual integration of the added techniques

`navier-stokes-techniques/quadratic_repair.py` imports the added
`CorrectionCycle` and `general_cone_representation` implementations.
For a complete signing a and every projective Boolean spin x, it computes
Q_a(x) exactly. The residual is max(0, Phi(a)-18) on order-12 inputs.
An edge flip e changes energy by delta_e(x)=-2 a_e x_i x_j.
For violating states, the cone solve uses directions
-sign(Q_a(x))*delta_e(x), together with negative coordinate slack directions,
to seek a nonnegative combination dominating |Q_a(x)|-18.
Its fractional coefficients rank twelve edges. Each single flip and each
pair among those twelve is scored on ALL 2,048 projective spin states.
Only a strict reduction of the exact norm is accepted, and the cycle
recomputes the complete residual. No fractional signing is accepted.

The inputs are the ten distinct one-vertex additions from the same aligned
target, restoring the original source's two mismatched edges. All 66 edges
may be repaired. This local experiment therefore does not impose the
two-repair budget of the separate fixed-path experiment below.

Soulkiller evaluated twenty independent source/method tasks with twenty
processes and one numerical-library thread per process. The baseline
uses single flips. Two additional controls use the same 66 pair proposals
per stage as the cone method: first twelve edge indices, and a fixed
seed-731 random choice of twelve edges. They use the same stopping rules.
This matches per-stage candidate budgets; trajectories and number of
accepted stages can differ.

| Added aligned vertex | Single flip | Cone pairs | Index pairs | Random pairs |
| --- | ---: | ---: | ---: | ---: |
| 11 | 22 | 22 | 22 | 22 |
| 12 | 22 | 20 | 22 | 22 |
| 13 | 22 | 22 | 22 | 22 |
| 14 | 18 | 18 | 18 | 18 |
| 15 | 22 | 22 | 22 | 22 |
| 16 | 22 | 22 | 22 | 22 |
| 17 | 22 | 22 | 22 | 22 |
| 18 | 20 | 20 | 20 | 20 |
| 19 | 22 | 20 | 22 | 22 |
| 20 | 22 | 22 | 22 | 22 |

The cone method escapes two single-edge stalls, also beating these two
specific controls. It does not increase the number reaching the target 18.
These ten prescribed starts and one random seed do not establish a general
advantage, speedup, optimality, or an all-orders contraction estimate.
The cone/baseline batch took 0.7013 seconds; the two control batches together
took 0.1930 seconds. These aggregate times are not a matched speedup ratio.
Raw stages and final signings are in `cone_result.json` and
`control_result.json`. NUKA independently replayed all 72 recorded stages
by direct integer pair sums over 2,048 spins; all passed.

## Separate exhaustive scheduling obstruction

`result.json` scores all 4,096 (added vertices, repaired old edges) states
for the archived fixed witness. Ten addition actions and two one-time
repair actions give 12! possible interleavings. Every added edge is fixed
by the aligned target; repairs are exactly old edges (3,8) and (5,10).
Each state is a genuine complete signing. Integer Gray-code enumeration
checks every projective spin, its zero mean, second moment binomial(n,2),
and a directly evaluated norm witness.

Dynamic programming uses cost Phi^2/n^3 and exact rational comparisons:
D(s)=min_predecessor max(D(predecessor),cost(s)). This induction covers
every permitted action ordering. Source cost is included.

| Permitted ordering | Exact optimum squared peak | Peak |
| --- | --- | ---: |
| All interleavings | 289/1024 | 0.53125 |
| Both repairs first | 961/3375 | 0.53361104 |
| Both repairs last | 1681/5832 | 0.53687737 |

The starting normalized norm is 17/11^(3/2)=0.46597208; the target is
44/21^(3/2)=0.45721844. Thus this witness cannot connect them while staying
below the starting level. This excludes that proposed finite bridge, not
other embeddings, extra repairs, other targets, or convergence.
Scoring took 3.0174 seconds on Soulkiller with 88 workers. NUKA independently
scored all 33 distinct recovered path states using binary enumeration and
direct pair sums, and checked three strict-threshold reachability
obstructions on the full stored table. No path variation optimum is claimed.

## Reproduction

The archived input member is
`original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json`
inside `../original_mo_growth_checkpoint_20260906/checkpoint.tar.gz`.
Its archive SHA-256 is
`f18159b14045fe5d64542cad749ada4a1ea76e66f578d5f88dfb4e04970d5ecd`.
Both output schemas retain the member hash. From the repository root:

```sh
g++ -O3 -fopenmp scripts/fixed_repair_scores.cpp -o /tmp/fixed-repair-scorer
python3 scripts/fixed_repair_path.py --archive evidence/original_mo_growth_checkpoint_20260906/checkpoint.tar.gz --scorer /tmp/fixed-repair-scorer --workers 88 --output /tmp/fixed-repair-result.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 navier-stokes-techniques/quadratic_repair.py --archive evidence/original_mo_growth_checkpoint_20260906/checkpoint.tar.gz --workers 20 --output /tmp/cone-repair-result.json
```

Choose worker counts within the execution host's available budget. Run
the second command with `--modes index_pair random_pair` and a fresh output
for controls. The two `scripts/verify_*repair*.py` files independently check
the saved outputs. These are finite tools, without autonomous goal selection.

The next unresolved implication is a construction with repair cost and
normalized defect controlled uniformly over increasing orders. This
experiment supplies neither the requisite family nor its uniform bound.
