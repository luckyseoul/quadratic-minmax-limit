# NS-technique port at the order-26 record — exact negative results (do not repeat as an L-route)

**Date:** 2026-09-12 (night session)
**Status of this artifact:** finite exact results; **no status change**; recorded so the
imported Navier–Stokes techniques are not re-run as a route toward E(1) or `L`.

## Why this exists

On 2026-09-11/12 the "Navier–Stokes" port (OpenAI blowup-paper techniques) was tried in
scoped scouts: covariance-cone majorizers (`evidence/navier_covariance_cone_probe_20260911/`),
relative-gauge connected-layer closure (`evidence/navier_relative_gauge_layers_20260911/`),
cone-ranked repair on cross blocks (`evidence/paley_two_block_cone_repair_20260911/`,
`evidence/paley_structural_lift_20260911/`, `evidence/fixed_repair_path_20260911/`). All were
negative or obstruction results. On 2026-09-12 the techniques were executed **at full strength,
exactly**, on the strongest live finite target: the verified order-26 record signing with
norm 61. The outcome is recorded here so no future session repeats this as a hoped-for path to
E(1): the transfer is bookkeeping (stress set / correction cycle / exact scoring), not an engine,
and its ceiling is finite local statements at fixed n.

## Target signing

- Source: `scratch/campaign/.../independent_winners/v100_n26_61.txt` (broad campaign
  2026-09-06; noted in `evidence/original_mo_broad_campaign_20260906/README.md`), also as JSON in
  `mo-nuka-source-20260906-C3PMKY/inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json`.
- Exact full-cube norm confirmed in this session: max |sum_{i<j} a_ij x_i x_j| = 61 over 2^25
  projective states. All Q(x) are odd (325 terms), so "undercut" means reaching <= 59.

## Method (the port, at strength)

1. **Stress set** = the complete tied-maximizer band: all 676 states with |Q| = 61 plus all
   2236 states with |Q| = 59 — nothing sampled (the q=13 lesson about incomplete stress sets).
2. **Exact scoring only**: every candidate is scored on the complete 2^25-state projective
   cube via half/quarter-cube statistics; no fractional or sampled acceptances.
3. **Exact k-flip reduction**: each flip changes Q(x) by exactly +-2, so with k flips only
   states with |Q| >= 60 - 2k can exceed 59. This makes "does any k-flip repair undercut 61?"
   an exact finite 0/1 problem: choose k edges so that at every live state
   sigma(x) * sum_{e in S} a_e x_i x_j >= tau(x) (tau = 1 iff |Q(x)| = 61). FEASIBLE gives a
   genuine undercut; INFEASIBLE proves no k-flip repair exists. Verified by CP-SAT.

## Exact results (complete scans / exact solves)

| k flips | result |
| --- | --- |
| 1 | all 325 single flips -> 63 (every one) |
| 2 | all 52,650 pairs -> 65 (every one) |
| 3 | all 5,668,650 triples -> minimum 67 (e.g. {(0,1),(0,2),(0,3)}) |
| 4 | CP-SAT INFEASIBLE (no repair <= 59) |
| 5 | CP-SAT INFEASIBLE (8 workers, 132 s) |
| 6 | CP-SAT INFEASIBLE (128 s) |
| 7 | CP-SAT INFEASIBLE (731 s) |
| 8 | CP-SAT INFEASIBLE (304 s; 300,144 live rows) |
| 9 | UNKNOWN at 805 s cap (615,264 live rows; 2 workers) |
| 10 | not resolved (same size class as k=9) |

- k = 1, 2 INFEASIBLE cross-validated by the independent brute-force scans (which exhaustively
  compute every flip set, not just feasibility at 59).
- Heuristic SA (200 restarts x 3000 swaps, k = 4..8): best live-max 65/67/65/67/69 —
  no undercut found.
- Band census (complete): |Q| in {35: 1005836, 37: 786500, 43: 315120, 45: 216684,
  51: 52260, 53: 28288, 59: 2236, 61: 676}; even values and 39/41/47/49/55/57 are absent.

## Structural reading

Every single flip damages by exactly +2; every pair by exactly +4; every triple is at least +6
(minimum attained). The maximizer family is a 3-covering system on the 325 edges: every 1-,
2-, and 3-subset of edges is simultaneously damaging at some tied maximizer. Hence no short
correction cycle can descend — the record is repair-rigid to depth 8 (exact).

## Limits (honest, per this program's standards)

- Local to this signing; no global statement about m_26 and no asymptotic claim.
- k >= 9 open (resource-capped). Heuristic search found nothing at k <= 8.
- This is a finite diagnostic; it is **not** a route to E(1), `L`, or the limit, and should
  not be re-run as one. It is also not evidence for or against convergence.

## Reproduction

Runs use the V100 venv: `~/.venvs/mo-exact/bin/python`. Work dir in scratch:
`/home/nick/scratch/ns_try_20260912/` (contains the larger `live43.npz` 615k-row band and
`pairS.npy`, regenerable from the scripts here).

```sh
python scripts/ns_undercut.py            # stage 0/1: norm check, 1-flip complete scan
python scripts/collect_live.py           # |Q| >= 43 band  -> live43.npz
python scripts/stage2b.py                # pairs (uses stage01.npz)
python scripts/triple_check.py           # all triples (exhaustive, ~24 s)
python scripts/sa_k.py                   # heuristic k = 4..8
python scripts/cp_sat3w2.py 5 600        # exact k-flip feasibility (workers=2; 88 workers hangs)
```

## File hashes (SHA-256)

```
collect_live.py   a990952b2c7572b1ea9afdb4b99ecd037320b3c64ef1a91d1373fa1a194bd710
cp_sat2.py        9d3263f1c986585cff94ba7f16ce30bc6a032805ecd7eb61728da469087b6e0d
cp_sat3.py        1a2663e0d2c662871e542abc2c1f091652cc3dc9891d6cf42c260edb60b51be1
cp_sat3w2.py      19a5612e118b9cbbc9cf1d182a21deb9c3c7db7dbcf3ceae8ae0aea28415633c
ns_undercut.py    8553d01ec24f9a4c9c903ff08b82af66006415cc8f0cc9f3db0caa6134347c66
sa_k.py           edaf1b08e2cbaecdb0232e7ffb85895f7f1befa8f11c5e7612be6d7725a28366
stage2b.py        6e01e80cfdd0e1bf99b235d5907280e857826a51f5b97bfb402c10b40fee5b02
stage2_pairs.py   6ff1526f5edb6be6a48e623972bbc065ae2d17e44f6bff90345dc31568918199
triple_check.py   b2f78fb0f1240737c6d4754056ee6cba0ef405f76c6d4b8edbfb3b77c89250a3
stage01.npz       886554641a82a4ecc59508a7b62d6d308f5976a9cc00032e571cb34f0cc9fe96
```

Verbatim run receipts: `receipts.txt`; structured numbers: `results.json`.
