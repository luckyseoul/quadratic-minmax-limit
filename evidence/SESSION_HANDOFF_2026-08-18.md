# Session handoff 2026-08-18

## Goal
Close leftovers so E(1) is Max+-free for all primes p≥5; then `e1_closed_general` / writeup L=1/2 only by import. Repo `/home/nick/quadratic-minmax-limit`. HEAD after this commit: 15.544. Scratch: `/tmp/grok-goal-a558c5f11751/implementer`.

## Status
Three leftovers open. Lemma D stays True. Aut-Schur / Gsum / pairing stay False. Soft-close forbidden. Do not unflip e1 as a ritual (live e1 is the old AND, True).

| Flag | Value | Note |
|------|-------|------|
| e1_closed_general | True | old AND; not acceptance |
| phi_F_ge_6 | False | leftover-1 |
| residual_ii k≥4p | False | leftover+splus empty all nF at p=5 k=20 (15.528) |
| type_I multilevel | False | Aut_e G>T open for p>5 |
| 15.272 / residual i | True | do not unflip |

## Shipped this arc
- 15.539 `n_free=c(p−1)⇔D=D_lattice(c)`; live pin only p=5,7
- 15.540 Type I `A_1d(r)=−4p³/(p−2)` on F_p, else 0
- 15.541 `c_eq` named Hoffman floor; fail ceil
- 15.542 half-net count = \|H+\| at p=5,7 only (census)
- 15.543 `μ_1d=κ/(p(p−2))` on \|κ\|=1
- 15.544 p=5 \|μ_full\|=1/p²; mix G>T; p=7 not a p-law

## Leftover-1 (binding)
Prove **NUM_SUM=16p A** for all odd p≥5 (A=2 n_1d−4 n_pp−3k), or a Max+-free name of Q_NL / Q_τ on 15.290 types. Fail-when-wrong must break the formula. Do **not** import phi_F until ⟨δ,ψ⟩≤2.

Dead / do not rerun: 57+57, n_bar0, occupancy, Aut involutions, Z[i] 48, raw ReJ/J/|J|² match of {32,48,1376,1496,1544} (0 hits), N* interpolant 8(A+n_pp−2p+2) (S0 D≠D_eq at p=13), 15.507 as a p≡3 close, more half-net census.

Live numbers: Q++=48/13, 1544/409; NUM_SUM=480, 21616=16p A at p=5,7. Q++=w Q_1d+(1−w)Q_NL; Q_1d named (15.356); Q_NL=64/15, 632/171 unnamed. 15.506 p=7 A-shifts die at p=11. 15.502: 8A/D + 2(p−4) N* shift is not a p≡1 theorem.

## Residual (ii)
p=5 k=20 leftover+splus S≥2 empty all nF (15.528). k=22 leftover+splus 3600s TLE (any-nF and nF 7–14). SK HiGHS k=24 7200s/24t still running (`sk_resii_k24.py`); no JSON yet. Official is s+≥2 not s+=2.

## Type I
Δ_conn / I / A_1d named; still den D. Aut_e G>T open for p>5. Next: name μ_full / A_full in p, or Aut_e on the mix.

## Remotes
JF 192.168.1.191 CPU when SK occupied (no NVIDIA). Orin 100.67.236.54 W=4. Do not duplicate SK MIPs. Do not put passwords in the repo.

## Do not
Commit 15.495 dirt / 15.496 / 15.530. `/goal clear`. Prize / Paata / ping in public files. Flip leftover flags without a fail-when-wrong unit.

## Next
1. Prove NUM_SUM=16p A (character-sum / 1D+NL mass) or name Q_NL.
2. Harvest k=24 when it exits.
3. Claude-referee primary on leftover-1, then openai-referee same slot. Branch only on BLOCK.

## Suggested skills
- agent-cost-optimization
- graph-engineered-completion
- use-available-compute
- claude-referee (primary; use often on stuck leftover-1 / after major units)
- openai-referee (secondary; same tool and claim after Claude)
- verification-before-completion
- handoff
- scientific-critique
- grill-me
- research
- arxiv
