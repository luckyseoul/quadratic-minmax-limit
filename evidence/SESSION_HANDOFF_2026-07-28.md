# Session handoff — 2026-07-28 (quota / shutdown)

**For:** next agent / tomorrow  
**Repo:** `/home/nick/quadratic-minmax-limit`  
**HEAD:** `c0284dd` (Props 15.45–15.47)  
**L = lim α_n:** still **OPEN** — do not soft-close Main Theorem.

## Goal

Settle existence of \(L=\lim\alpha_n\) with load-bearing proof (AI bounty path).  
Acceptance: settled + sandwich/ρ=1 intact + pytest green.  
**Forbidden:** soft-close from sandwich+denseness alone; \(m_n\ge\Phi-2\) from 15.40 alone (F13).

## Where you are

**Settlement chain still open:**
\[
m_n\ge\Phi-2 \text{ on }n=p^2+1
\;\Rightarrow\; \mathrm{E}(1)
\;\Rightarrow\; L=\tfrac12
\]
via denseness (Prop 6.2). Equivalence: \(m_n\ge\Phi-2\) iff no-descent on all gap-2 undercutters.

**Locked this arc (use, don’t re-prove):**
| Item | Ref |
|------|-----|
| Props 15.40–15.47 | `solution.md` |
| Bi-tight empty when \(g_{\min}>-(p-2)/(p(2p-1))\) | Prop 15.47 |
| Cert \(g_{\min},h_{\min}\) at \(p=5,7\) | `evidence/e1_bitight_gsum_obstruction.json`, `e1_gmin_p7.json` |
| Deep covers \(p=5\): small \(k\) infeas; large \(k\) spike | `e1_deep_sweep_p5.json`, `e1_deep_cover_phi.json` |
| Failure graph F1–F14 | `evidence/E1_FAILURE_GRAPH.md` |
| Pytest | 68 green last full run |

**Type I + deep-tight:** closed at \(p=5,7\) via bi-tight obstruction.  
**Still open for all \(p\ge5\):** uniform \(g_{\min}\) bound; deep non-tight gap-2.

## Unfinished job (interrupted)

**Task:** closed-form lower bound on disjoint-pair
\(g_{\min}=\min\mathbb E_+[f_ef_{e'}]\) so Prop 15.47 applies for every \(p\ge5\),
**or** prove deep non-tight always \(\Phi\ge\Phi-2\) / ND.

**Facts already computed:**
- \(p=5\): \(g_{\min}=-3/65\), \(h_{\min}=-6/65=2g_{\min}\)
- \(p=7\): \(g_{\min}=-436/11452\approx-0.03807\), \(h_{\min}=2g_{\min}\)
- Threshold bi-tight: \(-(p-2)/(p(2p-1))\) (equals \(-1/15\) at \(p=5\))
- Wedge \(G=\pm1/p\) exactly; sum of wedge \(G=0\) (from \(C^3\) diagonal 0)
- Avg disjoint \(G=E(n/2-1)/(2\,n_{\mathrm{disj}})>0\) (2-design only)

**Dead ends (do not re-run):**
1. 4-point boolean LP with pairwise \(E[y_ay_b]=C_{ab}/p\) → bound \(\sim-0.6\) (too weak)
2. Same + Paley patterns + \(C^2\) row bounds + projector → still \(\sim-0.6\)
3. Chebyshev / 1-outlier from disj mean+var → \(\ll-1\)
4. Claiming Wick \(-1/p^2\) as lower bound (false: actual more negative)
5. Soft-close L; serial MILP enum (F1); ignore failure graph (F14)

**Promising next edges:**
1. **Association scheme on edges** of conference/2-graph: \(G_{ee'}\) constant on refined relations (Johnson + \(C\)-signs); read off min disjoint class from eigenmatrix. Gram spectrum of \(G\) partially known at \(p=5\) (ev \(n/2,88/13,72/13,40/13,0\)).
2. **Character sum** via halfspace + \(\mathrm{P}\Gamma\mathrm{L}\) orbit for \(E[y_ay_by_cy_d]\) (construction in `minmax_quadratic.halfspace_boolean_vector` / `paley_conference_prime_power`).
3. **Deep path:** prove every deep two-sided \(F\) has \(\Phi(C\oplus F)\ge\Phi-2\) via Prop 15.46 (Max− 1-bit when \(s_-=-2\) and some \(\tau\le-(p+1)/2\)); force existence of such \(\tau\) from degree/\(S=-2\) level set. Large-\(k\) \(p=5\) covers already spike above \(\Phi(C)\).

## Concrete next steps (priority order)

1. Derive closed-form \(g_{\min}(p)\) (or bound \(>\) bi-tight threshold) — scheme or character sum.  
2. If (1) holds for all \(p\ge5\): Type I + deep-tight ND unconditional for all such \(p\).  
3. Kill deep non-tight gap-2 (spike lemma forall, or ND).  
4. **Only then:** set Main Theorem / HANDOFF to settled \(L=\tfrac12\), add tests+evidence, full pytest.  
5. Consult `E1_FAILURE_GRAPH.md` before any multi-minute compute; use ProcessPool from **files**, not stdin; \(W\approx nproc-2\).

## Commands

```bash
cd /home/nick/quadratic-minmax-limit
# p=5 Max± caches often under /tmp/maxplus_p5.npy, /tmp/e1_deep_p5/
# p=7 Max± under /tmp/e1_p7/maxplus.npy, maxminus.npy (rebuild: src/e1_gmin_p7.py, e1_bitight_p7.py)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -m pytest tests/test_minmax.py -q
```

## Suggested skills

- `use-available-compute` before heavy jobs  
- failure graph discipline (project `E1_FAILURE_GRAPH.md`)  
- no soft-close / verification-before-completion when claiming settle

## Do not

- Soft-close Main Theorem without the chain  
- Re-run p=5 bi-tight MILP (already known + Gsum proof)  
- Exhaust matchings / serial SA (F1/F2/F8)
