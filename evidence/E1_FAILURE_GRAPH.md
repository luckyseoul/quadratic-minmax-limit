# E(1) attack failure graph (do not repeat)

**Purpose:** Prevent re-running known dead loops. Update this file when a path dies.

**Machine:** 88 cores. Default: `full_workers ≈ nproc-2`. Never leave multi-minute work on 1 core without an explicit "inherently serial" note.

---

## Failure modes (BANNED unless new evidence reopens)

| ID | Mode | Why it failed / wastes time | Symptom |
|----|------|----------------------------|---------|
| **F1** | Sequential MILP no-good enum | One HiGHS solve → serial `maxR`/`phi_mitm`/clique-flip per cover | ~1–2 cores busy for hours; 88 idle |
| **F2** | Sequential SA seed loop for covers | Cover rate tiny; seeds independent → **must** ProcessPool | Single `for seed` |
| **F3** | Soft-close sandwich+denseness+ρ=1 | Explicitly forbidden; does not settle liminf=limsup | HANDOFF still OPEN correctly |
| **F4** | Treat continuous Γ bound as discrete lower bound | Cont ≥ 61 ∀z but discrete maxR on level can be 44; only **max**_z discrete matters | False “proof” |
| **F5** | Claim op=√41 for all Max-covers | MILP found covers with op≈6.50, still Φ=Φ(C), flip OK | Prop 15.36 overclaim |
| **F6** | GW/SDP lift cont→discrete | SDP·α < 60 on best Γ | Theory path dead at p=5 |
| **F7** | “Every S=1 y clique-flips” | On covers only ~24/168; many have E[Σ²]=1, Σ∈{±1} | Wrong quantifier |
| **F8** | Exhaust all PMs | 7.9e12 matchings | Impossible |
| **F9** | Re-derive BH / multipartite / pure SA as deliverable | MO re-audit: only load-bearing E(1) or permanent gap counts | Scope thrash |
| **F10** | Single-thread maxR without numba/process fan-out | 2^13×2^13 per matching on one core | Slow verify |
| **F11** | Poll/sleep on background jobs | Use completion notifications | Wasted turns |
| **F12** | Kill commands with `pgrep -f` matching the wrapper | Self-kill / no-op | Ops failure |
| **F13** | Claim \(m_n\ge\Phi-2\) from edge-minimal gap alone | Abstract 2-Lipschitz local min on hypercube can drop by \(2d\) at distance \(d\); Prop 15.40 only bounds edge-minimal undercutters | False E(1) |
| **F14** | Ignore this graph and re-run dead loops | Session thrash: serial MILP/SA, false m_n shortcut, soft-close, single-core re-census | User has to intervene; no new proof edge |
| **F15** | Plain Fréchet on Max+ conditional cov for \(g_{\min}\) | Cond means match Gaussian (Prop 15.50) but Fréchet only gives \(g_{\min}\gtrsim-0.4\) at \(p=5\), below bi-tight thresh | False hope of \(L(p)\) from 2-point Frechet |
| **F16** | Pin free modulus \(c\) by max \(g_{\min}\) under PSD of \(G(c)\) | At \(p=5\), PSD+rank hold on a continuum; max PSD \(g_{\min}\approx-0.040\) at wrong \(c\); true \(g_{\min}=-3/65\) needs \(\mathrm{Tr}(G^2)\)/spectrum | Spurious “better” gmin; wrong pin |
| **F17** | **Any** multi-minute CPU job on **1 core** (pytest without `-n W`, or research script with ProcessPool theater then serial main, or pure-Python `for quad in binom(n,4)` while 87 cores idle) | 88 cores idle; user has rebuked this **repeatedly** (including this session: m4_pseudo ~97% NLWP=1 after “F17 fixed”); destroys trust; wastes wall time | `ps`: one python `pcpu≈100` `nlwp=1`; load≪nproc; parent does heavy work after `as_completed` |
| **F18** | Character sums / moments on affine or PGL orbit of halfspace as full Max+ | Orbit size 60 of 260 at \(p=5\) (PGL+Frob+sign); affine gmin ≈ −0.6 ≠ −3/65 | Incomplete orbit; wrong \(g_{\min}\) |
| **F19** | **Moduli class-invariant thrash** after constancy already achieved | type6+CR already m4-constant at p=5 (26) and p=7 (48); coarse+CR / +κ / e4 / denser evec **do not** drop p=7 nullity below 2; local multi-param grid exceeds M_cand | New scripts `refine*`, `pin_extra`, more keys; no proof edge; user rebukes loops |
| **F20** | **GPU theater** — claim GPU while wall is CPU | CuPy m4 for 0.1–0.4s then hours of ProcessPool class-build/evec with `nvidia-smi` util 0%; evidence still says `use_gpu=True` | User: “haven’t seen you use the GPU a single time” |
| **F21** | Residual-(i) theater: ship a 15.xxx whose `gsum`/`type_I`/`e1` stay False | Does not settle L; LONG_HORIZON_GOAL forbids it | New `prop152xx.py` with `proved=True` on an identity that does not flip leftover |
| **F22** | Flip `e1_closed_general` on census \(p\le7\) | Envelope / max\(\lvert\mu\rvert\) hold at \(p=5,7\) only | Predicate True without Max+-free general hinge |
| **F23** | Infer \(k_\star=o(n^2)\) from edge-local optimality, switching-minimal cut inequalities, Max-Lipschitz, or second moments of \(A\circ C\) | For every fixed \(C\), random signings followed by greedy edge descent produce edge-local minima with both cube norms \(O(n^{3/2})\) but signed switching distance \(N/2-O(n^{3/2})\). The product/frame second moments are exactly independent of \(A\). | A locally stable, correct-scale signing is still \(\Theta(n^2)\) from Paley; global/cardinality minimality is indispensable. See `NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md`. |
| **F24** | Infer proximity to the Paley class from spectral defect \(\delta=o(1)\), or from \(\delta=0\) | \(\delta=0\) identifies the union of all conference classes. Inequivalent classes have identical operator spectrum, and Mathon's order-\(5r^2+1\) classes have \(\rho<1\) at every order without a known uniform gap. | A spectral-rigidity proof must also control the Boolean radius of the resulting conference class; Proposition 15.14 cannot silently replace “some conference matrix” by the chosen Paley matrix. |
| **F25** | Apply ordinary Delsarte/external-distance, a fixed-\(L^q\) moment norm, or any uniform fixed-temperature free-energy lower bound to the augmented cut code | Fixed moments lose the required \(n^{3/2}\) term. More decisively, Proposition 6.9 gives an infinite symmetric-conference counterfamily to the required free-energy bound for every fixed \(c>0\). | Only a growing-temperature criterion with error uniform in \(c_n\), or a different global mechanism, can reopen this route. See `NOTE_2026-09-01_FIXED_C_FREE_ENERGY_CONFERENCE_BARRIER.md`. |

### Residual-(i) dead mechanisms (2026-08-13; do not re-run)

Detail: `evidence/SESSION_HANDOFF_2026-08-13_residual_i.md`.

- Fréchet alone; joint Fréchet (two laws, opposite pairwise, same \(\mu\)) — only \(\lvert\mu\rvert\le1-2/p\).
- \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise — false at \(p=7\).
- \(\lvert\delta\rvert\le\mathrm{room}_\delta\) pointwise — false at \(p=5\).
- W0 (scheme dual) \(\in\mathrm{im}(\mathrm{Gsum})\) — resid 0.54 at \(p=5\).
- Type I + switching \(G=\pi(G)\) — rank \(132<600\) at \(p=7\).
- Low-degree Farkas \(y\) in \(f_e\) — 66 values, corr \(\sim-0.2\).
- 3-wise interpolant — \(\varphi\) not a function of the 6 \(K_4\) edges.
- \(\kappa_B\perp E_{\pm4p}\) — rel resid 0.997.
- Matching-PSD / one-edge \(\sum\mu\) — no pointwise.
- max-\(\lvert\mu\rvert\) at extreme \(\lvert\varphi\rvert\) — false at \(p=5\).
- Halfspace \(G\)-orbit as full \(\mu\) — differs \(\sim0.6\).
- \(\chi(\)never-zero quadratic\()\) in Max+ — \(0/300\) at \(p=5\).
- Unsigned Per / \(\sum\lvert\mathrm{per}\rvert\) — 15.231.
- Aut-line dim\(\le1\); 15.237 pair-span as a Type I kill.
- 4×4 Gram \(\{1,f_e,f_f,\chi_S\}\) — only \(\lvert m_4\rvert\le1-2/p\).
- Affine-quadratic level sets on AG(2,5) — only 60 linear halfspaces; 200 Max+ are not \(Q^{-1}(S)\).
- CR-class master+diamond LP — \(\max\lvert\mu\rvert=1\) at \(p=5,7\) (\(E_{\pm4p}\) survives).
- Single IP-valency \(K_4\) — regular at \(p=5\), not at \(p=7\).
- Global mix \(\mu=(1-t)\mu_{\mathrm{part}}+t f_4\) with \(t=t(p)\) — \(t\) depends on \((\kappa,\varphi)\).

**Live graph (update every turn):** `evidence/P0_ENGINEERING_GRAPH.md` — nodes, critical path, compute routing, change log. **Banner on that file was stale (claimed L CLOSED 2026-08-05); current L is OPEN — trust `HANDOFF.md` / `STATUS.md`.**

---

## Dependency graph (what actually closes L)

```
L = lim α_n
 └─ dense Paley / Prop 6.1–6.2
     └─ E(1) on n=p²+1 (ρ=1 family)  ← BLOCKER
         ├─ path A: k_⋆ = o(n^{3/2}) (or Max-Lip + k_⋆=o(n²))
         │    └─ needs structure of Φ-minimizers (path-cycle / Δ bound)  [OPEN]
         ├─ path B: matching non-undercut (p≥5) + matching dichotomy + non-matching control
         │    ├─ Prop 15.33 non-covers cannot undercut  [PROVED]
         │    ├─ Max-cover PM Φ≥Φ(C) ∀M, p≥5         [OPEN — F1/F4/F6/F7 failed as proofs]
         │    └─ dichotomy / k_⋆                         [OPEN]
         └─ path C: permanent relative gap → limsup < 1/2 along family  [no construction]
```

**Do not spend turns on F1–F17.** Next work must be a **proof edge** on path C residual (\(g_{\min}\) moduli + deep ND) or path A/B with new invariants — not serial census, F15/F16 reopens, or single-core full pytest.

### Full-suite compute (mandatory)

```bash
# ALWAYS (repo default pytest.ini uses -n 86; script recomputes W):
./scripts/pytest_full.sh
# or:
W=$(($(nproc)-2))
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py -n "$W" -q
```

`OMP_NUM_THREADS=1` is **per worker** (correct with xdist). The ban is **no real fan-out** on multi-minute jobs.

### F17 recurrence (2026-07-30) — treat as hard bug

**Filed:** `evidence/AGENT_BUG_F17_RECURRENCE.md` + session `feedback.jsonl`.

Anti-patterns that still count as F17 even if `require_workers()` was called:

1. **Pool theater:** `ProcessPool` for easy primes, then **serial** `p=7` (or any heavy leftover) on the parent.
2. **Unsharded loops:** one process walks `combinations(range(n),4)` / full Max+ census while W was only used for tiny tasks.
3. **Nested serial rebuild:** workers re-load and re-scan the full matrix alone with no shard.

**Mandatory self-check** after launching any >10s job:

```bash
ps -eo pid,pcpu,nlwp,cmd --sort=-pcpu | head -15
# FAIL if top job is python, pcpu>80, nlwp=1, and nproc>=8
```

On FAIL: kill, rewrite to shard, re-run. Do **not** wait for the user to notice.

Use `src/workers.py`: `require_workers()`, `pool()`, `assert_not_single_core_thrash()` (optional live check).

---

## Proved load-bearing nodes (keep)

- Sandwich, denseness, ρ=1 on n=p²+1
- Prop 15.33: non-cover ⇒ no undercut
- Props 15.30–15.32: spike / Γ / mod-4 tools (sufficient, not forall-M)
- Prop 15.31: clique-flip **sufficiency** (existence open)
- Prop 15.35: p=5 cover forces S=1, residue 1
- Prop 15.36 algebra: B=CD+DC blocks (proved); spectrum claims only certified samples
- 11+ SA/MILP covers: Φ=Φ(C), flip, maxR=60 (evidence, not forall)

---

## Open issues (graph)

| Node | Status |
|------|--------|
| Max-cover PM non-undercut ∀p≥5 | **open** |
| Continuous Γ ⇒ discrete maxR | **blocked** (F4, F6) |
| op=√41 for all covers | **false** as universal claim (F5) |
| Matching dichotomy / k_⋆ | **open** |
| Prop 15.38 n=10 two-sided k=5 | **proved** (only Δ=1 undercuts) |
| Prop 15.41 first-hit + dangerous-edge | **proved**; no-descent **open** (n=10 PM cert) |
| Prop 15.42 Max± dichotomy + tight no-descent | **partial**: dichotomy+tight+Type I large \(N_1\) proved; residual open |
| E(1) / lim α_n | **OPEN** |

---

## Compute policy (mandatory)

1. Run `compute-budget.sh` **and** check GPU (`nvidia-smi` / `src/gpu_budget.py`) before any >10s job; record `W` and whether GPU is free.
2. Independent irregular units = seeds, matchings, residue classes, p-values — fan out with `ProcessPoolExecutor(W)`, `OMP=1` per worker.
3. **Prefer GPU (CuPy / V100)** for large dense batches: m4 column products (`e1_gmin_m4_gpu.py`), Φ sampling, dense Rayleigh/power. One process owns CUDA — do not open 88 GPU contexts. Keep reductions on-device (argmax/atomic-style max), not full batch D2H.
4. **Atomic I/O when allowed** (`src/io_atomic.py`, same idea as Wieferich hunts):
   - **mmap** `np.load(..., mmap_mode='r')` for shared Max+ / large caches across ProcessPool workers (page cache shared; avoid 86× full copies).
   - **atomic writes**: write `*.tmp` → fsync → `os.replace` for evidence/checkpoints under concurrent jobs.
   - GPU hit/max collection via device reduce / `atomicAdd`-style patterns when scanning large spaces.
5. Serial only if algorithm is inherently serial; **state that in the log**, and still run other vectors in parallel or use GPU if it fits.
6. Verification: one parallel/GPU wave, one JSON, one pytest — not live-serial-per-cover.
7. Report in evidence: `workers`, `gpu.used`, `io` (mmap/atomic), wall time.

---

## Next allowed moves (only)

1. **Lift Prop 15.45:** prove \(g_{\min}>-1/p\) for all \(p\ge5\) (closes star force); prove bi-tight integral infeas for all \(p\ge5\); kill deep non-tight \(k>2p\).
2. **Proof draft** of Max-cover ⇒ clique-flip (or other spike) using design/AG(2,p) — no new SA.
3. **Proof draft** of k_⋆ / path-cycle for E(1) — no soft close.
4. **If** numeric cert: rewrite as `W`-worker batch, verify offline, one evidence file; mark incomplete enum as incomplete.

**Do not:** restart sequential MILP enum; re-run 20k random PM searches; re-derive sandwich; claim m_n≥Φ−2 from 15.40 alone (F13).
