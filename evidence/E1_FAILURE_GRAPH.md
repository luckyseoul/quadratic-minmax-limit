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

**Do not spend turns on F1–F14.** Next work must be a **proof edge** on path A or B (or a true parallel cert with W≥40 that writes durable JSON once), not another serial census.

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

1. Run `compute-budget.sh` before any >10s job; record `W`.
2. Independent units = seeds, matchings, residue classes, p-values — fan out with `ProcessPoolExecutor(W)`.
3. Serial only if algorithm is inherently serial; **state that in the log**.
4. Verification: one parallel wave, one JSON, one pytest — not live-serial-per-cover.

---

## Next allowed moves (only)

1. **Prove residual Prop 15.42 no-descent** (Type I large-k / deep non-tight k>2p) — closes m_n≥Φ−2. Done: dichotomy, tight S≡1/2, Type I large N_1.
2. **Proof draft** of Max-cover ⇒ clique-flip (or other spike) using design/AG(2,p) — no new SA.
3. **Proof draft** of k_⋆ / path-cycle for E(1) — no soft close.
4. **If** numeric cert: rewrite as `W`-worker batch, verify offline, one evidence file; mark incomplete enum as incomplete.

**Do not:** restart sequential MILP enum; re-run 20k random PM searches; re-derive sandwich; claim m_n≥Φ−2 from 15.40 alone (F13).
