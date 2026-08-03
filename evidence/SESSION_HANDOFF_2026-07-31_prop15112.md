# Session handoff — 2026-07-31 (Props 15.110–15.112)

## Goal / status
- **Goal:** close \(\lim\alpha_n\) (MO 413935) via Path C (or E(1)).
- **Status:** **L OPEN.** Residual structure advanced; **core inequality not proved** for general \(p\ge5\).

## What is actually open (one line)
Prove **\(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\)** (⇔ **\(\delta^2\le\rho_{\min}^2\)**) for all primes \(p\ge5\).  
For \(p\ge7\) this ⇒ orth≤room_hyp ⇒ 16N/bi-tight (with mult≥d−1). Then deep ND + Main Theorem.

## Real progress this arc (structure only)
| Prop | Proved | Not proved |
|------|--------|------------|
| 15.110 | ∑κ∏, e₄; ρ_min²<budget all p≥7; residual ⇐ δ²≤ρ_min² | δ²≤ρ_min² ∀p |
| 15.111 | α_κ, α_ρ, pair Schur; Φ−μ̄ driven only by δ | same |
| 15.112 | conf ‖κ‖²; design E[D²]=2n; ED4=ED4_flat+24δ² dictionary | ED4≤ED4_suf ∀p |

**Cert only:** residual / design moments at **p=3,5,7**.

## Key formulas (do not re-derive)
- \(\rho_{\min}^2=5n(p^2-1)(p^2+3)/(6p^2(p^2-5))\)
- \(\alpha_\kappa=(p^2+2)/(4p^2)\), \(\alpha_\rho=(7p^2+5)/(2p^2(p^2-5))\), pair\(=(p^2+11)/(4(p^2-5))\)
- On Z: \(E[(y^\top By)^2]=\bar\mu\|B\|^2+8\langle\delta,\kappa_B\rangle\); 16N ⇔ max⟨δ,κ_B⟩≤(n−10)/(n−6)
- ED4_suf = ED4_flat + 24 ρ_min²; ED4_suf < ED4_bud for p≥7

## Settled decisions
- Path C primary residual, not fresh p=5/7 census theater.
- **F19:** class_key is **not m₄-equitable at p=7** — do not close residual via moduli/class_key T.
- **F3:** never mark L closed from sandwich + denseness alone.
- Max+ is antipodal spherical 3-design skeleton when E[yyᵀ]=2P₊ (cert 3,5,7).

## Dead / weak (do not reopen)
- Unrestricted CS / crude LP / Gegenbauer on ED4 (too weak or false at p=7 for simple forms tried).
- E[Q₄]≤2Q₄(1)/N false at p=7.
- Aut-dim-1 Q₀ story is clean at **p=5 only**; not a general Aut-compressed proof.
- Jensen avg‖f_y−κ/p²‖² bound on ‖ρ‖² is uselessly loose.

## Next concrete steps
1. Bound **E[D⁴]** for boolean Cy=py antipodal 3-designs (weight enumerator / scheme-SDP / Gauss-sum), Max+-free if possible.
2. Or closed form for ED4 / δ² from Paley combinatorics.
3. If residual closes for all p≥5: ship bi-tight + update HANDOFF/graph; then N_DEEP; then Main Theorem only if full chain written.
4. Run `tests/test_prop15110.py`–`15112.py` after changes; keep `proved_*_for_all_p=false` until true.

## Paths
- Handoff: `HANDOFF.md` §0
- Graph: `evidence/P0_ENGINEERING_GRAPH.md`
- Solution notes: `solution.md` Props 15.110–15.112
- Code: `src/e1_gmin_m4_prop15110.py` … `prop15112.py`
- Evidence: `evidence/e1_gmin_m4_prop1511{0,1,2}.json`
- Max+: `/tmp/maxplus_p5.npy`, `/tmp/e1_p7/maxplus.npy`

## Suggested skills
`use-available-compute` · `agent-cost-optimization` · `goal-verifier` · `handoff` · optional `arxiv` for design/LP energy bounds

## Honest net
Better **attack surface** (one ED4 inequality + Schur channel). **No general residual proof.** lim α_n still open.
