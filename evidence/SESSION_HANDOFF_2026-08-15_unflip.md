# Session handoff — 2026-08-15 Aut-Schur unflip

**Status:** \(L\) OPEN. Residual (ii) CLOSED. Residual (i) OPEN. E(1) OPEN.

## What happened

A 15.270 Aut-Schur close was shipped: Jacquet “every irrep of \(\mathcal W_{++}^0\) has a \(U\)-invariant” was treated as “PSL-span of the \(k=3\) locked-triple space \(F\) equals \(\mathcal W_{++}^0\)”. Skeptic panel rejected it.

Witness: at \(p=5\), \(\dim\mathcal W_{++}^0=n(n-6)/8=65\), \(k=3\) Veronese rank \(61\) (15.212: \(k=1\cup k=3\) is \(65/65\)). Even-on-\(\Omega\) \(F\) has dim \(\sim p^2/4\); \(\dim\mathcal W_{++}^0\sim p^4/8\).

The harness then started verifying the leftover `solution.md` claim \(L=1/2\) (goal verification 1/6). That was the stale false close, not a new completion. Unflip restored OPEN and aborted that check.

## Live predicates (must stay False until a real hinge)

| flag | value | wiring |
|------|-------|--------|
| `psl_span_F_eq_Wpp0` | False | 15.270 failed-lift unit |
| `theorem_aut_schur` | False | same |
| `gplus_pd_proved_general` | False | AND includes Aut-Schur + span |
| `ker_sc_proved_general` | False | imports `gplus_pd` |
| `residual_i_closed_via_249` | False | `ker_sc` ∧ `cost_D` (cost_D True) |
| dual-eq empty | False | imports 15.249 |
| Type I / E(1) | False | import dual-eq |
| `gsum_disj_lb_proved_general` | False | unused |
| `solution.md` Main Theorem | OPEN | do not write \(L=1/2\) |

Singer/Fejer/\(\mu=0\)/Weil/Gershgorin/small DFT stay True as **structure** (PD of the \(k=3\) circulant on even-on-\(\Omega\)). They must not make `gplus_pd` True.

Honesty tests restored: `test_prop15170`, `15207`, `15216`, `15217`, `15218`, `15249`, `15268`, `15176`, `15270`. Pytest of that set passed.

## Still proved (do not reopen)

- Residual (ii) full ND (`residual_ii_full_closed`)
- Bi-tight at \(p=5\) (15.167)
- \(\nu=0\) on \(|\kappa|=1\) (15.268)
- 15.249 \(\mathrm{cost}_D<2-\alpha\) on scheme⊕cross (Weil)
- 15.207 reduction: \(\ker=\mathrm{sc}\) \(\Leftrightarrow\) \(G_+\succ0\) on \(\mathcal W_{++}^0\)
- \(p=5\) Veronese of \(k=1\cup k=3\) is \(65/65\) (15.212)

## Next hinge (one listed, Max+-free, all \(p\ge5\))

Do **not** flip predicates on structure.

1. **Preferred:** \(G_+\succ0\) on the **full** \(\mathcal W_{++}^0\) (all Max+ Veronese, Aut-SOS / \(P_\pm\) / \(G_+=B^*B\)). Then existing 15.207+15.249 wires residual (i).
2. **Alt:** \(N(\varphi)\neq0\) leftovers — D2 last parity when \(8\mid S\); D6NEQ \(p\equiv5\pmod{24}\), \(p\ge29\); GALB \(d\notin\{2,6\}\); even \(\hat P\neq0\). Scratch: `/tmp/grok-goal-838009eda84d/implementer/{D2_Q2S,D2_EEOO,D6_P5MOD24,GALB,GALB_D10}.md`.
3. **Dead:** envelope / reflection / \(\lvert\mu\rvert\le2/n\) (false at \(p=11\)); \(K_4\le\mathrm{Wick}_{hi}\) (false at \(p=11\)); \(\|m_4\|_2^2\le n(n-2)/4\) (false at \(p=13\) mix); Aut-Schur / Jacquet lift.

After a real bound: import only that flag, then two 15.170 dumps + honesty pytest, then `solution.md` \(L=1/2\). Soft-close forbidden.

## Files this session

- `src/e1_gmin_m4_prop15270.py` (new; Aut-Schur False)
- `src/e1_gmin_m4_prop15207.py`, `15216.py`, `15249.py`, `15170.py`, `e1_main_chain_status.py`
- `tests/test_prop15270.py` + honesty tests unflipped
- `solution.md` Main Theorem OPEN
- Do not treat scratch `ALT_HINGE.md` (“Proved for every prime”) as true

## Compute / process notes

- 88 cores, ~57 GiB, V100 idle. ProcessPool from module-level `.py` only; `OMP=1`. GPU unused for character sums.
- User: unlimited tokens this week; do **not** add `--budget`. Current goal is open-ended — do **not** `/goal clear` to “start fresh”.
- Do **not** re-search the MO/X thread after every compression.
- Referee MCP is optional (stuck / second opinion).
- Goal verification 1/6 was the harness checking the leftover false \(L=1/2\) page. Unflip is the correct abort.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee`
