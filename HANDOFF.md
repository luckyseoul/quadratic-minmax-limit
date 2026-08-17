> **LONG-HORIZON GOAL (binding):** settle MO 413935 — prove \(L=\lim\alpha_n\) exists and give its value, or prove non-existence. **Not done** until a terminal state in `LONG_HORIZON_GOAL.md`. Structure increments, honest OPEN, and “new prop shipped” are not completion.  
> **Now (2026-08-17 evening):** Leftover 1 still **OPEN**. 15.453–15.467 shipped. Coarse \(L_{ns}=\mu_+\mu_-\) at \(p=5,7\); named atoms miss \(L_{++}\); \(D\equiv3\pmod p\) dies at \(p=3\). Type-Fourier / \(S_k\) GJ / occupancy are dead. p=11 6-net ETA \(\ge580\) h — leave PID 712954, do not wait. Do **not** import \(\phi_F\). Residual (ii) \(k\ge4p\) and Type I multi-level still OPEN. Aut-Schur **false**. Live `e1` is the old AND. \(L=\tfrac12\) is **not settled**.  
> Afternoon: **`evidence/SESSION_HANDOFF_2026-08-17_qtau.md`**. Morning: **`evidence/SESSION_HANDOFF_2026-08-17_leftover.md`**. Status: **`STATUS.md`**. Package: **`evidence/share/denseness_path_package.md`**. Next `/goal`: **`GOAL.md`**.

# Research handoff: min-max ±1 quadratic form limit

**Status date:** 2026-08-17  
**Claim:** sandwich and Paley \(\rho=1\) are proved. **\(L=1/2\) is not.** E(1) still needs \(\lambda_{\min}(\Phi)\ge6\) (name \(Q_\tau\)), residual (ii) for even \(k\ge4p\), multi-level Type I, and Lemma D left checkable. Residual (ii) is closed only for affine + even \(k\le4p-2\). Aut-Schur remains **dead**. Gsum unused. residual/16N optional open.

**Problem:** [MO 413935](https://mathoverflow.net/questions/413935) · [X thread](https://x.com/PI010101/status/2081070728422752329).  
@PI010101 asked for GitHub, then posted that Grok Build 4.5 found a gap:  
https://grok.com/share/c2hhcmQtNA_218425aa-c1d1-4263-a3ea-9114ef04cd9c  
His bar: pass an **“AI test”** (any reasonable AI calls the writeup essentially correct in 2–3 repeated prompts) **before** he human-verifies. User cannot use MO (new account / no rep) — channel is **X + GitHub**, not MO answer.

**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**GitHub:** https://github.com/luckyseoul/quadratic-minmax-limit (`main`)  
**Repo description:** sandwich + Paley \(\rho=1\) proved; \(L=1/2\) OPEN pending four E(1) leftovers. Aut-Schur rejected; Gsum unused.

**Shipped close:**  
- **15.167:** bi-tight majorization algebra. `src/e1_gmin_m4_prop15167.py`  
- **15.170 / 15.171:** residual (i) Type I via dual-eq (15.216←249←207←270←272); residual (ii) ND. `gsum_disj_lb_proved_general()=False` (unused).  
- **15.179 / 15.236 / 15.237:** residual (ii) affine + (ii-b) + (ii-a) ND. `residual_ii_full_closed()=True`.  
- **15.272:** \(k=1\cup k=3\) spans \(\mathcal W_{++}^0\). Live units: built WWᵀ, DFT even-\(c\) rank \(m-1\), \(2\sigma=2\), \(M_3\) enum, isotypic dim fill. Pairing \(1^\top K^{-1}v\) stays False / unused.  
- residual/16N: still **OPEN** optional (not required)

---

## What’s left

**Done (2026-08-06):** gap named (Gsum disj LB / 15.158); claim retracted OPEN; short package `evidence/share/denseness_path_package.md`; load-bearing tests; 3 cold AI passes → NOT_PROVED + honest OPEN.

**Open (blocks \(L=\tfrac12\)):** (1) \(\lambda_{\min}(\Phi)\ge6\) — name \(Q_\tau\) by a **character-sum** formula on 15.290 types (15.466–467: named L fills miss Q++, coarse \(L_{ns}=\mu_+\mu_-\), named atoms miss \(L_{++}\), \(D\equiv3\pmod p\) dies at \(p=3\)); (2) residual (ii) for even \(k\ge4p\); (3) Type I when Max− is not two-level; (4) Lemma D already flagged True — do not cascade. Live `e1` is still True (wiring only).  

**Non-required:** Path-C / \(16N\).  

**Do not right now:** residual/16N thrash · X chrome · soft-close · ship a 15.xxx whose predicates stay False · flip e1 on census \(p\le7\) · re-thrash the residual-(i) dead list · `/goal resume` after `/goal clear` (start a **new** `/goal`).

**Census (not a close):** \(\nu\equiv0\) now **proved** (15.268), not just census. \(\mu\) a CR-class function at \(p=5,7\). Exact \(\sum_y\prod y_i\): \(p=3\) \(\to4\kappa\) (\(\mu=\mu_{\mathrm{part}}\)); \(p=5\) \(\to8\kappa-2\varphi\) (\(\mu=f_4\)); \(p=7\) \(\to240\kappa-10\varphi+96\kappa\varepsilon 1_{\lvert\varphi\rvert=10}\) (CR 9 vs 10). Envelope holds at \(p=3,5,7\). Switching / far-sum / Comm-repair dual (needs ker=sc) as before.

**Dead this continue (2026-08-15):** Aut-Schur / Jacquet “every irrep has a U-invariant \(\Rightarrow\) PSL-span of \(k=3\) \(F=\mathcal W_{++}^0\)”. Witness: \(p=5\) \(k=3\) Veronese rank \(61/65\); \(\dim F\sim p^2/4\) vs \(\dim\mathcal W_{++}^0=n(n-6)/8\sim p^4/8\). Singer Gram PD (Fejer / \(\mu=0\) / Weil / DFT) is only on even-on-\(\Omega\) \(F\), not a residual-(i) close. Do not re-import `gplus_pd` as True on that AND.  
**Dead 2026-08-14:** inversion-T as a \(p=5\) hinge (\(N=0\) on 2/3 K-chars; use Veronese); two-term \(M\) coded as if \(\Omega=\)squares (false zeros at \(p\equiv1\pmod4\)); Stickelberger digit-sum only \(\sim1/3\) of K-chars; affine-only / inv-T-alone kernels (together they span at \(7\le p\le23\)); envelope \(\lvert\mu\rvert\le2/n\) (false at \(p=11\)).  
**Dead earlier continue:** 4×4 / 5-point / K5 Gram (only \(1-2/p\)); Gsum 6×6 (\(|\mu|\le1\)); affine-quadratic level sets; CR-class master+diamond LP; IP-valency \(K_4\); global \(t(p)\) mix; 15.267 \(g^{-1}\) Aut slot (corrected in 15.268); unsigned Per even with \(|\mu|\le1-2/p\); Per-row CS; \(|\delta|\le\mathrm{room}_\delta\) from \(\mu_{\mathrm{part}}\) (false at \(p=5\)).

**Earlier dead (still dead):** Fréchet / joint Fréchet; \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise; W0\(\in\mathrm{im}(\mathrm{Gsum})\); Type I+switching \(G=\pi(G)\); low-degree Farkas \(y\); 3-wise interpolant; \(\kappa_B\perp E_{\pm4p}\); matching-PSD; max-\(\lvert\mu\rvert\) at extreme \(\lvert\varphi\rvert\); halfspace orbit as full \(\mu\); Aut-line dim\(\le1\); 15.237 pair-span as Type I kill.

**Hinge already shipped:** 15.272 \(k=1\cup k=3\) span of \(\mathcal W_{++}^0\) (not \(k=3\)-only \(F\)). Do not reopen Aut-Schur or \(N(\varphi)\neq0\). Envelope / reflection / \(\|m_4\|_2^2\le n(n-2)/4\) / \(K_4\le\mathrm{Wick}_{hi}\) / \(\lvert\mu\rvert\le2/n\) stay **dead**.

**Prior resumes:** `evidence/SESSION_HANDOFF_2026-08-14_residual_i.md`, `evidence/SESSION_HANDOFF_2026-08-13_residual_i.md`, `evidence/SESSION_HANDOFF_2026-08-12_residual_i.md`  
**P0 graph:** `evidence/P0_ENGINEERING_GRAPH.md` · **Failure graph:** `evidence/E1_FAILURE_GRAPH.md`  

**Suggested skills:** `agent-cost-optimization` · `graph-engineered-completion` · `scientific-critique` · `grill-me` · `self-refine-loop` · `handoff` · `session-handoff-packager` · `verification-before-completion` · `use-available-compute` · `research` · `arxiv` · `litreview` · `openai-referee`  · `claude-referee`

**Referee:** `claude-referee` (opus) **first**, then `openai-referee` on the **same** slot. Branch only on BLOCK. Do not run the two houses in parallel.  

**Compute:** 88 cores + idle V100; ProcessPool/GPU for Max+ contractions; no class_key thrash (F19); never soft-close L (F3). Token budget not a constraint this week — start a **new** `/goal` after `/goal clear` (cannot resume). Do not leave orphaned stdin `python3` eigsh jobs.

---

## 0. One-line status

Residual (ii) closed for affine + even \(k\le4p-2\), **OPEN** for even \(k\ge4p\). Floor leftover 1: \(Q_\tau\) still unnamed in \(p\) (15.467: name \(L_{++}\); ns average is \(\mu_+\mu_-\)). Multi-level Type I is OPEN. Aut-Schur stays false. Gsum unused. \(L=\tfrac12\) is **OPEN**.

### P0 status

**Bi-tight majorization (15.167):** mult≥d−1 + λ_min≥6 ⇒ L_* < 2d (algebra). In the E(1) AND at \(p=5\).

**E(1) not closed:** live `e1_closed_general` is still the old AND (wiring). Writeup must not treat that as a proof. Residual **(ii)** affine + \(k\le4p-2\) is ND; \(k\ge4p\) is open. Residual **(i)** two-level Type I via 15.272 (not `gsum_disj_lb`). **\(L=1/2\) is OPEN.**

**Hinge (every prime \(p\ge5\)):** 15.272 Johnson same-line hyperplane (live WWᵀ + DFT rank \(m-1\) + \(2\sigma=2\)) + \(k=3\) bad-\(\mu\) / complementary mixed / through-\(L_0\) + Singer PD on \(F\) for \(p\ge7\) + \(p=5\) rank 65. Pairing \(1^\top K^{-1}v\) unused. 15.271 \(k=3\)-only `fperp` stays False.

**Optional residual/16N (still OPEN, not required for denseness L):**
\[
E[s^4]\le \mathrm{Es4}_*(p)=\frac{4(3p^8-6p^6-148p^4-10p^2+129)}{p^4-8p^2-1}
\]
⇒ λ_max≤16 ⇒ 16N. Census p=5,7 only. Path-C form \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) also optional open.


**Do not re-derive:** 15.98 mult≥d−1; 15.102 δ-calculus; 15.110–15.116 reductions.  
**Dead for closing residual:** class_key T at p=7 (F19); unrestricted CS/LP on ED4; full \(E_{4p}\) energy of \(f_y\); soft-close from sandwich+denseness (F3).

**Modules:** `src/e1_gmin_m4_prop15110.py` … `prop15149.py` · evidence JSON same stems · tests `tests/test_prop1511*.py` `tests/test_prop1512*.py` `tests/test_prop1513*.py`  
**Max± caches:** `/tmp/maxplus_p5.npy` (260×26), `/tmp/maxminus_p5.npy`, `/tmp/maxplus_p7.npy` (11452×50), `/tmp/maxminus_p7.npy`, `/tmp/e1_p7/maxplus.npy`.
**Props 15.45–15.51 (load-bearing, do not re-derive):**
- Stars never bi-tight (wedge \(G^++G^-=0\)).
- Bi-tight level-2 empty if \(g_{\min}>-(p-2)/(p(2p-1))\) (Prop 15.47). **Certified \(p=5,7\)**.
- 1-bit spike (15.46). Deep at \(p=5\): small-\(k\) infeasible; large-\(k\) **spike** above \(\Phi(C)\).
- **Prop 15.48–15.49:** edge algebra; full CR classification of \(g_{\min}\) (cert p=3,5,7 in `e1_gmin_cr_classify.json`); **uniform LB candidate** \(g_{\min}\ge-(p-2)/(2p^2)\) algebraically beats the bi-tight threshold for all odd \(p>2\) and holds at certified \(p=5,7\) (`e1_gmin_uniform_lb.json`).
- **Prop 15.50 (proved):** Max+ conditional means given two coordinates equal the Gaussian frame interpolant \(\Sigma_{*S}\Sigma_{SS}^{-1}y_S\); conditional \(\mathbb E[y_ky_l\mid y_i,y_j]=\alpha+\delta y_iy_j\); disj average \(G=1/(p^2-2)\). Fréchet on cond cov **too weak** for \(L(p)\) (do not reopen). Evidence: `e1_gmin_cond_mean.json`.
- **Prop 15.51 (proved equiv form):** \(a=(1+pG)/(p+1)\); \(g_{\min}\ge T(p)\Leftrightarrow\min a(e')\ge1/(2p-1)\); disj sum deterministic on slice; residual Loewner cert p=5,7. Evidence: `e1_gmin_structure.json`.
- **Prop 15.52 (proved sum + moduli sketch):** \(\mathbf1^\top y=(p+1)y_\infty\) on Max+; m4 evec system sketch.
- **Prop 15.53 (proved pairing + cert moduli pin):** \(g_{\min}=-\max|m_4|\) on \(|\kappa|=1\); refined classes (37 at \(p=5\)); nullity-1; \(\mathrm{Tr}(G^2)\) pin recovers \(-3/65\); \(\mathbb E[\mathrm{dot}^2]\) from 2-design **proved**. `src/e1_gmin_moduli.py`.
- **Prop 15.54 (moduli \(c\)-calculus + wedge + a-slice + deep refresh):** wedge \(G=\pm1/p\); \(g_{\min}(c)\) pin; a_min p=5,7; deep covers spike at p=5.
- **Prop 15.55–15.56 (tight obstruction + star/cycle/Schur reduction):** If \(\lambda_{\max}(G)=n/2\) simple ⇒ no Max+-tight size-\(2p\) cover. **Proved:** stars map to \(\mathbf1\); \(\lambda_{\max}=\max(n/2,\lambda_{\mathrm{cycle}})\); \(\lambda_{\mathrm{cycle}}=2N\lambda_2(P\odot P)\); avg cycle eig \(<n/2\) for \(p\ge5\). **Certified gap** \(\lambda_2(P\odot P)\le d/(2N)\) at \(p=5,7\) (fails \(p=3\)). Evidence: `e1_gmin_spectral.json`.
- **Prop 15.57 (triangle bound):** \(|\mathrm{ft}|\le2p\), \(\mathbb E[\|By\|^2]\le4\), crude \(\lambda_{\max}(G)\le n\) (factor-2 short of \(n/2\)).
- **Prop 15.58 (Veronese residual form; 2026-07-30):** **Proved:** Max+\(\subset V_+\) (\(Cy=py\)); \((P\odot P)\mathbf1=\alpha\mathbf1\), \(\lambda_{\max}(P\odot P)=\alpha\) (Perron); gap \(\Leftrightarrow\|T(x)\|_F^2\le nN\|x\|^2\) for \(x\perp\mathbf1\) with \(T=\sum x_a y_ay_a^\top\); equiv. \(\lambda_2(W)\le\lambda_1(W)/2\) for \(W_{ab}=(u_a\cdot u_b)^2\). **Certified maximiser:** \(\mathrm{ft}=2p\), \(By\in V_+\), \(\|B_+\|_F^2=2\). Evidence: `e1_gmin_gap_probe.json`.
- **Prop 15.59 (center / rank / two-moment; 2026-07-30):** **Proved:** Max+ centrally symmetric \(\Rightarrow\sum y_a=0\Rightarrow P\mathbf1=0\). **Certified:** \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\) at \(p=3,5\); two-moment on \(W\) forces gap at \(p=7\) only. Evidence: `e1_gmin_veronese.json`.
- **Prop 15.60 (antipodal → projective ENTF; 2026-07-30):** **Proved:** \(T(x)=T(s)\); gap \(\Leftrightarrow\lambda_2(W^{(\mathrm{proj})})\le m/(2d)\); 2×sphere \(\Rightarrow\) gap for \(d\ge6\). Evidence: `e1_gmin_projective.json`.
- **Prop 15.61 (16N \(\Rightarrow\lambda_{\mathrm{cycle}}\le8\Rightarrow\) gap for \(p\ge5\); 2026-07-30):** **Proved equivalences:** \(Q(B):=\sum_y(y^\top By)^2\le16N\|B\|_F^2\) (Tr0 on \(V_+\)) \(\Leftrightarrow\lambda_2(W)\le4N/d^2\Leftrightarrow\lambda_{\mathrm{cycle}}\le8\). **Proved algebra:** \(\lambda_{\mathrm{cycle}}\le8\) and \(p\ge5\) (\(d\ge13>8\)) \(\Rightarrow\lambda_{\max}(G)=n/2\) simple \(\Rightarrow\) bi-tight empty (Prop 15.55). **Certified:** equality \(Q=16N\), \(\lambda_{\mathrm{cycle}}=8\) at \(p=3\); strict at \(p=5\) (\(Q/(16N)=11/13\)), \(p=7\) (\(\approx0.66\)). Evidence: `e1_gmin_16n.json`.
- **Prop 15.62 (typeA+wedge \(=6N\); \(Q=6N\|B\|_F^2+Q_4\); 2026-07-30):** **Proved:** for zero-diag \(B=P_+BP_+\), same-edge contrib \(=2N\|B\|_F^2\), wedge contrib \(=4N\|B\|_F^2\) (via \(\Sigma=2P_+\) on \(V_+\)), hence \(Q=6N\|B\|_F^2+Q_4\) with \(Q_4=Be^\top\mathrm{Gu}_{\mathrm{disj}}Be\) only. **Equivalence:** 16N \(\Leftrightarrow Q_4\le10N\|B\|_F^2\). **Certified** multi-seed p=3,5,7 (`e1_gmin_typeA_wedge.json`); p=3 saturates \(Q\equiv16N\) on whole space. Unrestricted \(\|\mathrm{Gu}_{\mathrm{disj}}\|_{\mathrm{op}}\) too large.
- **Prop 15.63 (H\(\Rightarrow\)16N; 2026-07-30):** **Proved algebra:** \(H(p):=(p+2)^2/d\le5\) for primes \(p\ge3\) (eq only \(p=3\)). **Proved:** \(\mathrm{ray}:=Q_4/(2N)\le H(p)\) \(\Rightarrow Q_4\le10N\) \(\Rightarrow\) 16N \(\Rightarrow\) bi-tight empty for \(p\ge5\). **Certified:** H holds p=3,5,7 with equality at p=3 and p=5 maximiser (\(\mathrm{ray}=49/13\)). Evidence: `e1_gmin_q4_ub.json`.
- **Prop 15.64 (dual \(\Phi\) form; residual reformulation; 2026-07-30):** **Proved:** \(\max Q=N\lambda_{\max}(\Phi|_{\mathcal Z})\); Wick baseline 8; H \(\Leftrightarrow\) residual\(\le(p+1)(p+7)/d\). **Certified** gen.eig p=3,5,7.
- **Prop 15.65 (κ spectrum; clean \(\lambda_2(P\odot P)\le4/N\); boolean essential; 2026-07-30):** **Proved:** \(16N\Leftrightarrow\lambda_2(P\odot P)\le4/N\). **Certified** κ spectrum + random-projector counterexample. Evidence: `e1_gmin_cumulant.json`.
- **Prop 15.66 (zero-diag freeness; pairing residual criterion; 2026-07-30):** **Proved:** \(\Phi\)-maximiser on \(\mathrm{Sym}_0\) automatically has ambient zero diag (\(\lambda_{\max}(\Phi|_{\mathcal Z})=\lambda_{\max}(\Phi|_{\mathrm{Sym}_0})\)). **Proved algebra:** for \(p\ge5\), \(|m_4-\kappa/p^2|\le(p-4)/(2p^2)\) on all \(|\kappa|=1\) 4-sets \(\Rightarrow g_{\min}\ge L(p)>T(p)\). **Certified:** \(g_{\min}\ge L\) at p=5,7; residual triangle holds p=7, **fails p=5** (large resid on small-\(|m_4|\) classes). Exact \(g_{\min}=-3/65,-109/2863\). Evidence: `e1_gmin_m4_residual.json`. **OPEN residual:** prove \(|m_4|\le(p-2)/(2p^2)\) on \(|\kappa|=1\) (or \(\lambda_2(P\odot P)\le4/N\)) for all primes \(p\ge5\). Bi-tight **not** closed for general \(p\).
- **Prop 15.67 (master m4 identity; same-sign Ext; census p=5,7; 2026-07-30):** **Proved:** \(\sigma_{sum}=4\kappa\) (K4 exhaustion); \(m_4=\kappa/p^2+\mathrm{Ext}/(4p)\) with \(\mathrm{Ext}=Tm_4\); same-sign \(|\mathrm{Ext}|\le2(p-4)/p\Rightarrow|m_4|\le L_{abs}\). **Certified W=86:** full \(|\kappa|=1\) census p=5,7 both satisfy direct \(|m_4|\le L\) and Ext criterion. Evidence: `e1_gmin_m4_proof.json`, `src/e1_gmin_m4_proof.py`, `src/workers.py`. **OPEN residual:** same as 15.66.4 for general primes \(p\ge5\).
- **Prop 15.68 (Tκ calculus; residual source on |κ|=3; resolvent reduction; 2026-07-30):** **Proved (conference C² + K4 exhaustion):** Tκ=0 on every |κ|=1 4-set; Tκ∈{±24} on |κ|=3; residual RHS Tκ/p² supported only on |κ|=3. **Proved algebra:** |m4|≤L ⇔ |ρ|≤(p-4)/(2p²) (same-sign) ⇔ resolvent gain |κ|=3→|κ|=1 ≤(p-4)/48; candidate |m4|≤(p-2)/(p(2p+3))≤L (sharp p=5). **Certified:** Paley degrees d3=p²-5, d1=3p²-7 at p=3,5,7. Evidence: `e1_gmin_m4_tkappa.json`. **OPEN:** resolvent gain / candidate for all primes p≥5.
- **Prop 15.71 (κ-stratum counts for any conference matrix; 2026-07-29):** **Proved:** ∑_S κ(S)² = n(n-1)(n-2)(n-5)/8 via C² wedge sum Σ=-n(n-1)(n-2) + K4 ratio-8 exhaustion ⇒ n1=n(n-1)(n-2)²/32, n3=n(n-1)(n-2)(n-6)/96 (Max+-free). **Certified W=86:** full κ census p=3,5,7,11 matches; d1,d3 constancy samples p=3,5,7,11. Evidence: `e1_gmin_m4_stratum.json`, `src/e1_gmin_m4_stratum.py`. Pins resolvent source size (# of |κ|=3 coords). **OPEN residual:** still |m4|≤M_mid/L_abs for all p≥5; lim α_n OPEN.
- **Prop 15.72 (resolvent-gain calculus; reverse degrees; type6 pin; 2026-07-29):** **Proved:** gain≤(p-4)/48 ⇔ same-sign |ρ|≤(p-4)/(2p²) ⇔ |m4|≤L on |κ|=1; Tκ/κ∈{±8} on |κ|=3; reverse degrees d1^(3)=3(p²-1), d3^(3)=p²-9 from n1,n3,d3 handshaking. **Certified W=86:** separate κ-weighted vanishing sk1=sk3=0 on all |κ|=1 at p=3,5,7; reverse deg constancy; type6 Max+-free resolvent |m|≤L and gain≪budget at p=5,7; true Max+ gains 0.0064, 0.036 < budget. Evidence: `e1_gmin_m4_resolvent_gain.json`. **OPEN:** prove gain≤(p-4)/48 for all primes p≥5; lim α_n OPEN.
- **Prop 15.73 (Max+ e4 identity; Paley ∑κ; multi-prime type6; 2026-07-29):** **Proved:** e4=-p(p-1)(p+1)(p+4)/12 from Prop 15.52 + boolean expansion ⇒ ∑_S m4(S)=e4. **Certified:** ∑κ=p²(p²-1)/4 at p=3,5,7,11,13; type6 Max+-free resolvent le_L and le_mid at p=5,7,11,13 with gain≪budget (type6 not exact upper bound — p=7 true m4>type6). Evidence: `e1_gmin_m4_e4_gain.json`. **OPEN residual:** true Max+ |m4|≤M_mid/L for all p≥5; lim α_n OPEN.
- **Prop 15.74 (candidate bound; true Max+ census; signed residual; 2026-07-29):** **Proved:** M_cand algebra; gain_cand=(p²-4p-3)/(24(2p+3)); 4p r=κ(Tρ) on |κ|=1. **Certified W=86 true Max+ (mmap, not type6):** max|m4|≤M_cand at p=5 (sharp 3/65) and p=7 (109/2863<5/119); same-sign gains 1/156 and ≈0.036. Evidence: `e1_gmin_m4_kernel.json`. **OPEN:** prove max|m4|≤M_cand for all primes p≥5; lim α_n OPEN.
- **Prop 15.75 (one-center σ_a=2·star_a; K4 Gram spectrum; GPU cand; 2026-07-29):** **Proved (any conference):** σ_a=∑_r C_ar κ(S_{a→r})=2·star_a on |κ|=1 via C²+Tκ=0; one-center residual form; K4 local G eigenvalues in {1±m4}∪{1±m4±2/p}∪{1±m4±2√2/p} ⇒ weak |m4|≤(p-2)/p. **Certified:** σ census p=3,5,7,11; GPU CuPy/V100 full |κ|=1 m4+cand (mmap Max+, device reduce, atomic JSON) p=5,7 le_cand in ~0.3s. Evidence: `e1_gmin_m4_onecenter.json`, `e1_gmin_m4_gpu.json`. **OPEN residual:** bound ∑C_ar ρ ⇒ |m4|≤M_cand all p≥5; lim α_n OPEN.
- **Prop 15.76 (one-center degrees d1^(1),d3^(1); residual split; 2026-07-29):** **Proved:** (3p²-7)/4 and (p²-5)/4 integral for odd p; 4 d1^(1)=d1. **Certified W=86:** one-center degrees constant on every |κ|=1 set/centre at p=3,5,7,11. **Proved:** abs bootstrap 4pρ≤d1 ρ+d3 R3 fails (4p−d1<0). **GPU residual:** R1,R3 on κ1/κ3 strata p=5,7; still le_cand. Evidence: `e1_gmin_m4_onecenter_deg.json`. **OPEN:** signed S1+S3 bound for M_cand all p≥5; lim α_n OPEN.
- **Prop 15.77 (star·S1≤0 structure; joint cand criterion; 2026-07-29):** **Proved:** star·(S1+S3)=pρ star−2/p²; star=+1 & S1≤0 ⇒ pρ≤2/p²+S3; same-sign ρ=(2/p²+S1+S3)/p at star=+1; cand ⇔ max joint ≤ p ρ_cand−2/p² (negative at p=5 ⇒ need strongly negative S1). **Certified GPU mmap+atomic:** star·S1≤0 on **all** |κ|=1 centres p=5,7 (strict); joint ⇒ ρ≤ρ_cand (sharp p=5); synthetic non-Max+ bump violates star·S1. Evidence: `e1_gmin_m4_S1_star.json`, `src/e1_gmin_m4_S1_star.py`. **OPEN:** prove star·S1≤0 + joint/S3 bound for all primes p≥5; lim α_n OPEN.
- **Prop 15.78 (moment form/GD; 4-set constancy; p=5 exact spectrum; 2026-07-29):** **Proved:** star·S1=E[φ]−E_Wick[φ] (φ=star f0 U1); GD ⇔ star·S1≤0; joint rewrite S1+S3=p m4−1/p−2/p² at κ=1,star=+1. **Certified W=86 pure-C:** star·τ1 constant on every |κ|=1 set p=3,5,7,11. **Certified GPU mmap+atomic:** star·S1 constant on every |κ|=1 Max+ set p=5,7; GD holds; star·S1≤0; **p=5 exact** star·S1∈{−2/65,−42/325} both <0 (⇒ star·S1≤0 **proved at p=5** by full census). Evidence: `e1_gmin_m4_S1_const.json`, `src/e1_gmin_m4_S1_const.py`. **OPEN:** prove GD for all p≥5; joint/S3 ⇒ M_cand; lim α_n OPEN.
- **Prop 15.79 (Aut-constancy of star·τ1 and star·S1; modular τ1; 2026-07-29):** **Proved:** PGL(2,p²) setwise stabilizer of any 4-set contains V4 acting regularly ⇒ transitive on the 4 points; star·τ1 Aut-equivariant ⇒ **constant on centres** (Max+-free); star·S1 Aut-equivariant (Max+ Aut-invariant) ⇒ **constant on centres**; τ1=2A−d1^(1) with d1 odd ⇒ star·τ1 odd. **Certified W=86:** constancy + t1=2A−d1 + star·τ1≡5 (mod 6) + #(values)=(p−1)/2 at p=3,5,7,11. Evidence: `e1_gmin_m4_S1_aut.json`, `src/e1_gmin_m4_S1_aut.py`. **OPEN:** prove star·S1≤0 (GD) all p≥5; joint/S3⇒M_cand; lim α_n OPEN.
- **Prop 15.80 (linear-form Wick; GD ⇔ E[ZU1]≤E_Wick[ZU1]; U1-special; 2026-07-29):** **Proved:** E[L²]=E_Wick[L²] for every linear L on Max+ (pairwise Σ only); GD ⇔ Cov(Z,U1)≤Cov_Wick with Z=star f0; GD fails for generic L (U1 κ1-support essential). **Certified GPU mmap+atomic p=5,7:** GD holds all |κ|=1 sets; U1/Wick ratio≡1; sum star·S1=−1128 (p5), −15271200/2863 (p7); sum star·τ1=ε n1 (ε=±1); generic L violation rate ~45–55%. Evidence: `e1_gmin_m4_S1_gd.json`, `src/e1_gmin_m4_S1_gd.py`. **OPEN:** prove E[ZU1]≤E_Wick all p≥5; joint/S3⇒M_cand; lim α_n OPEN.
- **Prop 15.81 (moduli-line GD criterion; p=5 complete; 2026-07-29):** **Proved form:** on nullity-1 line m=m★+c n, if max star·S1=α+βc with β>0 then GD⇔c≤c_GD=−α/β. **Certified p=5:** 37 constant-m4 classes, nullity 1, max star·S1 exact affine (β>0), c*≈−0.424<c_GD≈−0.296, max|m4|(c*)=M_cand=3/65, GD holds — **cand+GD at p=5 via moduli**. **Certified p=7:** pointwise evec holds; coarse classes 69/82 constant m4 (incomplete for nullity-1 line). Evidence: `e1_gmin_m4_S1_moduli.json`, `src/e1_gmin_m4_S1_moduli.py`. **OPEN:** refine classes all p≥5; prove c*≤c_GD generally; lim α_n OPEN.
- **Prop 15.82 (type6+CR refinement; m4 constant p=5,7; 2026-07-29):** **Certified GPU+W=86:** type6+CR makes m4 constant on all classes at **p=5 (26)** and **p=7 (48)** (coarse alone incomplete at p=7); coarse+CR also fully constant (48 / 130). **Moduli type6+CR:** p=5 nullity 1, c* safe-side of c_GD (β may be negative — orient by sign), cand+GD; p=7 nullity 2, true Max+ GD and max|m4|≤M_cand, multi-param pin OPEN. Evidence: `e1_gmin_m4_refine.json`, `e1_gmin_m4_refine_moduli.json`; `src/e1_gmin_m4_refine.py`, `src/e1_gmin_m4_refine_moduli.py`. **OPEN:** nullity≤1 or multi-param pin for all p≥5; lim α_n OPEN.
- **Prop 15.83 (resolvent-budget hierarchy \(M_{\mathrm{cand}}\) vs \(L\); Max+-free; 2026-07-30):** **Proved algebra:** \(M_{\mathrm{cand}}<M_{\mathrm{mid}}\le L_{\mathrm{abs}}<T_{\mathrm{abs}}\) for \(p>2\); \(0<\rho_{\mathrm{cand}}<\rho_L\); \(\mathrm{gain}_L-\mathrm{gain}_{\mathrm{cand}}=3(p-2)/(48(2p+3))>0\) so \(\mathrm{gain}_{\mathrm{cand}}<\mathrm{gain}_L\) for all primes \(p\ge5\). Thus any resolvent-gain bound \(\le\mathrm{gain}_{\mathrm{cand}}\) yields \(\max|m_4|\le M_{\mathrm{cand}}\le L\) and bi-tight empty. **OPEN residual:** prove \(\mathrm{gain}\le\mathrm{gain}_{\mathrm{cand}}\) (or \(\le\mathrm{gain}_L\)) for true Max+ for all \(p\ge5\). Evidence: `src/e1_gmin_m4_prop1583.py`, `evidence/e1_gmin_m4_prop1583.json`. **lim α_n OPEN.**
- **Prop 15.84 (GD⇒cand via \(S_3\) budget; diag-dom fails; Max+-free; 2026-07-30):** **Proved:** \(B_{\mathrm{cand}}=(p^3-4p^2-7p-6)/(p^2(2p+3))\); \(B_{\mathrm{cand}}(5)=-16/325<0\); \(B_{\mathrm{cand}}(p)>0\) for primes \(p\ge7\); GD+\(S_3\le B_{\mathrm{cand}}\)⇒cand; \(4p-d_1=-3p^2+4p+7<0\) (no abs bootstrap). **OPEN:** prove GD and \(S_3\le B_{\mathrm{cand}}\) for Max+ ∀p≥5. Evidence: `src/e1_gmin_m4_prop1584.py`, `evidence/e1_gmin_m4_prop1584.json`. **lim α_n OPEN.**
- **Prop 15.85 (\(Q_4\) mean/fluctuation split; Path C spectral; Max+-free; 2026-07-30):** **Proved:** \(\mu=e_4/C(n,4)\); \(S_1=0\) for zero-diag \(B=P_+BP_+\); unit-F \(\mathrm{ray}=-\mu-(\mu/2)S_w+\tfrac12 Be^\top\widehat G Be\); \(|\mu|\ll H(p)\); Gershgorin on \(\widehat G\) too weak. **OPEN residual:** control \(\widehat G\) / hypothesis H / \(\max|m_4|\) for ∀p≥5. Evidence: `src/e1_gmin_m4_prop1585.py`, `evidence/e1_gmin_m4_prop1585.json`. **lim α_n OPEN.**
- **Prop 15.86 (Wick mean \(\varepsilon(p)\); \(\tau_1\) spectrum; residual budgets; Max+-free; 2026-07-30):** **Proved form + certified W=86 pure C:** \(\sum_{|\kappa|=1}\mathrm{star}\cdot\tau_1=\varepsilon(p)\,n_1\) with \(\varepsilon(p)=(-1)^{(p-1)/2}\) (closes Prop 15.80.4 OPEN); mean Wick \(=\varepsilon/p^2\); GD⇒mean \(\mathrm{star}\cdot S_1\le0\); \(\tau_1\) value set AP of diff 6, size \((p-1)/2\), all \(\equiv5\pmod6\) (exact sets p=3,5,7,11); \(B_{\mathrm{cand}}\to1/2\), \(B_{\mathrm{cand}}/d_3=\Theta(1/p^2)\). **OPEN:** pointwise GD / \(S_3\le B_{\mathrm{cand}}\) / max|m4|≤M_cand / H for ∀p≥5. Evidence: `src/e1_gmin_m4_prop1586.py`, `evidence/e1_gmin_m4_prop1586.json`. **lim α_n OPEN.**
- **Prop 15.87 (K4 star; \(S_1\) pattern; GD reformulation; \(\mathbb E[U_1^2]\); 2026-07-30):** **Proved:** on every \(|\kappa|=1\) K4 labeling (48/64), \(\sum\mathrm{star}=0\), \(\prod\mathrm{star}=+1\), exactly two \(+\) stars; hence \(S_1(a)=g\cdot\mathrm{star}_a\) and GD\(\Leftrightarrow g\le0\); residual \(4p\rho=T\rho\) is tautological on \(\kappa1\) given \(\sum\mathrm{star}=0\); CS \(|E[ZU_1]|\le\sqrt{E[U_1^2]}=\Theta(p)\) too weak vs Wick scale. **Certified W=86:** \(E[U_1^2]\) near \(d_1\) at p=3,5,7. **OPEN:** pointwise \(g\le0\) for ∀p≥5. Evidence: `src/e1_gmin_m4_prop1587.py`, `evidence/e1_gmin_m4_prop1587.json`. **lim α_n OPEN.**
- **Prop 15.88 (pairwise sum; H-gap algebra; \(g\) via \(S_3\); 2026-07-30):** **Proved:** \(\sum_{i<j}y_iy_j=p\) on Max+ (from Prop 15.52); \(n/2-(3+H)=(p^4-8p^2-16p-21)/(2(p^2+1))>0\) for primes \(p\ge5\) (so H\(\Rightarrow\) bi-tight empty); \(g=p\rho\,\mathrm{star}-2/p^2-\mathrm{star}\cdot S_3\) (GD\(\Leftrightarrow\) star·S3 lower bound). **Certified p=5 full Max+:** \(g=g(\tau_1)\in\{-42/325,-2/65\}\). **OPEN:** H or pointwise \(g\le0\) for ∀p≥5. Evidence: `src/e1_gmin_m4_prop1588.py`, `evidence/e1_gmin_m4_prop1588.json`. **lim α_n OPEN.**
- **Prop 15.89 (Wick split of \(Q_4\); \(\kappa_C\cdot\kappa_B\); H residual form; 2026-07-30):** **Proved:** \(Q_4/N=8\sum m_4\kappa_B\); \(\sum\kappa_C\kappa_B=(n+1)/4\|B\|_F^2\) (cert p=3,5,7,11, Max+-free); unit \(Q_4/N=2+4/p^2+8\sum\rho\kappa_B\); H \(\Leftrightarrow\sum\rho\kappa_B\le(H-1-2/p^2)/4\). Wick part closed. **OPEN residual:** bound \(\sum\rho\kappa_B\) (Max+ \(\rho\)) for ∀p≥5. Evidence: `src/e1_gmin_m4_prop1589.py`, `evidence/e1_gmin_m4_prop1589.json`. **H not yet proved. lim α_n OPEN.**
- **Prop 15.90 (residual bound ≡ H; pointwise κ_B; orth form; 2026-07-30):** **Proved:** ∑ρ κ_B ≤ budget ⇔ ray≤H (same inequality). Pointwise ∑κ_B∏y = f²/8 − γ/2 + ‖B‖²/4 (cert p=3,5). Orth form E[‖By−(f/n)y‖²]≥2−(6+2H)/n ≡ H. Bound holds p=3,5 (eq), p=7 (strict). **Not an independent foothold.** Evidence: `src/e1_gmin_m4_prop1590.py`, `evidence/e1_gmin_m4_prop1590.json`. **H_proved=false. lim α_n OPEN.**
- **Prop 15.91 (independent dual forms of H; dim Z; sphere/harm; 2026-07-30):** **Proved:** dim Z=d(d−3)/2 (conference ETF Gram full rank). Orth / Φ / κ / harm forms ≡ H (no ∑ρ κ_B). Chain sphere<Wick<6+2H≤16 ∀p≥3. 2×sphere⇒16N/bi-tight for p≥5 (weaker than H). Evidence: `src/e1_gmin_m4_prop1591.py`, `evidence/e1_gmin_m4_prop1591.json`. **H_proved=false. Attack: orth LB or λ_max(κ|Z)≤(p+1)(p+7)/d or 2×sphere. lim α_n OPEN.**
- **Prop 15.92 (constant pairing sum; H/16N spectral reductions; 2026-07-30):** **Proved:** ∑_S κ∏y = n(n−1)(n−2)/8 constantly on Max+ ⇒ ∑ m4 κ_C = n(n−1)(n−2)/8. Reductions: 16N⇔λ₂(P⊙P)≤4/N; H⇔λ₂(P⊙P)≤(3+H)/(2N). W spectrum saturates H at p=3,5 (mult d top). Evidence: `src/e1_gmin_m4_prop1592.py`, `evidence/e1_gmin_m4_prop1592.json`. **H_proved=false. Attack: λ₂(P⊙P)≤4/N (16N) or H form. lim α_n OPEN.**
- **Prop 15.93 (Gu/FFT structure; 16N as λ_max(FFT|1⊥)≤8N; 2026-07-30):** **Proved:** FFT1=Nd1; Gu has eig Nd on C-edge vector; 16N⇔λ_max(FFT|_{1^⊥})≤8N⇔λ_max(Gu on Im)≤8N; H⇔≤N(3+H). Gu spectrum cert p=3 (16N eq), p=5 (H eq, 16N strict); non-Nd eigs =(N/2)·spec(Φ|Z). Evidence: `src/e1_gmin_m4_prop1593.py`, `evidence/e1_gmin_m4_prop1593.json`. **H_proved=false. Attack: λ_max(FFT|1⊥)≤8N ∀p≥5. lim α_n OPEN.**
- **Prop 15.94 (P⊙P annihilates range(P); gap criterion mult≥d; 2026-07-30):** **Proved:** P⊙P kills range(P) (third moments vanish by central sym); ∑P³=0. Gap criterion: mult(λ₂)≥d and ∑M²≤4d²(d+4) ⇒ λ₂≤d/(2N) ⇒ bi-tight empty for p≥5. Cert p=5: mult=d and criterion holds for gap (not 16N); p=3 correctly fails gap. Evidence: `src/e1_gmin_m4_prop1594.py`, `evidence/e1_gmin_m4_prop1594.json`. **OPEN: mult≥d ∀p≥5 and ∑M² bound (or 16N/H). lim α_n OPEN.**
- **Prop 15.95 (Wick≤thr; strengthened gap; C_diag; 2026-07-30):** **Proved:** Wick_hi=12n²+48n ≤ thr_gap=4d²(d+4) for all primes p≥5 (disc n²−16n−96>0); mult≥d + ∑M²≤Wick ⇒ gap ∀p≥5; C_diag=4n(11n−14)/p². **Certified** mult(λ₂)=d and gap_by_mult at p=5,7 (p=3 fails gap correctly); ∑M²≤thr at p=5,7. Evidence: `src/e1_gmin_m4_prop1595.py`, `evidence/e1_gmin_m4_prop1595.json`. **OPEN: mult≥d ∀p≥5 and/or ∑M²≤Wick (or 16N/H). lim α_n OPEN.**
- **Prop 15.96 (Wick–κ calculus; ‖κ‖²≤96n ⇔ ∑M²≤Wick; 2026-07-30):** **Proved:** yᵀΣy=2n; ⟨Wick,M⟩=12n²; ⟨Wick,κ⟩=−48n; ∑M²=12n²−48n+‖κ‖²; ∑M²≤Wick ⇔ ‖κ‖²≤96n; C_diag≤96n algebra. With mult≥d + ‖κ‖²≤96n ⇒ gap ∀p≥5. **Certified** p=3,5 (‖κ‖²=96n eq only p=3). Evidence: `src/e1_gmin_m4_prop1596.py`, `evidence/e1_gmin_m4_prop1596.json`. **OPEN: ‖κ‖²≤96n and mult≥d ∀p≥5 (or 16N/H). lim α_n OPEN.**
- **Prop 15.97 (Veronese mult ID; Ky Fan criterion; 2026-07-30):** **Proved:** mult(λ₂(P⊙P))=mult(λ_max(Γ|Sym₀)) via Veronese Gram; Aut-Schur; mult≥d ⇔ d orthonormal maximizers of Var(cᵀBc) on Sym₀ (Ky Fan). **Certified** mult=d and Ky Fan equality at p=3,5. Evidence: `src/e1_gmin_m4_prop1597.py`, `evidence/e1_gmin_m4_prop1597.json`. **OPEN: exhibit d maximizers for general p≥5; ‖κ‖²≤96n. lim α_n OPEN.**
- **Prop 15.98 (mult≥d−1 via PSL; gap with d−1; 2026-07-30):** **Proved for Paley Max+:** mult(λ₂)≥d−1 (PSL(2,p²) min irrep dim; PSL simple ⇒ faithful image). **Proved algebra:** mult≥d−1 + ‖κ‖²≤96n ⇒ gap for all primes p≥5. **Certified** p=3,5. Evidence: `src/e1_gmin_m4_prop1598.py`, `evidence/e1_gmin_m4_prop1598.json`. **OPEN for bi-tight: only ‖κ‖²≤96n ∀p≥5.** lim α_n OPEN (then deep ND).
- **Prop 15.99 (κ structure; min-distance; closed forms; 2026-07-31):** **Proved:** min $d_H\ge p+1$ for $y\ne\pm z$ in Max+; Wick/ρ budget algebra $\|\kappa\|^2\le96n\Leftrightarrow\sum\rho^2\le n(13p^2+3)/(6p^2)$; closed forms for $\sum\kappa,\sum\kappa^2,\sum m_4,\sum m_4\kappa,\sum\rho,\sum\rho\kappa$; master source $T\kappa$ on $|\kappa|=3$. **Certified** $\|\kappa\|^2\le96n$ at $p=3,5$ (eq only $p=3$). **OPEN:** $\|\kappa\|^2\le96n$ for all primes $p\ge5$. Evidence: `src/e1_gmin_m4_prop1599.py`, `evidence/e1_gmin_m4_prop1599.json`. **lim α_n OPEN.**
- **Prop 15.100 (dual-frame/flat≤Wick/κ_hyp; 2026-07-31):** **Proved:** $\|\kappa_{\mathrm{proj}}\|^2=64n(p^2-3)/(p^2-5)$; ED4_flat≤Wick; κ_hyp≤96n. **Certified** κ²=κ_hyp at p=3,5. **OPEN:** orth bound. Evidence: `src/e1_gmin_m4_prop15100.py`, `evidence/e1_gmin_m4_prop15100.json`.
- **Prop 15.101 (Fickus Gram residual / bulk variance; 2026-07-31):** **Proved:** PopP Fickus rank $\binom{d-1}{2}$; orth$=16N^2(\sum_{\mathrm{bulk}}\lambda^2-S^2/m)$; PSL level count; N-free ε-criterion (sufficient, not nec.). **Certified** bulk p=3,5. Evidence: `src/e1_gmin_m4_prop15101.py`.
- **Prop 15.102 (resolvent δ-calculus; 2026-07-31):** **Proved:** $T^2b=\mu^2b$ form $\mu^2=4(p^2+15)$; $\|\rho_{\min}\|_2^2=5n(p^2-1)(p^2+3)/(6p^2(p^2-5))$; $\kappa_{\min}=\mathrm{proj}$; $\|\kappa_{\mathrm{orth}}\|^2=24\|\delta\|_2^2$; target $\Leftrightarrow\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$. Evidence: `src/e1_gmin_m4_prop15102.py`.
- **Prop 15.103 (δ-bound cert p=3,5,7; 16N p=5,7; 2026-07-31):** **Certified** full Max+ $N=12,260,11452$: $\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24$ at $p=3,5,7$ (eq at 3,5; ratio $\approx0.187$ at 7). $\kappa^2=\kappa_{\mathrm{hyp}}$ only at $p=3,5$ (strict at 7). **16N** at $p=5,7$ ($\lambda_2(P\odot P)\le4/N$). Bi-tight empty at $p=5,7$ by 15.98+κ or 15.61+16N. **OPEN:** same δ-bound (or 16N) for all primes $p\ge5$. Evidence: `src/e1_gmin_m4_prop15103.py`, `evidence/e1_gmin_m4_prop15103.json`. **lim α_n OPEN.**

- **Residual to \(m_n\ge\Phi-2\Rightarrow L=\tfrac12\):** (1) prove $\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$ or $\lambda_2\le4/N$ for all $p\ge5$ (bi-tight); (2) deep non-tight ND. **F13/F16–F18** intact.

**g_min attack status (Prop 15.49–15.50):** \(g_{\min}=-\alpha_\star\) on constant-\(m_4\) \(|\kappa|=1\) CR classes. Values \(-\tfrac13,-\tfrac3{65},-\tfrac{109}{2863}\). Matching form \(\mathbf1_M^\top G\mathbf1_M\ge9.96>4\) at \(p=5\) (no tight matching). **Dead:** \(-3/\Phi\) general LB; 4-point LP; Chebyshev; Wick-as-LB; bare \(C\)-types; affine halfspace orbit; pure deg pigeon; min-norm \(V_+\) interpolation alone; Bose–Mesner alone (span too large); plain CLT residual (overestimates m4); **plain Fréchet on conditional cov** (Prop 15.50.4). **Need:** prove uniform LB for all \(p\ge5\) (character sum / bound \(\delta\) or \(m_4\) on min CR class) + deep residual.

**New (n=10 structure):** exact optima first appear at Hamming distance **5** from Paley \(C_{10}\), and the only 5-edge undercutters are **144 perfect matchings** (of 945). Absolute gap \(\Phi-m_{10}=2\) is consistent with E(1). See `evidence/N10_STRUCTURE.md`.

**New (n=10 classification, N10-C):** those 144 matchings are exactly one \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit, equivalently the matchings that drop every Paley maximizer to \(|Q|\le13\) (six \(+\) maximizers certify). See `evidence/N10_MATCHING_CLASSIFY.md`.

**New (E(2) interval formula):** for Paley \(q\equiv1\bmod4\), interval signing has exact
\(x^\top Cx=2-8\sum_{d\le(q-1)/2}d\chi(d)\) (`evidence/E2_INTERVAL_FORMULA.md`). Constructive
\(\rho_{\mathrm{int}}\) up to \(\approx0.978\) at \(n=1622\). **Not** a full proof of \(\rho\to1\).

**New (E(1) numerics):** exact-Φ SA at Paley \(n=14,18\): no undercut of \(\Phi(C)\) found
(`evidence/e1_paley_gap.json`); \(n=10\) recovers gap \(2\); \(n=14\) exact \(k\le4\) flips stay at \(\Phi=21\);
**\(n=26\) intensive** (172×10k SA + exact rescore + k-flips): **no undercut of \(\Phi=65\)**
(`evidence/e1_n26_intensive.json`, `E1_STATUS.md`). Maximizers of \(\rho=1\) Paley satisfy
\(\mathbb E[yy^\top]=I\) (2-design; \(p=3,5\)). **Not** a proof of E(1).

**Settlement path (shortest):** E(1) alone on the dense \(\rho=1\) family \(n=p^2+1\) \(\Rightarrow L=\tfrac12\) by Prop 6.2
(E(2) not needed for that subfamily). **E(1) is the blocking open problem.**

**E(1) reduction (2026-07-27):** Prop 15.20b edge-counting Lipschitz \(\Phi(A)\ge\Phi(C)-2k\) sharpens the sparse regime (vs Frobenius \(n\sqrt k\)). **E(1) \(\Leftrightarrow\) \(k_\star=o(n^{3/2})\)** for \(\Phi\)-minimisers after switching to Paley on \(n=p^2+1\). At \(n=10\), \(k_\star=5\) (matchings), though some other \(m_{10}\)-matrices have best-\(k\ge15\) (campaign). Rigidity \(k_\star=o(n^{3/2})\) still unproved in general. See `evidence/E1_EDGE_LIPSCHITZ.md`, `evidence/E1_RIGIDITY_ATTACK.md`.

**N10-C6 (2026-07-27):** exhaustive scan of \(\binom{45}{6}\) — all **360** Hamming-6 undercutters of Paley \(C_{10}\) are single **6-cycles** with \(\Phi=13=m_{10}\). Together with N10-S (144 matching undercutters at \(k=5\)), every undercutter at the two smallest cardinalities is a path/cycle graph (\(\Delta\le2\), \(k\le n\)). Sufficient for E(1) if generalised: **path-cycle dichotomy** \(k_\star\le n\) \(\Rightarrow L=\tfrac12\). Far optima need not *contain* a sparse undercutter as a subgraph; only existence of a sparse minimiser matters for \(k_\star\). See `evidence/N10_CYCLE_UNDERCUTTERS.md`, `src/n10_cycle_undercutters.py`.

**Prop 15.24 (2026-07-27):** maximizers of any ρ=1 conference are exactly the boolean ±p-eigenvectors. Paley counts of +p boolean evecs: 12, 260, 11452 at p=3,5,7 (`evidence/BOOLEAN_EVECS_MAX.md`). Ratio #/n^{3/2} increases, so k_⋆≤|Max| is not o(n^{3/2}).

**Prop 15.25 (2026-07-27):** recursive formula \(m_n=\min_B\max_x(|Q_B(x)|+|\sum x_i|)\) over Seidel \(B\) of order \(n-1\) (star-reduction). Certified for \(n\le11\). Does not by itself close E(1) (boost \(\max_L|s|\) often \(O(1)\)).

**Prop 15.26 (2026-07-27):** matching flips of any \(\rho=1\) conference **preserve local maximality** of boolean \(+p\)-eigenvectors (\(y_i(Ay)_i\ge p-2>0\)). At \(n=10\), for all 945 perfect matchings the stronger identity \(\Phi(C\oplus M)=\max_{\mathrm{Max}}|Q_{C\oplus M}|\) holds (cube max attained on boolean eigenset). Star / high-\(\Delta\) flips destroy this. Star-reduction probe: \(f(C')=\Phi(C)\), \(\Phi(C')=m_9\), but matching undercutters give \(d_H(B,C')=12\) (sparsity not preserved). See `evidence/E1_STAR_REDUCTION_PROBE.md`. Route to matching dichotomy still open.

**Prop 15.27 (2026-07-27):** for Paley \(\rho=1\), \(\mathrm{Max}_{+}\) is a tight frame: \(\mathbb E[yy^\top]=I+C/p=2P_{+}\). Consequences: (i) fractional Max-cover LP has value exactly \(p\); (ii) **Max-Lipschitz** \(\Phi(A)\ge\Phi(C)-2k/p\) after best switch (edge lip improved by factor \(p=\sqrt{n-1}\)). E(1) demand weakens from \(k_\star=o(n^{3/2})\) to \(k_\star=o(n^2)\) — still not free (worst-case \(k\sim n^2/4\)). Existence remains **OPEN**.

**Prop 15.28 (2026-07-27):** size-\(p\) Max-covers are tight (\(S_F\equiv1\) on \(\mathrm{Max}_{+}\)). **Proved for Paley:** \(\infty\)-stars whose leaves are nonsquare-direction affine lines of \(\mathrm{AG}(2,p)\) are tight Max-covers with \(\mathbb E_{-}[S^2]=2p-1\), hence \(\max_{\mathrm{Max}_{-}}S\ge1\) and \(\Phi\ge\Phi+2\) (cannot undercut). Square directions are non-covers. Certified: all covering \(p\)-stars (any centre) spike to \(\Phi+2p\) at \(p=3,5\); all 405 size-\(p\) covers at \(n=10\) have \(\Phi\ge17\). See `evidence/E1_SIZE_P_MAXCOVER.md`. **Does not close E(1)** (larger covers / \(k_\star\) still open).

**Prop 15.29 (2026-07-27):** (i) \(n/2\) odd \(\Rightarrow\) matching scores \(S_M\) always odd; non-covering PMs raise \(\Phi\) by \(\ge2\). (ii) Undercutters are two-sided Max-covers. (iii) **Correction:** Max-covering PMs **exist** at \(n=26\); the three SA-found covers are two-sided with Max\(\pm\) height \(63=\Phi-2\) but exact MITM \(\Phi=65=\Phi(C)\) (spike). 0 undercuts. See `evidence/E1_MATCHING_COVER_SPIKE.md`. **Does not close E(1).**

**Prop 15.30 (2026-07-27):** Matching spike criterion **proved**: \(S_M=-p\) \& \(Q_C\ge\Phi-2p\) \(\Rightarrow\Phi(C\oplus M)\ge\Phi(C)\). At \(n=10\), criterion holds iff matching does not undercut (complete over 945 PMs). At \(n=26\), all tested PMs satisfy criterion. Open: prove criterion for every PM when \(p\ge5\). See `evidence/E1_MATCHING_SPIKE_CRITERION.md`. **Does not close E(1).**

**Prop 15.31 (2026-07-27):** Clique-flip sufficiency for matching Max-covers **proved**. Only \(|F|=p\) works for covers; \(p=3\) undercutters block the arithmetic; \(p=5\) forces \(S_M=1\). Design: 390 Seidel \(p\)-sets, 60 Max\(_{+}\) extensions each, \(\ge236\) transversal per matching. All SA covers clique-flip to \(\Phi(C)\). Open: \((y,F)\) for every Max-cover matching, \(p\ge5\). See `evidence/E1_CLIQUE_FLIP.md`. **Does not close E(1).**

**Prop 15.32 (2026-07-27):** Γ-pairing reformulation of the spike criterion **proved**; coordinate product \(\pi\) constant on Max\(_{+}\) **proved**; \(S_M\bmod 4\) constant on Max\(_{+}\) for every matching **proved**. Case split: maximiser / 1-bit / clique-flip according to residue class. Census at \(p=5\): \(\max R\in\{60,70\}\) on all tested PMs (covers tight at 60); 0 counterexamples among random/SA-min-\(\max R\). Open: attainment of \(-p\) or \(-p+2\) in every residue class, and clique-flip on every cover. See `evidence/E1_GAMMA_PAIRING.md`, `e1_gamma_forall_census.json`. **Does not close E(1).**

**Prop 15.33 (2026-07-27):** **Non-covers cannot undercut** (proved): if \(\min_{\mathrm{Max}_{+}}S_F\le-1\) then \(\Phi(C\oplus F)\ge\Phi+2\). Matching non-undercut therefore reduces **only** to Max-covering matchings. Spike criterion is **not necessary**: certified matching with \(\max R=54<60\) but non-cover and MITM \(\Phi=75>\Phi\) (`e1_criterion_fail_no_undercut.json`). All SA Max-covers at \(p=5\) still satisfy criterion and \(\Phi=\Phi(C)\). Open: forall Max-cover matchings, \(p\ge5\). **Does not close E(1).**

**Prop 15.34 (2026-07-27):** Matching flip algebra \(D^2=I\), \(A=C-2D\), \(A^2=(n+3)I-2(CD+DC)\) **proved**. At \(p=5\), every tested Max-cover has \(\|A\|_{\mathrm{op}}=\sqrt{41}\) exactly; random matchings larger. Open: use spectral control of Max-cover flips to force \(\Phi(A)\ge\Phi(C)\). **Does not close E(1).**

**Prop 15.35 (2026-07-27):** At \(p=5\), every Max-cover PM **must** attain \(S_M=1\) with residue \(1\bmod4\) (proved by expectation). Census of **11** SA Max-covers: all two-sided, min+max as covers, \(\mathrm{op}=\sqrt{41}\), clique-flip, \(\Phi=\Phi(C)\). 0/20k random PMs are covers. Open: forall Max-cover \(M\) when \(p\ge5\). See `evidence/e1_maxcover_full_census.json`. **Does not close E(1).**

**Prop 15.36 (2026-07-27):** Matching flip **block algebra** proved: \(B=CD+DC\) always commutes with \(C,D\); \(B|_{V_+}=2p D_{++}\), \(B|_{V_-}=-2p D_{--}\); \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\). At \(p=5\), every SA Max-cover has \(\lambda_{\min}(B)=-6\) and \(\mathrm{op}=\sqrt{41}\), with **≥2** distinct \(D_{++}\) spectral types (both non-undercut, clique-flip, \(\max R=60\)). Random non-covers have more negative \(\lambda_{\min}(B)\). Open: prove \(\lambda_{\min}(B)=-6\Rightarrow\Phi\ge\Phi(C)\) for all Max-covers \(p\ge5\). See `evidence/E1_MAXCOVER_SPECTRUM.md`. **Does not close E(1).**

**Prop 15.37 (2026-07-27):** Continuous Γ-bound pattern: all 11 Max-covers have \(\min_{S=-p}\lambda_{\max}(\Gamma)\ge9.38758>120/13\), so cont bound holds on the **entire** spike level; discrete \(\max R=60\), clique-flip, \(\Phi=\Phi(C)\). GW insufficient for discrete (SDP\(\cdot\alpha<60\)). Residue-1 random PMs also clique-flip in samples. Open: prove cont bound / clique-flip forall Max-covers. See `evidence/E1_MAXCOVER_CONTINUOUS_BOUND.md`. **Does not close E(1).**

**Prop 15.38 (2026-07-28):** **Proved** by exhaustive \(\binom{45}{5}\): among two-sided Max-covers of size 5 on Paley \(C_{10}\), undercutters are **exactly** the 144 perfect matchings (\(\Delta=1\), \(\Phi=13\)); every two-sided cover with \(\Delta\ge2\) has \(\Phi\ge15\). Total 17154 two-sided \(k=5\) covers. Parallel cert: 80 workers. See `evidence/E1_N10_TWOSIDED_K5.md`. Supports low-\(\Delta\) undercutter pattern; **does not close E(1).**

**Prop 15.39 (2026-07-28):** On stored Max-cover PMs at \(p=5\), clique-flip pair count \(N_{\mathrm{flip}}\ge24\) (observed \(24,48,120\)). Open: prove \(N_{\mathrm{flip}}\ge1\) forall Max-cover PMs. See `evidence/E1_CLIQUE_FLIP_COUNT.md`. **Does not close E(1).**

**Prop 15.40 (2026-07-28):** **Proved:** edge-minimal undercutters satisfy \(\Phi(C\oplus F)\ge\Phi(C)-2\) (edge lip + minimality). Open step to E(1): upgrade to \(m_n\ge\Phi(C)-2\) for all Seidel \(A\) (no deeper far undercut). If that holds, gap \(O(1)\) \(\Rightarrow\) E(1) \(\Rightarrow L=\tfrac12\). **Does not close E(1) yet.**

**Prop 15.41 (2026-07-28):** **First-hit lemma proved** (first undercutting prefix on any add-edge chain has gap \(\le2\)). **Dangerous-edge criterion proved** (descent below \(\Phi-2\) requires rigid \(\sigma\)-alignment on all maximisers). **No-descent lemma OPEN** in general; if proved on the \(\rho=1\) family then \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Certified n=10:** all 144 PM undercutters have **0** dangerous edges; every add-1 returns \(\Phi\ge15\); multi-edge samples stay \(\ge13\). Evidence: `evidence/E1_NODESCENT.md`, `e1_n10_nodescent.json`. **Does not close E(1).** **F13:** do not claim \(m_n\ge\Phi-2\) from 15.40 alone.

**Prop 15.42 (2026-07-28):** **Max± dichotomy proved:** \(\Phi(A)<\Phi-2\) only possible for deep two-sided covers (\(s_+\ge2\), \(s_-\le-2\)); every \(s_+\le1\) or \(s_-\ge-1\) matrix has \(\Phi\ge\Phi-2\). **Counting freeness proved.** **Tight \(S\equiv1\) and \(S\equiv2\) no-descent proved** (frame mean \(1/p\)). **Type I freeness when \(N_1>N(p+1)/(2p)\) or \(k\le2p-2\) proved.** **Equivalence:** \(m_n\ge\Phi-2\) iff no-descent on all gap-2 undercutters (given 15.40+parity). **n=10 complete** (PM Type I strict freeness; C6 tight \(S\equiv2\)). Open: no-descent for Type I with large \(k\) and small \(N_1\), and deep non-tight \(k>2p\). **Does not close E(1).**

**Prop 15.43 (2026-07-28):** No-descent **proved** for Type I with freeness (strong) and tight deep \(S\equiv2\) (weak). Type I freeness-failure **isolated** to equality cases reducing to tight size-\(2p\) covers; at \(p=3\) 1-bit spike gives \(\Phi\ge\Phi-2\) for all-even-degree tight covers. **n=10 \(m_{10}=\Phi-2\) closed.** Residual for general \(p\): tight \(S\equiv2\Rightarrow\Phi\ge\Phi-2\) beyond \(p=3\); deep non-tight; \(k=3p-2\) boundary. **Does not close E(1) / \(L\).**

**Bi-tight / Prop 15.45 (2026-07-28):** Stars never bi-tight (all \(p>2\)). Integral bi-tight levels \(2,3,4\) **MILP-infeasible at \(p=5\)**; non-star size-\(p\) tight infeasible; deep two-sided \(k=10,12\) infeasible; size-\(2p\) Max+ covers have \(\max S_{-}\ge2\). Avg degree of level-2 bi-tight is \(4p/(p^2+1)<1\) for \(p\ge5\). Master lemma (Prop 15.44): tight Max+ covers either \(\max_{\mathrm{Max}_{-}}S\ge0\) or bi-tight. See `evidence/E1_BITIGHT.md`, `e1_star_bitight_obstruction.json`. Residual: lift to all \(p\ge5\); deep non-tight \(k>2p\).

**n=26 exact sparse MITM (2026-07-27):** shipped `phi_mitm` (exact \(\Phi\) for even \(n\le28\)). Random matchings, cycles \(C_4\)–\(C_{26}\), stars, and random \(k\le20\) flips of Paley \(C_{26}\): **0 undercuts of \(\Phi=65\)** (min observed 67 on single edges). Matching undercut of \(n=10\) does not lift. Consistent with \(k_\star=0\) at \(n=26\), not a proof. See `evidence/E1_N26_SPARSE_EXACT.md`.

**MO thread re-audit (2026-07-27):** full re-read of [MO 413935](https://mathoverflow.net/questions/413935) + comments via SE API — **0 answers; no hidden proof**. Author already knew liminf \(\ge2^{-5/2}\) and that such bounds do **not** settle existence; not interested in numerics. Starr’s Bowlin/multipartite pointer is already dead for existence (§9–§10). **Time savings:** stop BH re-derivation, stop multipartite reopen, stop pure-SA as deliverable; only load-bearing E(1) (or permanent relative gap / non-existence) counts. Details: `evidence/MO_THREAD_REAUDIT.md`.

**Existence of \(L\) remains OPEN** (2026-08-06). Soft-close \(L=\tfrac12\) CLOSED via Props 15.167–171 was **retracted 2026-08-06** (disj Gsum LB hinge; 15.158). Entries below through 15.166 are historical Path-C residual notes; residual/16N remains optional open.
- **Prop 15.105 (Norton/Fickus variance identity; 2026-07-31):** **Proved:** ∑(λ_α(Φ|Z)−μ̄)² = ‖κ_orth‖_F² (Fickus residual = Φ spectral variance); orth=0 ⇒ Φ scalar = 8(n−2)/(n−6); mult(λ_max)≥d + ‖κ‖²≤96n ⇒ λ_max(Φ)≤16 exactly (Fraction: μ̄+8(p²−9)/(p²−5)=16). **Certified** variance+16N+mult=d at p=3,5,7. **OPEN:** ‖κ‖²≤96n (or direct λ_max≤16) for all p≥5. Evidence: `src/e1_gmin_m4_prop15105.py`, `evidence/e1_gmin_m4_prop15105.json`. **L OPEN.**
- **Prop 15.106 (rest-average-8 + kurtosis residual; 2026-07-31):** **Proved:** mult(λ_max)≥d ⇒ (λ_max≤16 ⇔ rest-mean≥8); at λ_max=16 rest-mean=8 (Wick). ‖κ‖²≤96n ⇔ kurtosis κ₄(D)≤3+12/n. **Certified** p=3,5,7. **OPEN:** κ₄≤3+12/n or λ_max≤16 for all p≥5. Evidence: `src/e1_gmin_m4_prop15106.py`, `evidence/e1_gmin_m4_prop15106.json`. **L OPEN.**
- **Prop 15.107 (16N from mult≥d−1 + room_hyp; 2026-07-31):** **Proved:** mult(λ_max)≥d−1 (15.98) + orth≤room_hyp ⇒ λ_max(Φ)≤16 for all primes p≥3 (Fraction majorization (16−μ̄)²≥room_hyp·(m−m₁)/(m m₁)); Gegenbauer α=−6/(d+4), β=3/((d+2)(d+4)); 4-design ED4≤wick_hi. **Certified** orth≤room_hyp + 16N at p=3,5,7. **OPEN residual:** orth≤room_hyp for general p≥5. Evidence: `src/e1_gmin_m4_prop15107.py`, `evidence/e1_gmin_m4_prop15107.json`. **L OPEN.**
- **Prop 15.108 (residual-Gram/Schur dual; algebraic Thm A*; Parseval T_ρ; 2026-07-31):** **Proved:** Thm A gap = 128(p−3)(p+3)(p⁴−12p²−5)/[(p²−5)²(p²+1)²] ≥0 for all primes p≥3 (polynomial); λ_max(Φ)=4N λ₂(PopP); orth=16N²‖R‖_F² (Fickus residual Gram); orth≤room_hyp ⇔ ∑ρ²≤T_ρ(p) with T_ρ=ρ_min²+room_hyp/24 Max+-free closed form; m4 Parseval expansion. **Certified** ∑ρ²≤T_ρ at p=3,5,7. **OPEN residual:** ∑ρ²≤T_ρ for general p≥5. Evidence: `src/e1_gmin_m4_prop15108.py`, `evidence/e1_gmin_m4_prop15108.json`. **L OPEN.**
- **Prop 15.109 (Φ–m4 identity; Aut δ; PF+rank; 2026-07-31):** **Proved:** E[(yᵀBy)²]=6‖B‖²+8⟨m4,κ_B⟩ on Z; ∑κ_B²; PF+rank λ₂(PopP)<d/N p≥5; residual ⇔ c²≤room_hyp/24 **when** dim E_{4p}^{Aut}≤1. **Cert** p=5 equality. **Caveat:** class_key ≠ faithful m₄ partition at p=7 (15.112/F19). **OPEN:** general residual. Evidence: `src/e1_gmin_m4_prop15109.py`. **L OPEN.**
- **Prop 15.110 (closed Max+ IDs; ρ_min<budget p≥7; 2026-07-31):** **Proved:** ∑κ∏=n(n−1)(n−2)/8 and e₄=−p(p−1)(p+1)(p+4)/12 on Max+; ⟨m₄,κ⟩ closed; ρ_min²<room_hyp/24 for all primes p≥7 (polynomial); residual for p≥7 follows from δ²≤ρ_min²; c=Q_0(halfspace). **Cert** δ²≤ρ_min² at p=3,5,7. **OPEN:** δ²≤ρ_min² general p. Evidence: `src/e1_gmin_m4_prop15110.py`. **L OPEN.**
- **Prop 15.111 (pair Schur; closed α_κ,α_ρ; Φ=μ̄+8⟨δ,κ_B⟩; 2026-07-31):** **Proved:** zero-diag ∑κ_C κ_B trace identity; α_κ=(p²+2)/(4p²) on Z; pair=(p²+11)/(4(p²−5)); α_ρ=(7p²+5)/(2p²(p²−5)); channel recon from β_b=6/p, β_T=6(3p²+5)/p²; Φ excess is pure δ-channel so 16N⇔max⟨δ,κ_B⟩≤(n−10)/(n−6). **Cert** Schur scalarity of κ/p²,b,Tb,ρ_min on Z at p=3,5,7. **OPEN:** δ²≤ρ_min² general p. Evidence: `src/e1_gmin_m4_prop15111.py`, `evidence/e1_gmin_m4_prop15111.json`. **L OPEN.**
- **Prop 15.112 (design moments; conf ‖κ‖²; ED4 dictionary; 2026-07-31):** **Proved:** conference ‖κ‖₂²=(n p⁴/8)(n−6)+n(n−1)/2; Max+ antipodal; ED4=ED4_flat+24δ² so δ²≤ρ_min²⇔ED4≤ED4_suf; ED4_suf<ED4_bud for p≥7. **Cert** E[yyᵀ]=2P₊, E[D²]=2n, δ²≤ρ_min² at p=3,5,7. **Note:** class_key not m₄-equitable at p=7 (F19). **OPEN:** ED4≤ED4_suf general p. Evidence: `src/e1_gmin_m4_prop15112.py`. **L OPEN.**
- **Prop 15.113 (⟨f_y,Tκ⟩; ED4 via W; Q_δ criterion; 2026-08-01):** **Proved:** ⟨m₄,Tκ⟩=2p(p⁴−1); ⟨κ,Tκ⟩=0; E[W]=n/2, ED4=3n²+4E[W²]; Q_δ(y)≤ρ_min²⇒δ²≤ρ_min². **Cert** fy constancy + ED4≤ED4_suf at p=3,5,7. **OPEN:** Q_δ≤ρ_min² (or ED4≤ED4_suf) ∀p≥5. Evidence: `src/e1_gmin_m4_prop15113.py`. **L OPEN.**
- **Prop 15.114 (γ-calculus; Tf_y multiplicative; ∑γ,∑γ²; 2026-08-01):** **Proved:** Tf_y=(4p−2γ_y)⊙f_y; δ⊥γ⊙f; ∑γ=(6/p)C(n,4); ∑γ²=6C(n,4)+n(n−1)(n−2)/4 (adjacent edges cancel via Cy=py); ‖Tf_y‖² closed. **Cert** p=3,5. **p=7:** three ED4(y) types (not 2-point homogeneous); all ≤ED4_suf; pair-average ED4≤ED4_suf. **OPEN residual:** pair-average ED4≤ED4_suf ∀p≥5. Evidence: `src/e1_gmin_m4_prop15114.py`. **L OPEN.**
- **Prop 15.115 (Aρ=b; δ=P_{E_4p}m₄; spectral moments; 2026-08-01):** **Proved:** E[γ⊙f]=2κ/p (2-design)⇒Aρ=b for Max+ residual; κ⊥E_{4p} (4p∉{0,±μ})⇒δ=P_{E_4p}m₄=E_y P f_y; m1=4p−12/p, Var_spec=24(p²−3)(p²−4)/(p²(p²−2)). **Cert** resolvent+moments p=3,5; ED4≤suf p=3,5,7. Full E_4p energy of f_y too crude for residual. **OPEN:** coherent mass ‖E P_{E_4p} f_y‖₂²≤ρ_min² ∀p≥5. Evidence: `src/e1_gmin_m4_prop15115.py`. **L OPEN.**
- **Prop 15.116 (e4↔ED4↔δ dictionary; coherent mass; Aut-line; 2026-08-01):** **Proved:** e₄(s)=s⁴/24+((4−3n)/12)s²+n(n−2)/8; ∑m₄²=ED4/24+n(4−3n)/6+n(n−2)/8; ⟨κ,ρ_min⟩=n(n−1)(n−2)(n−6)/(2p²(p²−5)); flat Pythagoras identity; coherent mass=δ²; Aut-line criterion when dim E_4p^{Aut}≤1. Min-dist ED4 envelope **too weak** for p≥5. **Cert** δ²≤ρ_min² p=3,5,7. **OPEN:** coherent mass ≤ρ_min² ∀p≥5. Evidence: `src/e1_gmin_m4_prop15116.py`. **L OPEN.**
- **Prop 15.117 (Path C hyp residual primary; ρ_min pairings; 2026-08-01):** **Proved:** Path C residual ⇔ δ²≤room_hyp/24 with closed form 4(p²−9)(p²−1)²/(3(p²−5)(p²+1)); ρ_min²≷hyp by p (gt at 5, lt for p≥7); slack κ_hyp−κ²=24(hyp−δ²); ⟨b,f_y⟩=⟨ρ_min,b⟩=2(p⁴−1)/p; ⟨ρ_min,m₄⟩ closed; E⟨b,γ⊙f⟩=0; if pointwise then ⟨ρ_min,f⟩=4(p⁴−1)/(3(p²−5)). **Cert** bgf=0 at p=3; hyp residual at p=3 (and p=5,7 when Max+ cache present). **OPEN:** δ²≤room_hyp/24 ∀p≥5. Evidence: `src/e1_gmin_m4_prop15117.py`. **L OPEN.**
- **Prop 15.118 (pointwise bgf=0; T²κ; ρ_min·f closed; 2026-08-01):** **Proved:** ⟨κ,γ⊙f⟩=p(p²+1)(p²−1)(p²−4)/4; ⟨T²κ,m₄⟩=8p²(p⁴−1); pointwise ⟨b,γ⊙f_y⟩=0 on Max+ (via ⟨T²κ,f_y⟩=8p²(p⁴−1)); ⟨ρ_min,f_y⟩=4(p⁴−1)/(3(p²−5)). **Cert** full Max+ p=3,5. **OPEN:** δ²≤room_hyp/24 ∀p≥5 (need m₄/ED4 pin). Evidence: `src/e1_gmin_m4_prop15118.py`. **L OPEN.**
- **Prop 15.119 (residual budget dictionary; weight enum; halfspace pin; 2026-08-01):** **Proved:** ED4_flat/ED4_bud/EW2_flat/EW2_bud/m4f_flat/m4f_bud closed rational forms; residual ⇔ ED4≤ED4_bud ⇔ E[W²]≤EW2_bud (and ⇔ ⟨m₄,f⟩≤m4f_bud when Q_δ constant); Max+ dots ≡2 (mod 4), |D|≤p²−2p−1, crude 2n D_max² too weak. **Cert** weight enum + saturation ED4=ED4_bud at p=3,5; orth·N=147456 at p=5; p=7 ratio δ²/room=124875/669124 (prior). **OPEN:** δ²≤room_hyp/24 ∀p≥5 (need independent E[W²]/⟨m₄,f_hs⟩ UB). Evidence: `src/e1_gmin_m4_prop15119.py`. **L OPEN.**
- **Prop 15.120 (pointwise E[W²] factorization; Pythagoras; majorization; 2026-08-01):** **Proved:** E_z[W_y²]=EW2_flat+6 Q_δ(y) for every Max+ y (disj-edge residual); ∑m₄²=F(p)+δ² with F=‖κ/p²+ρ_min‖₂²=m4f_flat Max+-free; PSD majorization ED4≤2n³ (too weak vs bud for p≥5); CS-γ / LP / maj documented dead. **Cert** EW2 constant + factorization at p=3,5. **OPEN:** δ²≤room_hyp/24 ∀p≥5 (need weight enumerator / Gauss sums). Evidence: `src/e1_gmin_m4_prop15120.py`. **L OPEN.**
- **Prop 15.121 (spectral residual dictionary; Frobenius form; 2026-08-01):** **Proved:** ED4=4n²+4N⁻²‖FFT|_{1^⊥}‖_F²; E[W²]=‖FFT‖_F²/N²; EW2_flat=(n²+T²/m)/4; residual ⇔ ‖FFT|_{1^⊥}‖_F²≤N²(EW2_bud−d²) ⇔ ∑(λ_α−μ̄)²≤room_hyp (Φ variance); 16N is op-norm on same FFT (not Frob); H/16N·Tr majorization too weak for ED4_bud. **Cert** spectral identities p=3,5. **OPEN:** δ²≤room_hyp/24 ∀p≥5. Evidence: `src/e1_gmin_m4_prop15121.py`. **L OPEN.**
- **Prop 15.122 (Max+ disagreement u∈V₊; Aut-line; λ_max(T); 2026-08-01):** **Proved:** for y,z∈Max+, u=(y−z)/2 satisfies uᵀCu=pk and u∈V₊ (ternary V₊ characterisation of dots); λ_max(T)<4p⇒δ=0 (p=3); λ_max=4p at p=5,7; δ∈E_{4p}^{Aut}, residual⇔c²≤room_hyp/24 if dim≤1; FFT budget B(p)=EW2_bud−d² Max+-free; LP/PGL still weak (F18). **Cert** identity+spectra p=3,5. **OPEN:** δ²≤room_hyp/24 ∀p≥5 (ternary weight enumerator / Gauss Q_0). Evidence: `src/e1_gmin_m4_prop15122.py`. **L OPEN.**
- **Prop 15.123 (switching; conference srg regular sets; dual Krawtchouk residual; 2026-08-01):** **Proved:** after switch C'=D_yCD_y, Max+↔V₊∩{0,1}ⁿ via w=(1−z')/2; G=srg(n,p(p−1)/2,μ−1,μ) with μ=((p−1)/2)²; supports are regular sets, k∈{0,n}∪even[p+1,p(p−1)]; B_i=W_i; A'₄=∑m₄²; residual⇔A'₄≤m4f_bud; χ_S−(k/n)1∈V₊∩1^⊥ two-valued. **Cert** srg+W_k p=3,5. **OPEN:** closed W_k / A'₄≤m4f_bud ∀p≥5. Evidence: `src/e1_gmin_m4_prop15123.py`. **L OPEN.**
- **Prop 15.124 (closed weight moments j≤3; E[k⁴] partition; residual as R₄; 2026-08-01):** **Proved:** E[k]=n/2, E[k²]=n(n+2)/4, E[k³]=n²(n+6)/8 Max+-free (E[s²]=2n, E[s³]=0); exact_≤3 partition of E[k⁴]; ED4=ed4_from_exact3+16 R₄; residual⇔R₄≤R₄_bud⇔A'₄≤m4f_bud; Hamming Delsarte max A'₄ saturates only at p=3, weak for p≥5; W_{p+1}=d (Hoffman) cert p=3,5. **Cert** moments+partition+saturation p=3,5. **OPEN:** closed W_k / A'₄≤m4f_bud ∀p≥5. Evidence: `src/e1_gmin_m4_prop15124.py`. **L OPEN.**
- **Prop 15.125 (perfect 2-colorings; 4-design defect; closed R₄ budget; 2026-08-01):** **Proved:** W_k=# of τ-equitable bipartitions of conference srg; Max+ is spherical 2-design in V₊, residual=4-design defect (ED4 vs 3n⁴/(d(d+2))); ed4_from_exact3=−(p²+1)(p⁶+3p⁴−25p²+13); R₄_bud closed; Hamming Delsarte+E[kʲ] j≤3 still weak p≥5; antipodal A'_j=A'_{n−j}. **Cert** p=3,5. **OPEN:** closed W_k / R₄≤R₄_bud ∀p≥5. Evidence: `src/e1_gmin_m4_prop15125.py`. **L OPEN.**
- **Prop 15.126 (geometric Hoffman seed; 1-design; simplex bound; 2026-08-01):** **Proved:** subfield line F_p∪{∞} is Hoffman coclique after hs-switch (cert p=3,5,7); Aut⇒1-design algebra; equal-∩ simplex Gram ⇒ b≤d (eq at p=3). **Cert** W_{p+1}=d at p=3,5 only. **Superseded for W=d general:** see 15.127. Evidence: `src/e1_gmin_m4_prop15126.py`. **L OPEN.**
- **Prop 15.127 (closed W_{p+1}; inversive plane; W=d false; 2026-08-01):** **Proved:** F_p-sublines = miquelian inversive plane S(3,p+1,p²+1); **closed** \(W_{p+1}=((1+\chi_4)/2)d+((1-\chi_4)/2)(3p+1)/2\) with \(\chi_4=(-1)^{(p-1)/2}\) (i.e. \(d\) if \(p\equiv1\bmod4\), \((3p+1)/2\) if \(p\equiv3\bmod4\)). **Cert** full coclique enum p=3,5,7,11 (W=5,13,11,17). **Counterexample:** W_{p+1}=d fails at p=7 (11≠25). Max cocliques = regular sublines at these p. **OPEN:** full W_k / 4-design defect ∀p≥5. Evidence: `src/e1_gmin_m4_prop15127.py`. **L OPEN.**
- **Prop 15.128 (full W_k census p=3,5,7; ED4 at p=7; 2026-08-01):** **Certified** complete weight enumerators p=3,5,7 (p=7 via 2²⁵ Max+ gen, N=11452); W_{p+1} matches 15.127; E[kʲ] j≤3 match 15.124; structural zero W_10=0 at p=7; **ED4(p=7)=12835984/409 < ED4_bud** (strict hyp residual). **OPEN:** closed general W_k / defect bound ∀p≥5. Evidence: `src/e1_gmin_m4_prop15128.py`. **L OPEN.**
- **Prop 15.129 (Jensen coherent-mass; Hoffman r̄; 2026-08-01):** **Proved:** δ²=‖E P f_y‖² ≤ E‖P_{E_4p} f_y‖² (Jensen; eq iff P f constant); residual dictionary census-checked p=3,5,7; r̄=W_{p+1}(p+1)/n; **p=7 Hoffman layer is not a 1-design** (r̄=44/25∉ℤ); p=5 Hoffman pair geometry (30 disjoint pairs = W_12/2). **OPEN:** E‖P f‖²≤room_hyp/24 or closed W_k ∀p≥5. Evidence: `src/e1_gmin_m4_prop15129.py`. **L OPEN.**
- **Prop 15.130 (P m₄=δ; ρ_min-sufficient for p≥7; 2026-08-01):** **Proved:** P m₄=δ; Jensen refined; **room−ρ_min²** closed form >0 for all primes p≥7 so δ²≤ρ_min²⇒residual for p≥7; census δ²≤ρ_min² at p=3,5,7; Aut-line program restated. **OPEN:** δ²≤ρ_min² (or E‖P f‖²≤room / Q_0(hs)) for general p≥5. Evidence: `src/e1_gmin_m4_prop15130.py`. **L OPEN.**
- **Prop 15.131 (pair-avg vs basepoint ED4; p=7 Q_δ spectrum; 2026-08-01):** **Proved:** ED4(y)=flat+24Q_δ(y), ED4_pair=flat+24δ²; when Q non-constant, ED4(hs)≠ED4_pair. **Certified p=7 (W=86):** exactly 3 Q_δ types (cnt 2352/8400/700; Q=−124800/4499, 82176/4499, 200448/4499); **true δ²=19180800/1840091** (prior W-based 82176/4499 was Q_δ(hs), not δ²); max Q≤ρ_min² with slack 812048/220451; Var(Q)>0. Residual holds at p=7 (pair + pointwise). **OPEN:** max Q_δ≤ρ_min² (or δ²≤ρ_min²) for general p≥7. Evidence: `src/e1_gmin_m4_prop15131.py`. **L OPEN.**
- **Prop 15.132 (Max+-free dictionary; Aut δ; γ-parity; dead envelopes; 2026-08-01):** **Proved:** residual ⇔ ∑m₄²≤m4f_bud (RHS Max+-free); δ Aut-invariant ⇒ Q_δ constant on Aut-orbits (≥3 orbits at p=7); γ∈{−6,…,6} even, formal 4p-fiber γ=0. **Certified:** γ=0 mass constant 4350 at p=5; 3-valued at p=7 sample; moment-LP and pole+D_max envelopes **dead** (exceed ED4_suf) for primes p=5..19. **OPEN:** Max+-free max Q_δ/δ² bound for all p≥7. Evidence: `src/e1_gmin_m4_prop15132.py`. **L OPEN.**
- **Prop 15.133 (class_key Bose–Mesner; F19 quantitative; Aut-line; 2026-08-01):** **Proved:** Aut-line form (dim E_{4p}^{Aut}≤1 ⇒ residual ⇔ c²≤room). **Certified class_key T-spectrum:** dim E_{4p}^{ck}=**0,1,0** at p=3,5,7 (λ_max=4p only at p=5). **F19 quantitative:** at p=7 nul_ck=0 but δ²=19180800/1840091>0 ⇒ δ∉V^{ck}. **CR dead:** PGL cross-ratio orbits not Aut(C) (κ not constant). Aut-line works at p=5; fails as carrier at p=7. **OPEN:** true Aut(C) orbits / Gauss Q_0 / max Q_δ. Evidence: `src/e1_gmin_m4_prop15133.py`, `evidence/e1_gmin_m4_prop15133_classkey_spectrum.json`. **L OPEN.**
- **Prop 15.134 (strict Aut(C) Bose–Mesner; residual projection; 2026-08-01):** **Proved:** strict Aut G = affine square-semilinear (∞ fixed), |G|=p²(p²−1); inversion/PGL not ≤Aut(C). **Certified:** G-orbits 9/42/128 at p=3,5,7 (κ-const, T-equitable; finer than class_key at p=5,7); **dim E_{4p}^{G}=0,2,7** (Aut-line dim≤1 **fails** for G at p≥5); residual proj δ=P m₄ recovers δ²=1536/65 (p=5) and 19180800/1840091 (p=7). G **carries** residual at p=7 (unlike class_key). **OPEN:** Gauss m₄ on G-orbits for general p. Evidence: `src/e1_gmin_m4_prop15134.py`. **L OPEN.**
- **Prop 15.135 (coherent-mass spectral form; hs char sums; G·hs dead; 2026-08-01):** **Proved:** residual ⇔ ∑c_j²≤room on ONB of E_{4p}^{G}; f_hs Max+-free F_p indicators with ∑f_hs=e₄; P_G(1)=P_G(κ)=0 so order-3 moments do not pin c_j. **Certified:** G·hs ⊊ Max+ (60/260 at p=5, 168/11452 at p=7); δ² from m₄^{Ghs}≫room (**dead** UB). **OPEN:** full Max+ Gauss sums for c_j / m₄ on G-orbits. Evidence: `src/e1_gmin_m4_prop15135.py`. **L OPEN.**
- **Prop 15.136 (Max+-free flat on G-orbits; m4=flat+δ; 2026-08-01):** **Proved:** flat=κ/p²+A⁺(Tκ/p²) on G-orbits is **Max+-free for general p**; ‖ρ_min‖²=ρ_min², ‖flat‖²=m4f_flat; m₄=flat+δ with δ∈E_{4p}^{G}. **Certified:** geom type (∞,κ,fp_dim) does **not** determine m₄; residual = free c_j only. Partial p=5 ∞-orbit m₄ ∈{−1/5,−21/65}. **OPEN:** Gauss for c_j over full Max+. Evidence: `src/e1_gmin_m4_prop15136.py`. **L OPEN.**
- **Prop 15.137 (c_j over Max+ G-orbits; p=5 two-type formula; 2026-08-01):** **Proved:** Q_j G-equivariant on Max+; c_j=∑_t w_t Q_j(y_t) on H_+ orbits; **p=5 exact** c_j=(3/13)Q_j(hs)+(10/13)Q_j(y_*); Q_j(hs) Max+-free. **Certified census:** p=3 r=1; p=5 r=2 (w=3/13,10/13); p=7 r=5 types in H_+. **OPEN:** Max+-free y_t and Q_j(y_t) for non-hs types, general p. Evidence: `src/e1_gmin_m4_prop15137.py`. **L OPEN.**
- **Prop 15.138 (Max+-free non-hs y_*; p=5 residual Max+-free; 2026-08-01):** **Proved:** y_*=hs⊙z'_S via lex-first norm circle S={u:N(u−t)=c} Hoffman coclique of C'=D_hsCD_hs (Max+-free). **p=5:** (t,c)=(0,3) gives non-hs type; two-type character sum closes residual **without Max+ census** (∑c_j²=1536/65). **Certified:** construction works p=5,7,11; p=13 empty; p=7 covers only 2 of 5 H_+ types. **OPEN:** remaining orbit types for all p≥7. Evidence: `src/e1_gmin_m4_prop15138.py`. **L OPEN.**
- **Prop 15.139 (affine halfspaces + double Seidel–norm-circle; p=7 all size classes; 2026-08-01):** **Certified:** affine halfspaces |S|=(p+1)/2 are Max+ at p=5,7 (all S). **p=7 AP dichotomy:** 4-AP S⇒|O|=84; non-4-AP⇒|O|=56 (21+14). **Double switch** y=y0⊙z_nc(C0) covers 588 (hs⊙nc), 1176 (y56⊙nc), 294 (y_nc⊙nc). **All five H_+ size classes Max+-free at p=7.** **OPEN:** Max+-free weights + Q_j on every orbit (incl. mult-4 of 1176) for residual at p=7 without census; general p≥5. Evidence: `src/e1_gmin_m4_prop15139.py`. **L OPEN.**
- **Prop 15.140 (weights |G|/|Stab|; character-sum Q_j; p=7 residual form; 2026-08-01):** **Certified:** w_t=|O_t|/|H_+| via stab (sizes 56/84/294/588/1176×4, stabs 42/28/8/4/2); |H_+|=5726; c_j=∑w_t Q_j(y_t) recovers **δ²=19180800/1840091≤room=3072/55**. **7/8** orbits Max+-free geometric Q; **one size-1176 OPEN** for Max+-free y. **OPEN:** that type + general p. Evidence: `src/e1_gmin_m4_prop15140.py`. **L OPEN.**
- **Prop 15.141 (size-12 Seidel partner; p=7 residual Max+-free; 2026-08-01):** **Proved/cert:** y_♯=y0⊙z_T with y0 affine S={2,3,4,5}, T size-12 field set — last size-1176 type Max+-free. **All 8 H_+ types Max+-free.** Free c_j recovers δ²≤room — **p=7 residual Max+-free closed.** Bi-tight empty at p=5,7 via residual+mult≥d−1 form. **OPEN:** general p≥5. Evidence: `src/e1_gmin_m4_prop15141.py`. **L OPEN.**
- **Prop 15.142 (uniform affine all-S; k-AP split; fourths/size-12; p=11 sample; 2026-08-01):** **Cert:** all S of size (p+1)/2 give Max+ affine halfspaces at p=5,7,11. k-AP split: p=5 all AP; p=7 21/14→84/56; p=11 multi-way (orbits 132/330/660). **Fourths-coset partners work only at p=5** — not uniform size-12 law. p=7 size-12 fibre has 84 T-sets; explicit T still p=7-only. p=11 Max+-free samples (ystar 3630, dbl 7260…). **OPEN residual p≥11 / general p.** Evidence: `src/e1_gmin_m4_prop15142.py`. **L OPEN.**
- **Prop 15.143 (p=11 affine 6-orbit census; double-switch LB; 2026-08-01):** **Cert:** all 462 affine S form **exactly 6** H_+ G-orbits (132×1, 330×2, 660×3) with constructive samples. Double-switch: ystar 3630; nc density class-dependent (132-rich, 660-empty); |H_+|≥28182 LB (sizes 132/330/660/3630/7260). **Affine types complete at p=11; full H_+ type list OPEN.** Residual p≥11 OPEN. Evidence: `src/e1_gmin_m4_prop15143.py`. **L OPEN.**
- **Prop 15.144 (p=11 free orbits; type-enum DEAD; residual redirect; 2026-08-01):** **Cert:** free H_+ orbit |O|=|G|=14520 via affine-132 dbl path (33,3)→(69,9); size-2420 via ystar chain (22,4)→(91,4)→(25,4)→(95,5). Deep dbl produces many free orbits — **type-list residual DEAD for p≥11**. Redirect: type-free Max+-free δ²≤ρ_min² / pointwise Q. Residual closed only p=5,7. Evidence: `src/e1_gmin_m4_prop15144.py`. **L OPEN.**
- **Prop 15.145 (type-free residual package; 2026-08-01):** **Proved Fraction:** δ²≤ρ_min² ⇔ ‖ρ‖²≤2ρ_min² ⇔ ‖m₄‖²≤m4f_suf ⇔ ED4≤ED4_suf with closed m4f_suf=(p²−1)(p²+1)(3p⁴+37p²+60)/(24p²(p²−5)), ED4_suf=4(p²+1)(3p⁶−3p⁴+7p²−15)/(p²(p²−5)); ρ_min²/room→5/8; type-free targets (E‖P f‖² / Q_δ / weight enumerator). **Cert** δ²≤ρ_min² at p=5,7 only. **OPEN:** type-free δ²≤ρ_min² ∀p≥7. Evidence: `src/e1_gmin_m4_prop15145.py`. **L OPEN.**
- **Prop 15.146 (type-free R₄/μ₄ residual channel; 2026-08-01):** **Proved Fraction:** ED4_from_exact3=−n⁴+28n²−40n; δ²=(2/3)(R₄−R4_flat); δ²≤ρ_min²⇔R₄≤R4_suf=R4_flat+(3/2)ρ_min² ⇔μ₄≤μ4_suf⇔Ē[∏₄w]≤R4_suf/(n)₄; spectral mass w*=(p²−3)(p²−4)/(p²(p²−1)) gives ‖P f‖² UB too weak for Jensen residual. **Cert** R₄ channel p=5,7. **OPEN:** R₄≤R4_suf ∀p≥7. Evidence: `src/e1_gmin_m4_prop15146.py`. **L OPEN.**
- **Prop 15.147 (inclusion-density residual; ULC near-miss; 2026-08-01):** **Proved:** R₄=E[k^{underline 4}]; d₁=1/2, d₂=(p²+1)/(4p²), d₃=(p²+3)/(8p²) Max+-free; residual⇔d₄≤d4_suf; U=d₃²/d₂ < d4_suf all primes p≥5 (P(p²)=(p²−1)(p⁸−…)>0) so ULC would close residual; **ULC fails** at p=5,7 (d₄/U≈1.036,1.019) while residual holds. **OPEN:** d₄≤d4_suf ∀p≥7. Evidence: `src/e1_gmin_m4_prop15148.py`. **L OPEN.**
- **Prop 15.149 (size-bias residual; independence excesses; 2026-08-01):** **Proved:** residual⇔E_μ[k]≤k_suf for μ∝k^{underline 3}dW (triple-covered regular sets); k_suf=3+8 R4_suf/(n(n−2)(n+2)); k_flat shift→3⁻ closed; d₂−1/4=1/(4(n−1)), d₃−1/8=3/(8(n−1)), d4_flat>1/16; Gauss program via λ_τ on srg triples. **OPEN:** E_μ[k]≤k_suf ∀p≥7. Evidence: `src/e1_gmin_m4_prop15149.py`. **L OPEN.**
- **Prop 15.150 (srg triples; λ_e; Max+-free π_e; 2026-08-01):** **Proved:** n_e closed; under Aut+affine, λ_e=N(p+3−2e)/(8p) (cert p=5,7); π_e=n_e(p+3−2e)/Tot Max+-free; residual⇔∑π_e m_e≤k_suf with m_e=mean regular-set size through type-e triple. **OPEN:** bound m_e. Evidence: `src/e1_gmin_m4_prop15150.py`. **L OPEN.**
- **Prop 15.151 (m_e covariance formula; 2026-08-01):** **Proved:** regular-set t-identities; E[t_e]=n_e(p+3−2e)/(8p) Max+-free; m_e=n/2+8p Cov(k,t_e)/(n_e(p+3−2e)); residual⇔Cov(k,C(k,3)) budget; exact m_e at p=5; CS/Popoviciu too weak. **Note:** Thm B constancy-on-weight overstated for p=7 — corrected in 15.152. **OPEN:** bound Cov/m_e. Evidence: `src/e1_gmin_m4_prop15151.py`. **L OPEN.**
- **Prop 15.152 (free-param t₃; multi-orbit; residual≡R₄; 2026-08-01):** **Proved:** t-vector = affine(k)+γ t₃ (one free param); p=5 closed t₃(α)=0 (α≤3), 3α²−21α+40 (α≥4); **p=7 multi-type** at k∈{16…34} (3–8 orbits/k; pure t_e(k) **dead**); Cov(k,C(k,3))=(R₄+(3−n/2)E₃)/6 with R₄=E[k^{underline 4}] ⇒ residual ∑π m ≡ R₄ channel; Cov(k,t_e) reduces to Cov(k,t₃)+weight moments (γ=(−1,3,−3,1)); E[t₃]=n₃(p−3)/(8p). **Cert** W=86 full Max+ t-census p=5,7. **OPEN:** R₄≤R4_suf / char-sum m_e / Cov(k,t₃). Evidence: `src/e1_gmin_m4_prop15152.py`. **L OPEN.**
- **Prop 15.153 (switched μ₄ residual dictionary; 2026-08-01):** **Proved:** switched e₁=e₃=0, e₂=1/(n−1); d₄=(1+6/(n−1)+μ₄)/16; residual ⇔ μ₄≤μ4_suf with **μ4_flat=(3p²+17)/(p²(p²−2)(p²−5))**, **μ4_suf=(3p⁴+37p²+60)/(p⁴(p²−2)(p²−5))**; exact m_e at p=7. **Cert** μ₄≤μ4_suf p=5,7. **OPEN:** Paley/Weil bound μ₄≤μ4_suf. Evidence: `src/e1_gmin_m4_prop15153.py`. **L OPEN.**
- **Prop 15.154 (avg(χκ)=3/(n−3); η residual; 2026-08-01):** **Proved Max+-free:** avg(χ κ)=3/(n−3)=3/(p²−2) by conference+halfspace C₂ row-sum algebra; μ₄=κ_main+η with **κ_main=3/(p²(p²−2))**; residual ⇔ **η≤η_suf** with η_suf=4(13p²+15)/(p⁴(p²−2)(p²−5)); κ_main < μ4_flat < μ4_suf. **Cert** η≤η_suf p=5,7. **OPEN:** Weil/Aut bound on η=avg(χ Ext)/(4p). Evidence: `src/e1_gmin_m4_prop15154.py`. **L OPEN.**
- **Prop 15.155 (Aut-line e₄/Tχ/Q; η=c₁R₄+c₀; 2026-08-01):** **Proved:** e₄(s)=s⁴/24+((−3n+4)/12)s²+n(n−2)/8; **Tχ=χ(4p−2σ_z)**; **Q=(p/4)[s²(n−4)+n(6−n)]** on C₂-evecs; avg(χ Ext)=E[(4p e₄−2Q)/C(n,4)]; **η=c₁ R₄+c₀** with c₁=16/(n)_4, c₀=−(p⁴+4p²−9)/(p²(p²−2)). Residual = pure E[s⁴] of switched Max+. Crude E[s⁴]≤2n³ **dead**. **OPEN:** Weil/Paley E[⟨z,y⟩⁴] or 3-design defect. Evidence: `src/e1_gmin_m4_prop15155.py`. **L OPEN.**
- **Prop 15.156 (κ₄ residual dictionary; 2026-08-01):** **Proved:** residual ⇔ **κ₄≤κ4_suf** with κ₄=E[s⁴]−12n²; **κ4_flat=16(p²+1)(p²+3)/(p²−5)**, **κ4_suf=4(p²+1)(9p⁴+22p²−15)/(p²(p²−5))**; bridge **κ₄=(n)_4 η−16n**; spherical 4-design is LB only (wrong direction); crude 2n³ + moment-LP on allowed k **dead**. **Cert** κ₄≤κ4_suf p=5,7 (ratios 0.90, 0.66). **OPEN:** Weil/3-design-UB on κ₄. Evidence: `src/e1_gmin_m4_prop15156.py`. **L OPEN.**
- **Prop 15.157 (Gegenbauer defect μ_G4; 2026-08-01):** **Proved:** on spherical 2-designs E[s⁴]=n⁴(a₀+a₄ μ_G4) with **a₀=3/(d(d+2))**, **a₄=(d²−1)/((d+2)(d+4))**, μ_G4=E[Q₄(s/n)]≥0; residual ⇔ **μ_G4≤μ_G4_suf** with closed **μ_G4_suf=4(21x³+19x²+35x−75)(x+9)/[x(x−5)(x+1)³(x+3)(x−1)]** (x=p²). **Cert** p=5,7 (ratios 0.94, 0.83); not a 3-design. Dead: μ≤1, μ≤1/h₄ (false p=5), μ≤d/h₄ (too weak p≥7). **OPEN:** Weil/Aut bound on μ_G4. Evidence: `src/e1_gmin_m4_prop15157.py`. **L OPEN.**
- **Prop 15.158 (Q₄ closed; Max+ non-scheme; 2026-08-01):** **Proved:** **Q₄(t)=[(d+2)(d+4)t⁴−6(d+2)t²+3]/(d²−1)**, Q₄(0)=3/(d²−1); Max+ **1-homogeneous** with Gram spectrum **2N (×d)+0 (×N−d)**; **not** an IP-association scheme (p=5); pole decomposition μ_G4≤2/N+P(E)+Q₄(0)P(Eq); conditional residual if N≥N_*=⌈2/μ_G4_suf⌉ and nonpositive-Q₄ support — N_* holds p=5,7 but support fails; Chebyshev split **dead**. **OPEN:** Weil/Aut-coherent config/closed W_k. Evidence: `src/e1_gmin_m4_prop15158.py`. **L OPEN.**
- **Prop 15.148 (relaxed-ULC residual calculus; 2026-08-01):** **Proved:** C_act=p²(d₄/U−1)=C_flat+κδ²; C_max=Q(p²)/((p²−5)(p²−2)(p²+3)²); residual⇔C_act≤C_max; for p≥7, C_max≥C₇=79923/87373 and d₄≤U(1+C₇/p²)⇒residual (uniform target). C_max→1⁻. Census C_act≈0.90–0.91. **OPEN:** prove relaxed ULC. Evidence: `src/e1_gmin_m4_prop15148.py`. **L OPEN.**
- **Prop 15.159 (Φ|Z spectrum + dual gap; 2026-08-03):** **Proved/cert:** exact Φ spectra p=5,7 (mult top=d); design thr 16(d−2)/d algebra; dual gap G=(d/32)(16I−Φ) ≽I at p=5,7; 16N chain mult≥d+‖κ‖²≤96n. **OPEN:** G≽I / residual general p. Evidence: `src/e1_gmin_m4_prop15159.py`. **L OPEN.**
- **Prop 15.160 (H vs dual-gap algebra; 2026-08-03):** **Proved:** H−thr_ray=(p−5)(p+1)/(2d) so H≤thr_ray for p≥5; H⇒G≽I⇒16N; H⇒16N (15.63). Cert ray=H at p=5. **OPEN:** Hypothesis H (ray_max≤H(p)) all p≥5. Evidence: `src/e1_gmin_m4_prop15160.py`. **L OPEN.**
- **Prop 15.161 (Φ-frame; mult≥d+κ₄≤48n; 2026-08-05):** **Proved:** ‖v_y‖²=n(n−2); ⟨v_y,v_z⟩=(y·z)²−2n; 16N budgets under mult≥d. Cert p=5,7. **OPEN:** mult≥d and κ₄≤48n general p. Evidence: `src/e1_gmin_m4_prop15161.py`. **L OPEN.**
- **Prop 15.162 (maximizers in Z; mult≥d−1; 2026-08-05):** **Proved for all p≥5:** maximizers of Γ lie in Z; **mult(Φ)≥d−1** (PSL min irrep); E[s⁴]=C₀+R type expansion; 16N⇔mult≥d+m₄-mass. **OPEN:** mult≥d; m₄-mass/κ₄. Evidence: `src/e1_gmin_m4_prop15162.py`. **L OPEN.**
- **Prop 15.163 (Wick m₄ mass; Aut₀; H_C; 2026-08-05):** **Proved:** Wick ∑(m₄^W)² closed; η-room n(19p²−3)/(6p²); Aut₀ V₊≅1⊕(d−1) (cert p=5); H_C split 1+(d−1) (cert p=5,7). **OPEN:** η-bound / mult≥d general p. Evidence: `src/e1_gmin_m4_prop15163.py`. **L OPEN.**
- **Prop 15.164 (16N via mult≥d−1+Es4_*; 2026-08-05):** **Proved:** two-level majorization; Es4_* budget under mult≥d−1+λ_min≥6 (no mult≥d needed); κ4_*/η_* equivalents. Cert 16N p=5,7. **OPEN residual:** E[s⁴]≤Es4_*(p) all p≥5. Evidence: `src/e1_gmin_m4_prop15164.py`. **L OPEN.**
- **Prop 15.165 (exact Es4; closed Es4_*/η_*; GoG↔Φ; 2026-08-05):** **Proved:** E[s]=0, E[s²]=2n; GoG spectrum ↔ Φ; closed Es4_*=4(3p⁸−…)/(p⁴−8p²−1), η_*; m₄ is C-eigen. Exact Es4 census p=3,5,7 (spectrum not W_CENSUS). **OPEN:** Es4_* general p. Evidence: `src/e1_gmin_m4_prop15165.py`. **L OPEN.**
- **Prop 15.166 (16N⇔Q₂ thr; Wick C-eigen; 2026-08-05):** **Proved:** Max+ unit 2-design; Wick m₄ C-eigen (does not pin η); λ_max(Φ)=4d(d−1)/N·λ_max(Q₂) ⇒ 16N⇔λ_max(Q₂)≤4N/(d(d−1)). Cert p=5,7. residual_closed_general=false. Evidence: `src/e1_gmin_m4_prop15166.py`. **L OPEN.**
- **Prop 15.167 (bi-tight majorization; 2026-08-05):** **Proved for all primes p≥5:** mult(λ_max Φ)≥d−1 + λ_min≥6 ⇒ L_*=(p⁴+24p²−1)/(2(p²−1)) < 2d ⇒ λ_cycle < d ⇒ **bi-tight empty**. No residual/16N. Evidence: `src/e1_gmin_m4_prop15167.py`.
- **Prop 15.169 (Type I k=3p−2 reduction + deep multi-s; 2026-08-05):** **Proved:** freeness-fail k=3p−2 Fraction structure; Φ 2-Lipschitz edge flip; gap-2 Type I forces s_−=−1; multi-s auto-freeness k≤p(s+1)−2. **OPEN:** s_−≤−1 impossible for freeness-fail Type I (cert p=5 MILP only); deep freeness-fail k≥3p ND. E1/L OPEN. Evidence: `src/e1_gmin_m4_prop15169.py`.
- **Prop 15.168 (E(1) structure; 2026-08-05):** Tight obstruction from 15.167; Type I freeness ND (prior); Type I fail k=2p−1 ND; deep auto-freeness k≤3p−2; deep fail-eq k=3p−1 ⇒ tight L3 empty. **OPEN:** Type I k=3p−2 boundary; deep freeness-fail k≥3p. **E1/L OPEN** (no soft-close). Evidence: `src/e1_gmin_m4_prop15168.py`.

## 1. Exact quantity (do not restate incorrectly)

For \(n\ge2\),

\[
m_n
=
\min_{a_{ij}=\pm1}
\max_{x\in\{\pm1\}^n}
\Biggl|
\sum_{1\le i<j\le n}a_{ij}\,x_i x_j
\Biggr|,
\qquad
\alpha_n=\frac{m_n}{n^{3/2}},
\qquad
L\stackrel{?}{=}\lim_{n\to\infty}\alpha_n.
\]

**Existence:** \(L=1/2\) (sandwich + Paley denseness + E(1) via 15.272). Historical per-prop “L OPEN” lines below are leftover from the search and are not the current claim.

### Equivalent matrix form (`solution.md` §1)

Associate the symmetric zero-diagonal Seidel matrix \(A\) with \(A_{ij}=a_{ij}\). Then

\[
\sum_{i<j}a_{ij}x_ix_j=\tfrac12 x^\top Ax,
\quad
\Phi(A)=\max_x\bigl|\tfrac12 x^\top Ax\bigr|,
\quad
m_n=\min_A\Phi(A).
\]

Same \(\alpha_n\). Tests: `test_form_Q_matches_half_xAx`, `test_equivalence_m_vs_half_phi_on_optimum_n6`.

### Factorization (`solution.md` Prop 15.6)

\[
\Phi(A)=\tfrac12\,n\,\|A\|_{\mathrm{op}}\,\rho(A),
\quad
\rho(A)=\frac{\max_{x\in\{\pm1\}^n}|x^\top Ax|}{n\,\|A\|_{\mathrm{op}}},
\quad
r(A)=\rho(A)\cdot\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}.
\]

Conference matrices uniquely minimize \(\|A\|_{\mathrm{op}}\) (floor \(\sqrt{n-1}\)). Asymptotic optimality of \(m_n\) along conference orders \(\Leftrightarrow\) \(\min_A r(A)=\rho(C)+o(1)\).

---

## 2. Proved results (with pointers)

| Item | Statement | Where |
|------|-----------|--------|
| **Sandwich** | \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) | `solution.md` Main Theorem, Prop 5.2, Thm upper via Paley/conference |
| **Dual-Gaussian LB** | For every Seidel \(A\): \(\Phi(A)\ge n\sqrt{n-1}/\pi\); hence \(m_n\ge n\sqrt{n-1}/\pi\) | `solution.md` Prop 5.2; `dual_gaussian_lower_bound` in `src/minmax_quadratic.py` |
| **BH floor (weaker)** | \(\liminf\alpha_n\ge2^{-5/2}\) | `solution.md` Prop 5.1 |
| **Cut-code identity** | \(m_n=\binom n2-2\rho(D_n)\) with \(D_n=\{\pm(x_ix_j)_{i<j}\}\) | `solution.md` Prop 1.2 / cut-code section |
| **Monotonicity / steps** | \(m_n\le m_{n+1}\le m_n+n\); consecutive \(\alpha\) gaps \(O(n^{-1/2})\) | `solution.md` §3 |
| **Denseness / Paley reduction** | Paley orders \(n_k=q_k+1\) (\(q\equiv1\bmod4\)) satisfy \(n_{k+1}/n_k\to1\); existence of \(\lim\alpha_n\) \(\Leftrightarrow\) convergence along Paley orders alone | `solution.md` Prop 6.1–6.2 |
| **Conference spectral UB** | If \(C^2=(n-1)I\), then \(\Phi(C)\le\tfrac12 n\sqrt{n-1}\) | `solution.md` §2; `spherical_half_bound` |
| **Limsup via Paley \(\rho\)** | \(\limsup\alpha_n\le\tfrac12\limsup_k\rho(C_k)\le1/2\) | `solution.md` Prop 15.8 |
| **Unique min-op / tr\(A^4\) / \(\mathbb E[Q^4]\)** | Conference uniquely minimizes op-norm, \(\mathrm{tr}(A^4)\), and \(\mathbb E[Q^4]\); \(L^2\) mass of \(Q\) universal | `solution.md` Props 15.5, 15.10–15.13 |
| **Exact optimality \(n=6\)** | \(m_6=\Phi(C)=5\) | `solution.md` Cor 15.15; `test_exact_optimality_n6_via_Q4_gap` |
| **\(\rho=1\) for \(n=p^2+1\)** | Paley over \(\mathbb F_{p^2}\): halfspace boolean evec \(Cx=px\), \(\rho(C)=1\), \(\Phi(C)=\tfrac12 n\sqrt{n-1}\) | `evidence/PROOF_rho_eq_1.md`; `solution.md` Main Theorem corollary; shipped `paley_conference_prime_power` + `halfspace_boolean_vector` |
| **Limsup \(\rho=1\) along dense family** | On \(n_k=p_k^2+1\), \(\rho=1\) for all \(k\), \(n_{k+1}/n_k\to1\) | Corollary of \(\rho=1\) theorem |

Soft multipartite / Hadamard / annealed / rank-one blow-up **cannot** force \(\liminf=\limsup\) (`solution.md` §9–§10). Dead: global Q4 path for large \(n\) (Props 15.16–15.19).

---

## 3. Open blockers for settling existence

| ID | Blocker | Notes |
|----|---------|--------|
| **E(1)** | Asymptotic optimality: \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) along Paley (or \(\min r=\rho(C)+o(1)\)) | **Exact** optimality fails at \(n=10\): \(m_{10}=13<\Phi_{\mathrm{Paley}}=15\). Must be asymptotic. Local edge-opt of Paley at small orders is not a proof. |
| **E(2)** | \(\rho(C_k)\to\rho_*\) for all Paley (or full conference) orders | On \(n=p^2+1\), \(\rho\equiv1\) already. For general Paley, \(\rho\) is increasing on small exact orders; interval constructions give constructive \(\rho_{\mathrm{int}}\gtrsim0.99\) at large \(n\) (heuristic/constructive, not full \(\rho\to1\)). |
| **Non-existence path** | Two subsequences with unequal \(\lim\alpha\) | Denseness forces any oscillation to appear along Paley orders too. No certified construction. |
| **Thm E (conditional)** | E(1)+E(2) along general Paley \(\Rightarrow L=\rho_*/2\) | Open. **Shortcut:** E(1) on the dense \(\rho=1\) family alone \(\Rightarrow L=\tfrac12\) (Prop 6.1–6.2: \(\alpha_{n_k}\to\tfrac12\Rightarrow\lim\alpha_n=\tfrac12\)). |
| **Thm F (Stolz)** | If \(\delta_n/\sqrt n\to\ell\) then \(\alpha_n\to\tfrac23\ell\) | Extension-cost regularity open. |

**Do not claim existence from sandwich + denseness alone.**

---

## 4. Numerics inventory

### 4.1 Exact \(m_n\) (load-bearing)

| \(n\) | \(m_n\) | \(\alpha_n\) | Certification |
|------:|--------:|-------------:|:--------------|
| 2 | 1 | 0.354 | live `exact_m` |
| 3 | 3 | 0.577 | live `exact_m` |
| 4 | 4 | 0.500 | live `exact_m` |
| 5 | 4 | 0.358 | live `exact_m` |
| 6 | 5 | 0.340 | live `exact_m`; = Paley \(\Phi\) |
| 7 | 9 | 0.486 | live `exact_m` |
| 8 | 10 | 0.442 | live `exact_m` |
| 9 | 12 | 0.444 | recorded multi-worker Gray (constants in tests + `evidence/exact_m_table.json`) |
| 10 | **13** | 0.411 | recorded multi-worker Gray; **\(m_{10}<\Phi_{\mathrm{Paley}}=15\)** |
| 11 | \(\le17\) (claimed \(=17\)) | \(\le0.466\) | UB witness in prior work; external cut-code package claims exact 17 (Esmaeili–Zaghian counterexample). Treat exact equality as **external claim** unless re-verified. |

**Evidence paths:**

- Live \(n=2..8\): `evidence/exact_m_table.json` (also under `{SCRATCH}/evidence/exact_m_table.json`), regenerated by `regen_exact_table.py` / pytest session fixture.
- Recorded \(m_9,m_{10}\): `tests/test_minmax.py` (`M9_RECORDED`, `M10_RECORDED`); `evidence/exact_m_table.json`.
- Dual-Gaussian vs exact: pytest `test_dual_gaussian_lower_bound_holds_for_exact_m`.

### 4.2 \(\rho=1\) family (certified)

| \(p\) | \(n=p^2+1\) | \(\Phi=\frac12 np\) | \(\rho\) | Check |
|------:|------------:|--------------------:|---------:|:------|
| 3 | 10 | 15 | 1 | halfspace evec |
| 5 | 26 | 65 | 1 | halfspace evec |
| 7 | 50 | 175 | 1 | halfspace evec |

**Evidence:** `evidence/rho1_verify.json`, `evidence/PROOF_rho_eq_1.md`, pytest `test_rho_eq_1_paley_prime_power`, `{SCRATCH}/evidence/rho1_verify.log`.

### 4.3 Heuristic / non-certified (do not treat as exact \(m_n\))

| Claim | Status | Trap |
|-------|--------|------|
| SA “best_UB” from `phi_local` | **Lower** bound on \(\Phi(A)\) only | Minimizing `phi_local` does **not** upper-bound \(m_n\) unless finalists are re-scored with **exact** \(\Phi\) |
| Interval \(\rho_{\mathrm{int}}\ge0.99\) at large \(n\) | Constructive lower bound on \(\rho\) for specific sign patterns / intervals | Supports limsup \(\rho\) near 1; not full E(2); not exact \(\Phi\) |
| Local-search \(\rho\) “dips” | Uncertified | Not non-existence |
| Prior SA at \(n=26\): no certified beater of Paley \(\Phi=65\) | Suggestive local min | Hamming distance of true optimizers may be large (at \(n=10\), 1–2 edge flips from Paley do **not** reach \(m_{10}=13\)) |
| **n=10 structure (certified)** | Hamming-5 threshold; 144 perfect-matching optima; 1-edge local opt | `evidence/N10_STRUCTURE.md`, `n10_structure.json`, `n10_matching_optima.json`; tests `test_n10_*` |
| **n=10 classification N10-C** | Maximizer-drop \(\Leftrightarrow\) \(\Phi=13\); single \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) orbit of size 144 | `evidence/N10_MATCHING_CLASSIFY.md`, `n10_matching_classify.json`; `src/n10_matching_classify.py`; `test_n10_matching_classify_*` |
| n=26 matching probe | 86 random perfect-matching flips of Paley: \(\Phi\ge73>65\) | Matchings do **not** undercut at \(n=26\) (`n26_matching_probe.json`) |
| n=26 SA + exact rescore | 86×5k SA finalists, exact \(\Phi\): none \(<65\); certified \(m_{26}\le65\) | `n26_sa_exact_rescore.json` — `phi_local` alone is **not** a UB |
### 4.4 Independent external artifacts (context only)

- Sol / Codex writeup: sandwich \(1/\pi\)…\(1/2\) + cut-code ([Robby955/mo-413935-ai-attempt](https://github.com/Robby955/mo-413935-ai-attempt)) — aligns with our Prop 5.2.
- Curtis cut-code \(M_{11}=17\): [antipodal-cut-code-k11](https://github.com/CurtisAccelerate/antipodal-cut-code-k11) — finite; does not settle the limit.
- X thread: mostly AI hallucinations (Parisi/SK, graphon uniqueness, Wick holomorphy). See prior triage notes if present; ignore as proofs.

---

## 5. Resume playbook

### 5.1 Next attacks (ranked)

**Stance (2026-08-03):** residual is **localized**; do **not** prove δ²≤room by inventing Prop 15.159. Ask what **structure forces** residual (or kills free δ). Full map: `evidence/STRATEGY_REFRAME_2026-08-03.md`.

1. **Structure ID (preferred over norm thrash)** — Identify residual operator with one of: Bose–Mesner / coherent config (not Max+ IP-scheme — **blocked 15.158**; try **G-orbits / edge BM**, 15.134), **SDP dual** for 16N/orth, **G-irrep / character** of E_{4p}^G, **Weil m₄ on G-orbits**, or (speculative) signed flag algebra. Equality at p=5 is a clue.  
2. **E(1) on \(\rho=1\) family** — \(m_n=\Phi-o(n^{3/2})\) or permanent gap. First-hit + no-descent framework (15.40–15.43); bi-tight empty if \(g_{\min}>T(p)\) (15.47, cert p=5,7). Need general \(g_{\min}\) or deep non-tight ND / \(k_\star\).  
3. **Limit objects of all Φ-minimisers** — Conference as *consequence*, not assumption. Graphon USC already **fails**; need non-hallucinated rigidity.  
4. **E(2) analytic** — \(\rho\to1\) for large Paley (not required for lim=1/2 if E(1) on ρ=1 family).  
5. **Non-existence** — Only via two dense subsequences with **proved** unequal α; denseness (Prop 6.2) mandatory.
### 5.2 Traps to avoid

| Trap | Why |
|------|-----|
| Soft multipartite / Hadamard “\(c_k\to0\)” existence | Cross-block error is \(\Theta(N^{3/2})\); abstract sequences can oscillate inside the sandwich |
| Q4 / fourth-moment shell for large \(n\) | Shell vacuous when \(\Delta_*/3>\max\delta\) (`solution.md` Props 15.16–15.19) |
| Treating `phi_local` min as \(m_n\) UB | Local max underestimates \(\Phi\) |
| Uncertified \(\rho\) dips as non-existence | Need certified \(\Phi\), not SA LB |
| Exact conference optimality for all \(n\) | **False** at \(n=10\) |
| Claiming \(L=1/2\) from \(\rho=1\) alone | Needs **E(1)** on that family (Prop 6.1–6.2 then force full limit); dual-Gauss alone only gives liminf \(\ge1/\pi\) |
| ProcessPool via stdin / single-core pegging | Use script files; `W=ncpu-2`, `OMP_NUM_THREADS=1` per worker |
| Parisi / graphon / Wick “existence proofs” from X | Hallucinations |

### 5.3 Compute policy

- Workers: `full_workers = nproc - 2` when machine idle.  
- BLAS: `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` per process.  
- Scripts as files (never `python -c` + ProcessPool from stdin).  
- Scratch for ephemeral runs: goal `{SCRATCH}` only; ship durable evidence under `evidence/` and tests.

### 5.4 Key code entry points

| Function | File | Role |
|----------|------|------|
| `exact_m(n)` | `src/minmax_quadratic.py` | Exact \(m_n\) for \(n\le9\) (n=9 heavy) |
| `paley_conference_matrix(q)` | same | Prime \(q\equiv1\bmod4\) |
| `paley_conference_prime_power(p)` | same | Order \(p^2+1\) |
| `halfspace_boolean_vector(p)` | same | \(\rho=1\) evec |
| `dual_gaussian_lower_bound(n)` | same | Prop 5.2 floor |
| `phi` / `phi_local` | same | Exact / local \(\Phi\) |

### 5.5 Files map

| Path | Role |
|------|------|
| `HANDOFF.md` | **This document** — resume entry point |
| `solution.md` | Full writeup + obstruction analysis |
| `README.md` | Short status |
| `src/minmax_quadratic.py` | Shipped library |
| `tests/test_minmax.py` | Load-bearing tests |
| `evidence/` | Durable numerics + \(\rho=1\) proof note |
| `evidence/PROOF_rho_eq_1.md` | Full \(\rho=1\) proof |
| `evidence/N10_STRUCTURE.md` | Theorem N10-S: matching undercutters of Paley-\(10\) |
| `evidence/N10_MATCHING_CLASSIFY.md` | Theorem N10-C: maximizer criterion + \(\mathrm{P}\Gamma\mathrm{L}\) orbit |
| `src/n10_structure.py` | Maximizer balance + \(k\)-flip + SA structure campaign |
| `src/n10_matching_optima.py` | Perfect-matching census (144/945) |
| `src/n10_matching_classify.py` | Classification: criterion + \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) orbit |
| `src/n26_matching_probe.py` / `n26_sa_exact_rescore.py` | \(\rho=1\) family probes at \(n=26\) |
---

## 6. Conditional landscape (if E(1)+E(2) close)

- Along \(\rho=1\) family: \(\Phi(C_n)/n^{3/2}\to1/2\).  
- If E(1): \(\alpha_{n_k}\to\tfrac12\) along that dense family \(\Rightarrow\lim\alpha_n=\tfrac12\) by Prop 6.1–6.2 (both liminf and limsup).  
- If E(2) gives \(\rho_*<1\) and E(1) on general Paley, conditional Thm E says \(L=\rho_*/2\).

None of this is established. **Leave Main Theorem as sandwich + OPEN.**

---

## 7. Verification commands

```bash
cd /home/nick/quadratic-minmax-limit
# FULL SUITE — never single-core (F17). Default: pytest.ini -n 86 / scripts/pytest_full.sh
./scripts/pytest_full.sh
# targeted single test OK without -n; full suite must use W=nproc-2
python3 -c "from src.minmax_quadratic import dual_gaussian_lower_bound; print(dual_gaussian_lower_bound(10))"
# expect ~ 9.5493 = 10*3/pi
```

Expected: all tests pass; existence not claimed settled in `HANDOFF.md` or `solution.md` Main Theorem.
