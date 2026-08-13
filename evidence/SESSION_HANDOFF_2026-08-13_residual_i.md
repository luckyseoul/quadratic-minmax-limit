# Session handoff — residual (i) 2026-08-13 (night)

## Verification stamp (2026-08-13 post-15.268)

| Check | Result |
|-------|--------|
| `residual_ii_a/b/full_closed` | **True** / **True** / **True** |
| `gsum` / `residual_i` / `type_I` / `e1` | all **False** |
| `nu_zero_on_kappa1_proved_general` | **True** (15.268) |
| `residual_i_closed_via_268` | **False** |
| envelope / reflection hyp | **False** (criteria proved; hyp open) |
| pytest 15.170/176/217/218/266–268 | exit 0 |
| soft-close | residual (i) / E1 / L honestly **OPEN** |

**Campaign NOT met.** After `/goal clear` start a **new** `/goal` (cannot resume). Token budget is not a constraint this week.

**15.268:** pairing-pole \(\prod(s-r)\) is a square; correct Aut is \(m_4^+(S)=\varepsilon(g,S)m_4^+(g(S))\); \(g=m_\sigma\circ\tau\in\mathrm{PSL}\) has \(\varepsilon=+1\); \(\nu=0\) on all \(|\kappa|=1\). Certified all 9900 finite \(|\kappa|=1\) at \(p=5\). Even \(\delta\) leftover (\(\mu=f_4\neq\mu_{\mathrm{part}}\) at \(p=5\); within-type even span at \(p=7\)).

## Goal and status

**Goal:** Settle MO 413935 (\(L=\lim\alpha_n\)). Binding leftover: residual **(i)** only.

**Status: OPEN — goal NOT met.** Residual (ii) ND is closed. Residual (i) / E(1) / \(L\) still OPEN. Soft-close forbidden. Do not flip `e1_closed_general` on census \(p\le7\). Do not ship a 15.xxx whose predicates stay False.

HEAD: 15.238–268 + STATUS/HANDOFF/LONG_HORIZON (this commit). Residual (i) still OPEN. \(\nu=0\) closed. Last residual-(ii) on `main` before this stack: `4d89353`.

| Flag | Value | Source |
|------|-------|--------|
| `residual_ii_a_ND_closed` | **True** | 15.237 |
| `residual_ii_b_ND_closed` | **True** | 15.236 |
| `residual_ii_full_closed` | **True** | 15.179 ∧ 15.237 ∧ 15.236 |
| `deep_s2_freeness_fail_k_ge_3p_ND_closed` | **True** | wires to 15.193 full |
| `gsum_disj_lb_proved_general` | False | 15.170 (hardcoded) |
| `residual_i_dual_eq_empty_proved_general` | False | 15.216 |
| `type_I_k_3p_minus_2_closed_general` | False | 15.170 |
| `e1_closed_general` | False | type_I ∧ deep_s2 ∧ bi-tight |

`/goal` hook is the live session plan (`goal/plan.md` under the session dir), **not** a repo `.md`.

## Residual-(i) hinge (unchanged)

Preferred: \(|\mu|\le 1/(2p)\) on \(|\kappa|=1\) (⇔ Gsum\(\ge-1/p\) ⇒ 15.176 Farkas).  
Sufficient: \(|\mu|\le 2/n\); envelope \(\mu\in\mathrm{conv}\{0,\mu_{\mathrm{part}},f_4\}\); dual-eq empty; \(G_+\succ0\) on \(\mathcal W_{++}^0\) (ker=sc) then Comm-repair dual; \(K_4\le n(15n-22)\).

## Exact 4-point sums (census; GPU + int16)

\(\nu\equiv 0\) and \(\mu\) constant on each CR class at \(p=5,7\). Integer sums \(\sum_{y\in\mathrm{Max}+}\prod_{i\in S}y_i\):

| \(p\) | \(\lvert\mathrm{Max}+\rvert\) | non-split sum | \(\mu\) |
|------|-------------------------------|---------------|--------|
| 3 | 12 | \(4\kappa\) | \(=\mu_{\mathrm{part}}\) (\(\varphi=-2\kappa\)) |
| 5 | 260 | \(8\kappa-2\varphi\) | \(=f_4=(4\kappa-\varphi)/(pn)\) |
| 7 | 11452 | \(240\kappa-10\varphi+96\,\kappa\,\varepsilon\,1_{\lvert\varphi\rvert=10}\) | CR 9 vs 10 split at extreme \(\varphi\) |

At \(p=7\): \(\varepsilon=+1\) on CR 9 (\(\mu=\pm109/2863\)), \(\varepsilon=-1\) on CR 10 (\(\pm61/2863\)). Mean of the split is the linear part \(240\kappa-10\varphi\). Interval envelope holds on **all** \(\lvert\kappa\rvert=1\) sets at \(p=3,5,7\).

Do **not** flip predicates on these formulas. Need a Max+-free \(F(\lambda)\) or a general envelope proof.

## Dead this continue (do not re-run)

| Attempt | Why dead |
|---------|----------|
| 4×4 Gram \(\{1,f_e,f_f,\chi_S\}\) | only \(\lvert m_4\rvert\le1-2/p\) |
| Affine-quadratic level sets on AG(2,5) | only the 60 linear halfspaces; 200 Max+ are not \(Q^{-1}(S)\) |
| CR-class master+diamond LP | 10 / 18 classes; HiGHS \(\max\lvert\mu\rvert=1\) (\(E_{\pm4p}\) survives) |
| Single IP-valency \(K_4\) | regular at \(p=5\) (valencies \(1,13,20,36,60\) + negs); **not** regular at \(p=7\) |
| Global mix \(\mu=(1-t)\mu_{\mathrm{part}}+t f_4\) | \(t\) depends on \((\kappa,\varphi)\) (\(p=7\): \(\approx0.20\) on \(\varphi=\pm2\), \(\approx0.04\) on \(\varphi=\pm6\)) |

## Earlier dead (still dead)

Fréchet / joint Fréchet; \(\lvert\mu\rvert\le\lvert f_4\rvert\) pointwise; \(\lvert\delta\rvert\le\mathrm{room}_\delta\) (false at \(p=5\) on \((\kappa,\varphi)=(-1,-6)\): \(\lvert\delta\rvert\approx0.039>\mathrm{room}=0.028\) while \(\lvert\mu\rvert=\lvert f_4\rvert\) still safe); W0\(\in\mathrm{im}(\mathrm{Gsum})\); Type I+switching \(G=\pi(G)\); low-degree Farkas \(y\); 3-wise interpolant; \(\kappa_B\perp E_{\pm4p}\); matching-PSD / one-edge \(\sum\mu\); max-\(\lvert\mu\rvert\) at extreme \(\lvert\varphi\rvert\); halfspace orbit as full \(\mu\); \(\chi(\)never-zero quadratic\()\); unsigned Per; Aut-line dim\(\le1\); 15.237 pair-span as Type I kill.

## Still viable

1. **Hull / envelope** on \(\lvert\kappa\rvert=1\) (15.240–241: ⇒ \(\lvert\mu\rvert\le2/n\le1/(2p)\)). Prefer hull \(\mu\in[\min(0,\mu_{\mathrm{part}},f_4),\max(0,\mu_{\mathrm{part}},f_4)]\) (census-stable; needs 0-vertex at p=7). Max-abs envelope is the weaker L∞ form. Census 0 viol at p=5,7.
2. **\(\lambda_*\) floor** on \(\mathcal W_{++}^0\) (15.242–243): certified min≥λ_* at p=3,5,7 with **mult n** at bottom (p=5 sharp). **Proved** ∑κ_C κ_B=(n+1)/4‖B‖² Max+-free; **proved** E[q²]≥λ_* ⇔ ∑ρ κ_B≥−6/n−1/(2p²). Attack the ρ lower bound by Aut-isotypic / orbit SOS (not scalar CS). Free-e_sc budget already Max+-free once ker=sc.
3. **Max+-free \(F(\lambda)\)** matching exact Max+ 4-point sums (p=7: two values at \(\lvert\varphi\rvert=10\)).
4. **\(K_4\le n(15n-22)\)** / Wick_hi / Path C \(\delta^2\).
5. Dual-eq Farkas **not** low-degree in \(f_e\).

## Rayleigh spectrum (census 15.242)

| p | dim | min Rayleigh | mult(min) | vs λ_* |
|---|-----|--------------|-----------|--------|
| 3 | 5 | 16 | 5 | ≫ |
| 5 | 65 | 80/13 | n=26 | **=** |
| 7 | 275 | 3072/409 | n=50 | > |

## Dead this continue (2026-08-13 post-compaction)

| Attempt | Why dead |
|---------|----------|
| Halfspace-only character sum as residual-(i) proof | m4_H exceeds 1/(2p) at p=5,7; full Max+ cancellation required (15.241 C) |
| CR-LP with Per eigenconstraint (probe) | Too expensive as written; prior CR master+diamond LP already allowed \(\lvert\mu\rvert=1\) without stronger SOS |
| Pointwise \(\lvert\delta\rvert\le\lvert g\rvert\) alone for envelope | Insufficient when μ_part and f4 have opposite signs |

## Geometry / IP (census)

- \(p=5\) Max+: 60 linear halfspaces + 200 other; 4 AG(2,5) line-intersection types. Other type is not affine-quadratic.
- \(p=5\) IP-regular; \(p=7\) not (row histograms differ).
- \(T\) definition (15.68, Max+-free): \((Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r})\).

## Scratch / caches

Scratch: `/tmp/grok-goal-51d7ac45c1de/implementer/` (`resi_cr_formula.py`, `resi_cr_lp.py`, `resi_mu_closedform.py`, `resi_fit_quad.py`, `resi_ip_int.py`, `resi_gram4.py`). May be deleted when the goal ends.  
Caches: `/tmp/maxplus_p5.npy` (260×26), `/tmp/maxminus_p5.npy`, `/tmp/maxplus_p7.npy` (11452×50), `/tmp/maxminus_p7.npy`.  
ProcessPool from `python -` stdin **fails** — write real scripts. 88 cores + idle V100; `W=86`.

## Shipped 15.244 (this resume)

- **Proved Max+-free** ∑_S φ κ_B = −n/4 ‖B‖² (star reduction Z=xxᵀ⊙B + CB=pB).
- **Proved** ∑ μ_part κ_B = (p²+3)/(4(p²−5)) ‖B‖²; γ=(3p²+5)/(2p²(p²−5))>0.
- **Proved** λ_* ⇔ ∑(m4−μ_part)κ_B ≥ −6/n−1/(2p²)−γ.
- Predicates still **False**. Residual (i) still OPEN.

## Shipped 15.245 (Z-frame)

- **Proved Max+-free** Z_y=yyᵀ−2P₊∈W++0, ‖Z‖²=n(n−2), Op=E[Z⊗Z], scheme=ker Op, avg Rayleigh=8(n−2)/(n−6)≥λ_*.
- **OPEN:** Op≽λ_*I. Dead: CS on ρ; s_max K₄; ρ linear only.
- Predicates **False**.

## Shipped 15.246 (edge Cov / R-path)

- Edge features f of Max+: constant-weight sphere; Cov=Op/2 on W++0.
- R-path ∑η² budget strictly weaker than Wick_hi for all p≥5.
- Predicates **False**. Prefer Cov gap or R-path L² over Wick_hi.

## Shipped 15.247–248 (m4_part + Comm dual forms)

- **15.247:** m4_part=μ_part+z star; R≤2p ⇔ ‖δ‖₂²≤room_δ^R; (3/2)(n−1) dual budget.
- **15.248:** We=½−1/(2p²(p²−2)); sum_ne^Comm closed and <(3/2)(n−1) Max+-free.
- **OPEN:** full D(C) inflation; ker=sc; or ‖δ‖₂²≤room; or Cov/|μ|.
- Predicates **False**.


## Shipped 15.249 (algebraic dual D_alg)

- **Proved Max+-free (any conference):** degree regular after Comm(scheme); Comm(diag) zero-diag with α=2p²/(p²−1); degree preserved; ⟨S,C⟩=n.
- **Closed forms:** We_alg=(p⁴−4p²+1)/(2p²(p²−2)); sum_ne¹=p²(p²+1)(p²−2)/(p⁴−4p²+1)−1; star=(p²+1)/den>0.
- **Paley/Weil:** far Sp·C·factor∈{1,p²,2p²−1}; Q=∑χ(γ(γ−α)(γ−β)), |Q|≤2p; m_p≥1−2p; t_ub=2(2p−1)/den.
- **Proved:** cost_D<2−α for free-e over sc, all primes p≥5.
- **OPEN for residual_i:** ker=sc **or** ‖δ‖₂²≤room_δ^R / |μ|≤1/(2p).
- Predicates still **False**. Soft-close forbidden.


## Shipped 15.250 (R-path as Max+ fourth moment)

- **Proved Max+-free:** Max+ closed under y↦−y ⇒ odd multilinear moments vanish.
- **Proved:** E[(y·z)⁴] = 9n²−10n + 24‖m4‖₂² (partition expansion + 2-design second moments).
- **Proved:** R-path ‖m4‖₂²≤n(n−2)/4 ⇔ E[s⁴]≤15n²−22n (⇔ residual-i via 15.217).
- **Census:** inequality holds at p=5,7 (Es4 room positive).
- **OPEN:** prove E[s⁴]≤15n²−22n for all primes p≥5.
- Predicates still **False**.


## Shipped 15.251 (Cy-identity; |μ|≤2/n path)

- **Proved Max+-free:** Per(C[S,S])=1 on |κ|=1; Cy-expansion (p⁴−1)m₄⁺+2φ=Ext on |κ|=1 (Paley).
- **Proved:** |m₄⁺|≤2/n and |m₄⁻|≤2/n on |κ|=1 ⇒ residual-i (via 15.176).
- **Census p=5,7:** m₄⁺=m₄⁻ on |κ|=1; max|μ₄|≤2/n (p=5: 3/65; p=7: 109/2863).
- **OPEN:** |μ₄|≤2/n for all primes p≥5 (control Ext, or hull, or Es4, or ker=sc).
- Predicates still **False**.

## Shipped 15.252 (Master/T²/Ext criteria)

- **Proved Max+-free:**
  - T²μ = 16(p²μ − κ) pointwise; 16κ+T²μ = 16p²μ (triangle alone vacuous).
  - |ρ|≤L_abs=(p−2)/(2p²) on |κ|=1 ⇒ |μ|≤1/(2p).
  - |Ext|≤2p²−4p+6 on |κ|=1 (with |φ|≤2(p−2)) ⇒ |μ|≤2/n.
- **Census p=5,7 (fresh):**
  - Envelope |μ|≤max(|μ_part|,|f4|) OK (0 viol); pure |f4| fails p=7 (117600 viol).
  - max|Ext|: p=5 24.8≤36; p=7 71.37≤76 (room ~4.6).
  - max|ρ|: p=5 0.0554≤0.060; p=7 0.0177≤0.051.
  - |16κ+T²μ| max: p=5 18.46 (=16p² M_cand sharp); p=7 29.85 ≤ 31.36 (2/n thr).
  - IP: p=5 regular (10 distances); p=7 not row-regular; s_max=(p−1)²−2 census-sharp both.
  - Es4: p=5 9261.5≤9568; p=7 31196≤36400.
- **Dead this resume:** crude |T²|≤16(n−2)/n as general maj (false p=5; and triangle is vacuous with master rewrite).
- Predicates still **False**. Soft-close forbidden.

## Shipped 15.253 (Wick-reflection — preferred hinge)

- **Proved Max+-free:** |ρ_f4|≤L_abs for all primes p≥5
  (maj (3p²−8p+1)/(p²n); g(p)=p³−8p²+17p−4, g(5)=6, g'>0 on [5,∞)).
- **Proved criterion:** |ρ|≤|ρ_f4| ⇔ μ∈[f4, f4^♯] with f4^♯=2κ/p²−f4
  ⇒ |ρ|≤L_abs ⇒ |μ|≤1/(2p) ⇒ residual-(i) Farkas.
- **Proved:** hull ⇒ |ρ|≤L_abs (triple max(1/p²,|ρ_part|,|ρ_f4|)).
- **Census p=5,7:** reflection viol=0; t=(μ−f4)/ρ_f4 ∈[−1.703,−0.098]⊂[−2,0].
  On Aut classes without CR-split, t is constant (e.g. p=7 (−1,2): t=−1393/818).
- **Strictly weaker than |μ|≤|f4|** (false at p=7); Wick-centred interval is the right target.
- Predicates still **False**. Soft-close forbidden.

## Shipped 15.254 (T m₄± forms; Paley C∼−C)

- **Proved Max+-free (any conference + π):**
  - \(T m_4^+=4p\,m_4^+-4\kappa/p\), \(T m_4^-=-4p\,m_4^-+4\kappa/p\)
  - \(T\mu=2p(m_4^+-m_4^-)\), \(\nu=\frac12(m_4^+-m_4^-)\)
- **Proved Paley:** \(D P^\top C P D=-C\) (P=mult-by-nonsquare, \(D_\infty=-1\));
  \(\mathrm{Max}- = D\cdot(\mathrm{Max}+\circ\pi^{-1})\);
  \(\mu(S)=\frac12\bigl[m_4^+(S)+\chi_D(S)\,m_4^+(\pi^{-1}S)\bigr]\).
- **(κ,φ) under π on |κ|=1:** ∞∉S ⇒ fixed; ∞∈S ⇒ (κ,φ)↦(−κ,−φ); \(f_4(S)=\chi_D f_4(\pi^{-1}S)\).
- **Census:** m₄⁺=m₄⁻ and Tμ=0 on all |κ|=1 at p=5,7.
- **OPEN:** prove m₄⁺=m₄⁻ on |κ|=1 generally (⇔ \(m_4^+(S)=\chi_D m_4^+(\pi^{-1}S)\)), and/or reflection (15.253).

## Next concrete steps

1. **Prefer:** reflection |ρ|≤|ρ_f4| on |κ|=1 (15.253), or m₄⁺=m₄⁻ + further bound.
2. Alts: hull / |Ext| maj / Es4 / ker=sc.
3. Wire `gsum` / `type_I` / `e1` **only** via real imports if proved.
4. Do **not:** soft-close; |μ|≤|f4| as general; census-only flip.

**Do not:** ship a False-predicate 15.xxx; flip e1 on census; soft-close; re-use |μ|≤|f4| as general; re-thrash dead list.

## Suggested skills

`agent-cost-optimization` · `graph-engineered-completion` · `goal-verifier` · `verification-before-completion` · `handoff` · `session-handoff-packager` · `scientific-critique` · `grill-me` · `self-refine-loop` · `research` · `arxiv` · `litreview` · `use-available-compute` · `openai-referee`
