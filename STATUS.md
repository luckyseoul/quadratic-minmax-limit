# Status (2026-08-06)

**\(L=\lim\alpha_n=\tfrac12\) is OPEN** (candidate denseness path blocked by a named hinge).

| Claim | Status | Reference |
|-------|--------|-----------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) | CLOSED | `solution.md` |
| \(\rho=1\) on Paley \(n=p^2+1\) | CLOSED | `PROOF_rho_eq_1.md` |
| Bi-tight majorization algebra (15.167) | CANDIDATE | `prop15167.py` |
| Type I freeness-fail ND (residual i) | OPEN | needs Gsum disj LB |
| Deep freeness-fail ND (residual ii) | **CLOSED** | freeze-to-tight 15.179 |
| E(1) / \(L=\tfrac12\) | **OPEN** | residual (i) hinge |
| Path-C residual / 16N | OPEN optional | not required for denseness path |
| Prize acceptance | OPEN | X + GitHub; Paata AI-test |

### Fatal gap (one sentence)

Residual **(i)** dual-equality Farkas (e∉G, Type I \(k=3p-2\)) needs pointwise disj \(\mathrm{Gsum}>\mu_*=(6/p-4)/k\) (Prop **15.176**); sufficient \(\mathrm{Gsum}\ge-1/p\). Residual **(ii)** closed by freeze-to-tight (**15.179**), no Gsum.

### Remainder progress (15.172–181)

- **15.179 residual (ii) CLOSED:** dual two-level freeness-fail affine \(\Rightarrow S_H\equiv3\Rightarrow k=3p-1\); impossible for \(k\ge3p\); fail-eq empty under bi-tight.  
- Avg disj Gsum \(=2/(n-3)\); \(G_0\) PSD; **15.176** μ_* / −1/p sufficiency.  
- **15.177–178:** |μ₄| hinge form; star identity; dual-eq \(n_d\le1\) kill; \(n_d=2\) wedge kill for \(p\ge7\).  
- **15.180:** dual-eq Q_pairs=30−6p−24/p<0; open dual-eq core is \(n_d\ge2\) after PSD/score filters.  
- **15.181:** vertex-star \(\sum_{e\ni i}f_e=\pm p\) ⇒ \(\sum_{e\ni i}\mathrm{Gsum}_{ef}=2\); Max+-free κ-counts through edge \(n_3=(n-2)(n-6)/8\), \(n_1=3(n-2)^2/8\); matching \(n_d=2\) PSD floor −4 (harder than wedge); p=5 dual-eq k-sparse linear box empty (census).  
- **15.182:** dual-equality normal form \(x=\alpha\mathbf{1}-2e_*+\kappa\) with \(\alpha=6/(pn)\), \(\kappa\in\ker(G_+)\cap\ker(G_-)\), \(\kappa_e=2-\alpha\), box on \(\alpha+\kappa\); particular solution Max+-free; p=5 ker-box LP infeasible (census). General ker-box obstruction open.  
- **15.183:** Max+⊥Max− (any symmetric conference); \((G_+G_-)_{ee}=n/(2p^2)\); ker(Gsum)⊥1 automatic; matching \(n_d=2\) PSD \(a+b\ge-2\sqrt{2+c}\); binary dual-eq form. p=5 max \(\kappa_e\) under ker+box ≈0.811<2−α.  
- **15.184:** Max+-free \(T^2\kappa=-24\varphi+48\kappa\) on \(|\kappa|=1\) (C² reduction + 48-labeling); \(|T^2\kappa|=24|\varphi-2\kappa|\).  
- **15.185–187:** Paley \(|\varphi|\le2(p-2)\) all odd primes (Auer–Top supersingular\(\Rightarrow\)double residue + Hasse ladder). Global \(T^2\kappa=-24\varphi+48\kappa\) on **all** 4-sets (any conference, 64-labeling). \(T^2\varphi=4(p^2+3)(\varphi-2\kappa)\) (Paley census p=3,5,7). \(\mu_{\mathrm{part}}=[(p^2-1)\kappa-2\varphi]/(p^2(p^2-5))\) solves master when \(T^2\varphi\) holds; majorant \(\le1/(2p)\) for \(p\ge5\).  
- **15.188:** Target correction — \(|\mu_{\mathrm{actual}}|\not\le|\mu_{\mathrm{part}}|\) pointwise (p=7 census: ~29k violations). Viable sufficient bound \(|\mu|\le2/n=2/(p^2+1)\) (\(\le1/(2p)\) for \(p\ge5\)). p=5: actual\(=(4\kappa-\varphi)/(pn)\), max \(3/65\). p=7: max \(109/2863<1/14\) and \(<2/n\).  
- **15.189 (Max+-free):** \(1^\top y=(p+1)y_\infty\) on Max+; \(E_\pm[y_iy_j]=\pm C_{ij}/p\); tight frame; adjacent Gsum=0 from π; G+ 6×6 PSD ⇒ \(|\mu|\le1-2/p\) on \(|\kappa|=1\) (**too weak** for residual i).  
- **15.190 (Max+-free):** scheme-ker \(\kappa_{ij}=f_i+f_j\) (\(\sum f=0\)) lies in \(\ker(Gsum)\); scheme-ker max \(\kappa_e=\alpha(n-2)/2=3(n-2)/(pn)<2-\alpha\) for all \(p\ge3\). Full dual-eq ker-box empty at p=5,7 (census: max \(\kappa_e\approx0.811,0.593<2-\alpha\)); **feasible at p=3**. 3-point moments vanish (certified p=3,5).  
- **15.191 (Max+-free partial):** Derangement permanent of \(C[S,S]\) equals 1 on \(|\kappa|=1\) (64-exhaust); star-sum \(\sum_s\prod_{i\neq s}C_{is}=0\) on \(|\kappa|=1\); Cy-expansion size1+size2 \(=-2\varphi\) (any conference \(C^2\) + Paley \(\pi\)); envelope \(|4\kappa-\varphi|/(pn)\le2/n\le1/(2p)\) for \(p\ge5\). **Correction:** \(|\mu|\le|f_4|\) fails at p=7 (many classes; f4 not a pointwise majorant); viable target remains \(|\mu|\le2/n\) (census p=5,7) or \(\le1/(2p)\).  
- **15.192 (Max+-free):** Gsum diag\(=2\); row sum\(=n\); avg disj Gsum\(=2/(n-3)\). Aut_e averaging: dual-eq feasible iff Aut_e-invariant dual-eq feasible. \((3/2)\cdot\)scheme-max \(<2-\alpha\) for all \(p\ge5\). **Census Aut_e ker-box:** p=3 feasible (\(\max\kappa_e=14/5\)); p=5,7 empty (\(\max=369/455\), \(11736/19775\); ratios to scheme \(41/28\), \(163/113\), both \(<3/2\)).  
- **Still open (residual i):** Max+-free \(\max\kappa_e\le(3/2)\cdot\)scheme-max (or any bound \(<2-\alpha\)) for all \(p\ge5\), **or** \(|\mu|\le2/n\) (or \(\le1/(2p)\)) on \(|\kappa|=1\).  
- `gsum_disj_lb_proved_general()=false`; residual (ii) closed; **E1/L OPEN**.

### Short package

`evidence/share/denseness_path_package.md`

### Required opens (denseness prize path)

1. **Math:** Prove disj Gsum \(\ge-1/p\) (or any \(\mu>\mu_*\)) for residual **(i)** only (15.176–178); residual (ii) done (15.179).  
2. **Predicates:** Flip `gsum_disj_lb_proved_general` → residual (i) → E1 → L only after (1).  
3. **AI-test:** Re-run on short package only after L predicates CLOSED.  

**Non-required:** Path-C / 16N / Hypothesis H.  
**Out of agent control:** Ping Paata (user).  

**Current:** residual (ii) closed; residual (i)/E1/L open; claim **not** asserted.
