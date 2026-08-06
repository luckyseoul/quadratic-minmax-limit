# Status (2026-08-06)

**\(L=\lim\alpha_n=\tfrac12\) is OPEN** (candidate denseness path blocked by a named hinge).

| Claim | Status | Reference |
|-------|--------|-----------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le1/2\) | CLOSED | `solution.md` |
| \(\rho=1\) on Paley \(n=p^2+1\) | CLOSED | `PROOF_rho_eq_1.md` |
| Bi-tight majorization algebra (15.167) | CANDIDATE | `prop15167.py` |
| Type I / deep freeness-fail ND | OPEN | needs Gsum disj LB |
| E(1) / \(L=\tfrac12\) | **OPEN** | gap below |
| Path-C residual / 16N | OPEN optional | not required for denseness path |
| Prize acceptance | OPEN | X + GitHub; Paata AI-test |

### Fatal gap (one sentence)

Dual-equality Farkas needs pointwise disj \(\mathrm{Gsum}>\-2/p\) (Prop **15.172** threshold); only wedge-zero + avg \(2/(n-3)\) are proved Max+-free. Candidate \(-12/(pn)\) fails at \(p=3\) and is not scheme-justified (**15.158**).

### Remainder progress (15.172–173, not a close)

- Avg disj Gsum \(=2/(n-3)\); triangular \(G_0\) spectrum/PSD proved.  
- Farkas threshold: \(\mathrm{Gsum}>\-2/p\) iff \(H>\-1/p\) (\(H=\mathrm{Gsum}/2\)).  
- **Prop 15.173 vector structure (Max+-free proved):** \(\langle\xi,v_e\rangle=1/p\), star simplices on \(w_e\), \(K=H-J/p^2\succeq0\), \(K\mathbf1=0\).  
- Weak CS \(H\ge2/p^2-1\) proved but **too weak** for Farkas.  
- Census: \(p=3\) tight \(H=-1/p\); \(p=5\) min \(H=-3/65>-1/p\).  
- **Still open:** \(H_{ab}\ge-1/p\) for all disj pairs, all primes \(p\ge5\).  
- `gsum_disj_lb_proved_general()=false`; **L OPEN**.

### Short package

`evidence/share/denseness_path_package.md`

### Required opens (denseness prize path)

1. **Math:** Prove \(H_{ab}\ge-1/p\) for disj edges, all primes \(p\ge5\) (Max+-free; Prop 15.173 frame), **or** alternate residual (i)/(ii).  
2. **Predicates:** Flip `gsum_disj_lb_proved_general` → residual i/ii → E1 → L only after (1).  
3. **AI-test:** Re-run on short package only after L predicates CLOSED (expect “essentially correct”).  

**Non-required:** Path-C / 16N / Hypothesis H.  
**Out of agent control:** Ping Paata (user).  

**Current:** (1)–(3) open; claim **not** asserted.
