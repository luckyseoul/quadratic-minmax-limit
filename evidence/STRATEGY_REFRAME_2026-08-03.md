# Strategy reframe (2026-08-03) — structure, not Prop 15.159

**Source:** external cold read of manuscript + handoff (user message).  
**Status:** Adopted as P0 research stance. **L still OPEN.**  
**Do not:** ship Prop 15.15x that only re-encodes δ²≤room.

---

## Diagnosis (agree)

Stuck **not** for lack of compute. Stuck because natural analytic attacks are already compressed into **documented obstructions**. Soft analysis exhausted; Paley/conference pushed far; many equivalent residual forms; dozens of inequalities killed. What remains is **one structural statement** (localized residual), which is better than most open problems.

**Wrong question (current thrash mode):**  
> How do I prove δ² ≤ room_hyp/24?

**Right question:**  
> What larger structure forces δ² ≤ room_hyp/24 (or makes free residual mass impossible except equality cases)?

δ should not be a free parameter to bound; it should be **rigidified away** (Hoffman/Delsarte/design/scheme equality-case style).

**Stop:** adding lemmas with diminishing returns unless they introduce a **new structural identification**.

---

## Six attack vectors — map to project state

### 1. Association schemes / Bose–Mesner

**Intent:** If P⊙P (or T / Norton product) lives in a BM algebra, λ₂ is representation-theoretic.

| Already known | Status |
|---------------|--------|
| Max+ **not** IP-association scheme (p=5) | **Blocked** pure R_s (15.158) |
| class_key BM | **Dead carrier at p=7** (nul=0 but δ≠0; F19) |
| Strict Aut(C) G BM on 4-sets | **Open carrier**: dim E_{4p}^G = 0,2,7 at p=3,5,7; residual lives in G-space at p=5,7 (15.134) |
| Conference graph / Johnson / signed Johnson for T | Partially used; full BM pin of PopP bulk spectrum **not** closed |

**Leverage:** Do **not** re-try IP-scheme of Max+. Do try: (i) coherent configuration finer than IP but coarser than free Aut; (ii) BM of **edge** / **pair** operators where Φ and P⊙P act; (iii) known eigenvalues of Schur squares of design projectors in the literature (Fickus/ETF residual-Gram pattern already cited 15.101 — extend, don’t re-name).

### 2. SDP dual certificates

**Intent:** λ₂(P⊙P) or Φ|Z ≤ 16 as SDP optimum → exhibit dual feasible certificate.

| Already known | Status |
|---------------|--------|
| GW / rank-1 SDP for E(1) deloc | **Dead** (F6; SDP·α short of discrete) |
| SDP_+ ≥ n√(n−1) for all Seidel | **False** (counterexample in solution.md) |
| Dual certificate for residual orth / 16N | **Not shipped** as general-p dual |

**Leverage:** Formulate residual as explicit SDP (or SOS on 4-set functions / matrices on Z). Search for **closed dual** (polynomial in p, conference C only) that proves λ_max(Φ|Z)≤16 or orth≤room_hyp. Avoid redoing GW for maxcut.

### 3. Representation theory (S_n / Specht / PSL)

**Intent:** T or Veronese/Norton acts by scalars on irreps → bound is character arithmetic.

| Already known | Status |
|---------------|--------|
| PSL min irrep ⇒ mult(λ₂)≥d−1 | **Proved** (15.98) — used for majorization, not full spectrum |
| Aut-Schur on maximizer locus of Γ | Structure (15.97); mult=d not closed |
| Full Specht decomposition of 4-set space under S_n | **Not done** (probably too large); restrict to G=Aut(C) isotypics |

**Leverage:** Decompose **E_{4p}^G** (dim 2 at p=5, 7 at p=7) into G-irreps; compute T / m₄ pairing as characters. If residual subspace is a single known irrep with closed matrix coefficient, δ² becomes formula not estimate.

### 4. Flag algebras (signed)

**Intent:** Global 4-set averages under local configuration constraints.

| Already known | Status |
|---------------|--------|
| Hamming Delsarte / moment LP on W_k | **Too weak** for residual (15.124) |
| Spherical LP / μ_G4≤1, 1/h₄ | **Dead or false** (15.157) |
| Signed flag algebra | **Not attempted** |

**Leverage:** Speculative. Only if scheme path fails. Objects are ±1 Seidel / Max+, not ordinary graphs — need signed flag machinery. High setup cost; no prop mill without a toy inequality first (e.g. recover known p=5 equality).

### 5. Finite geometry / character sums

**Intent:** Residual is a theorem about AG(2,p), lines, quadrics, χ after translation.

| Already known | Status |
|---------------|--------|
| Halfspace ρ=1 from F_{p²} linear form | **Proved** |
| Affine+Frob orbit ⊊ Max+ (p=5: 60/260) | **F18** — incomplete as full Max+ |
| G-orbits on 4-sets; Gauss m₄ on orbits | **Named open** (15.134.E) — not closed |
| Interval ρ character sum (E2) | Partial; not residual |

**Leverage:** Highest continuity with 15.134: **Max+-free m₄ on G-orbits via Weil estimates**, then project to E_{4p}^G. That is “geometry forces residual,” not another CS bound.

### 6. Forget conference: limit objects of all optimal A

**Intent:** Characterize graphon / operator limits of Φ-minimisers; conference as consequence.

| Already known | Status |
|---------------|--------|
| Graphon USC for Φ | **Fails** (recorded; random ±1) |
| Soft multipartite / Hadamard existence | **Dead** (solution §) |
| E(1) as product ρ·‖A‖_op rigidity | Open; deloc barrier |

**Leverage:** Still valid dual settlement (Prop 6.2 needs lim along dense family). Do not use fake graphon uniqueness. Possible: operator-norm + boolean cut rigidity without graphons (e.g. only conference achieve near-min Φ asymptotically).

---

## Research protocol (next month stance)

1. **No new Prop 15.15x** unless it **identifies** the residual operator with a named structure (BM element, SDP dual, irrep, geometry sum).  
2. Spend effort on **identification**, not inequality inventing.  
3. Equality case at p=5 (δ² = room_hyp/24) is a **clue**, not a nuisance — any structure must force saturation only there or explain the drop at p=7.  
4. Parallel **E(1)** (no-descent / k_⋆) remains alternate Main path; same rule: structure over prop count.

---

## P0 graph update (nodes)

| Node | Status |
|------|--------|
| Residual inequality forms | **Localized** (many equivalent) |
| Bound δ by analysis | **Diminishing returns / thrash risk** |
| Structure identification | **Open — preferred** |
| BM (Max+ IP) | **Dead** |
| BM (G on 4-sets) | **Open carrier** |
| SDP dual for 16N/orth | **Open** |
| G-irrep / character | **Open** |
| Signed flag | **Untried / low priority** |
| Weil m₄ on G-orbits | **Open (named)** |
| Limit object of all optima | **Open alternate** |

---

## Suggested skills for next agent

`graph-engineered-completion` · `use-available-compute` · `agent-cost-optimization` · `goal-verifier` · `handoff`  
Optional: lit/arxiv search for BM eigenvalues of design Schur squares / signed flag algebras — **before** coding.
