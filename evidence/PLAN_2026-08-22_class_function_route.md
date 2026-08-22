# Route of attack: the class-function / character-table reduction for R1

Date: 2026-08-22.  Target: R1 = `‖δ‖² ≤ n/12` ⟺ `Φ_δ ⪰ −(2n+20)/(n−6)·I`
⟺ `λ_min(Φ) ≥ 6`, which closes leftovers 1 AND 3 (15.595).
**Nothing here is proved yet.  This is a plan plus two diagnostics.**

## Why a new route is needed (what the dead ends have in common)

Killed so far: linear 4-point LP; linear 6-point LP; level-4 SoS; Delsarte
2-design + min distance; degree escalation of the contraction kernel;
`(12+ε)n²` majorants; uniform `C/p⁴` scaling; pointwise SOS; equal-density
shortcut; Aut-invariant master + `|m₄|≤1`.  **Every one is a moment
relaxation** — it replaces Max+ by a moment sequence, and the feasible
region always contains fake sequences that meet all constraints without
being sets of ±1 eigenvectors.  Adding moments does not help: degree 6
gained nothing, established twice from independent directions (15.594 and
the 6-point box LP).  A closing argument must use something no moment
relaxation sees.

## The reformulation

For W ∈ Z, `⟨B̃_y, W⟩ = yᵀWy` (15.597), and `π(g)B̃_y = B̃_{g·y}`, so

    A_y(g) := ⟨B̃_y, π(g)B̃_y⟩ = s(y, g·y)² − 2n,     s(y,z) := y·z

is a positive-definite function on G, and the constituent energies are its
Fourier coefficients:  `‖P_c B̃_y‖² = (deg_c/|G|)·Σ_g χ_c(g)·A_y(g)`.

**Diagnostic 1 (done, decisive).**  Is the spread pointwise or only on
average?  `scripts/constituent_energy_diag.py`:

| p | per-vector `‖P_c B̃_y‖²` | verdict |
|---|---|---|
| 5 | CONSTANT across all 260 vectors (160, 288, 176 — exact integers) | pointwise holds |
| 7 | VARIES: min 0, max 1333 on the mult-100 constituent | **pointwise FAILS** |

p=5 is special (Max+ is one G-orbit there; at p=7 it cannot be, since
|G|/|Max+| = 117600/11452 ∉ ℤ).  So a single-vector argument is DEAD —
averaging over the ensemble is essential.  Total energy per y is constant
= n(n−2) at both primes, as proved (15.593 A).

## The route: Γ is a class function, so R1 is ~q numbers, not C(n,4)

Define the **ensemble correlation**

    Γ(g) := E_y[ s(y, g·y)² ] − 2n .

Then `λ_c = (1/|G|)·Σ_g χ_c(g)·Γ(g)` (using mult_c = deg_c, Z
multiplicity-free, 15.589 B).

**Key: Γ is a CLASS FUNCTION.**  Max+ is G-invariant, so for any h,
Γ(hgh⁻¹) = E_y[s(y, hgh⁻¹·y)²] = E_y[s(h⁻¹y, g·h⁻¹y)²] = Γ(g).
PSL(2,q) has ≈ q+4 conjugacy classes, so **Γ is determined by ≈ q numbers
instead of C(n,4) ≈ n⁴/24 four-set moments** — a compression by ~n³.

Then `λ_min(Φ) ≥ 6` becomes: for every irreducible χ_c appearing in Z,

    Σ_K (|K|/|G|)·χ_c(K)·Γ(K)  ≥  6 ,

a finite inequality against the **classical, fully explicit character table
of PSL(2,q)** (trivial, Steinberg, (q−3)/2 principal series, (q−1)/2
discrete series).

## Why this can work where the relaxations could not

* It uses the group action on Max+, which no moment LP encodes.
* The unknowns collapse from ~n⁴ to ~q, and the coefficients (characters)
  are known exactly in closed form for all q.
* Γ(K) plausibly depends only on the fixed-point structure of g on P¹(F_q)
  (elements are classified by trace: split / non-split / unipotent /
  central), which is exactly how Paley/Max+ structure is organised —
  so a closed form for Γ(K) is a realistic target, not a fit.
* `Γ(g) = Σ_{x,x'} d_x(g) d_{x'}(g) m₄(x, g(x), x', g(x'))` expresses Γ
  through m₄ but only via ~q aggregates, so the δ-part enters as ~q
  numbers rather than an unbounded function — the object that made every
  LP feasible-but-wrong.

## Concrete next steps, in order

1. **Compute Γ(K) at p=5, 7** for every conjugacy class, from the exact
   ensembles (cheap; both fit in memory).  Verify it is genuinely a class
   function (invariance check, not assumed).
2. **Reconstruct λ_c from Γ(K) + the PSL(2,q) character table** and check
   it reproduces the measured spectra (6.1538/11.0769/13.5385 at p=5;
   7.511/8.2152/8.9193/9.8582/10.5623 at p=7).  This validates the whole
   pipeline before any general-p claim.
3. **Hunt a closed form for Γ(K)** as a function of the class parameter
   (trace / fixed-point count).  Fit at p=5,7 then PREDICT p=11 and check
   against the stored 37.4M ensemble — out-of-sample, per the
   pre-asymptotic-fit discipline.  A form that survives p=11 is evidence;
   one that fails kills the route cheaply.
4. If a closed form holds, **prove `Σ_K (|K|/|G|)χ_c(K)Γ(K) ≥ 6`** using
   the explicit character values.  This would close R1, hence leftovers
   1 and 3 simultaneously.

## R2 (leftover 2)

Max− Walsh containment, exact at p=3,5,7,11 (15.596).  Documented trap:
at p=11 `rank(B_U) = 60`; dim H is **60, not n/2=61** (15.598 census).
A general-p proof via "B_U has full rank n/2" is FALSE.

**15.598 (proved):** square-direction affine lines \(L=a+\mathbb F_p b\),
\(\chi(b)=1\), force \(\sum_{S=\{\infty\}\cup L} y=0\) on Max−, by the
elementary sum \(\sum\chi_p(x(x+\delta))=-1\) and
\(\sum_{i\in L}C_{ij}=-\chi(b)\). Pair-slice \(U\) is the xor cut of
\(\mathbb F_2^n\). Walsh ∀p reduces to affine_span(\(U\))=H∩{ℓ=c}.

**15.599 (proved pin, Walsh still open):** rank(SSᵀ)=n/2−1,
rank(S)≤n/2, so rank∈{n/2−1,n/2}; equality n/2 at p=3..37.
Antipodes fix the p=11 dim-60 half-ensemble. Aut_e-irreducibility
is false (referees: do not reopen). residual_ii still False.

**15.600 (proved):** rank(S)=n/2 for every odd prime. dim H0=n/2
is a theorem. Walsh remaining: U spans the xor-hyperplane of H⊂H0.

---

# Step 1 EXECUTED (2026-08-22) — the compression is real and Γ_δ is quantized

`scripts/gamma_class_function.py`, p=5 (|G⁺| = 31,200 enumerated by BFS on
signed permutations, |Max+| = 260):

* **Γ is a class function** — 0/300 conjugation violations (tested, not
  assumed).
* **Γ takes only 14 distinct values** over 31,200 group elements.  Γ(e) =
  624 = n(n−2) ✓.  It does NOT reduce to the fixed-point count alone
  (fix=0 → 3 values, fix=1 → 4, fix=2 → 4), so the genuine conjugacy class
  is needed — consistent with ≈ q+4 classes.

## Theory sharpened while running it

    Γ(g) = tr(Φ·π(g))        ⟹      Γ = Σ_c λ_c χ_c ,

so **the λ_c are exactly the Fourier coefficients of Γ** and λ_c = ⟨Γ,χ_c⟩
by orthogonality.  Splitting with Theorem A* (Φ = λ̄I + Φ_δ):

    Γ(g) = λ̄·ψ_Z(g) + Γ_δ(g),      Γ_δ(g) := tr(Φ_δ π(g)) = Σ_c(λ_c−λ̄)χ_c(g)

and ψ_Z — the character of Z — is **known in closed form** (15.589 B:
Z = W_e ⊕ (q−9)/8 principal series).  So the explicit half of Γ is already
determined for all p, and every open bit sits in Γ_δ.

## Γ_δ is quantized (new, exact)

At p=5, Γ_δ takes **12 distinct values, all integer multiples of a single
unit** u := ‖δ‖²/24 = 64/65:

    Γ_δ / u  ∈  {−10, −9, −6, −4, −1, 0, 2, 6, 9, 10, 12, 24}

with Γ_δ(e) = 0 — matching the proved corollary tr(Φ_δ) = 0 — and
max Γ_δ = 24u = ‖δ‖² = 1536/65 exactly.  Γ_δ is nonzero on 74.8% of G⁺.

This is the structure a closed form needs: the open content of R1 is a
class function taking ~12 integer multiples of one scale, and
λ_min(Φ) ≥ 6 ⟺ λ̄ + min_c ⟨Γ_δ, χ_c⟩ ≥ 6 — a finite pairing of an
integer-valued class function against the classical PSL(2,q) character
table.  **Nothing proved yet**; the quantization is one prime (p=5), and
per the pre-asymptotic discipline it must be reproduced at p=7 before any
closed form is proposed.

## Revised next steps

1. Repeat at p=7 (|G⁺| ≈ 117,600, |Max+| = 11,452): confirm class-function
   property, count distinct Γ_δ values, and test whether the quantization
   unit is again ‖δ‖²/24.  **This is the gate** — if the integrality fails
   at p=7 the closed-form hunt is off.
2. If it holds, identify which conjugacy classes carry which multiples
   (by trace / split / non-split / unipotent type).
3. Then pair against the explicit character table and bound ⟨Γ_δ, χ_c⟩.

---

# Step 1 GATE RESULT (p=7): quantization is DEAD; the compression SURVIVES

`scripts/gamma_class_p7_gate.py`, p=7 (3,254 distinct elements sampled by
random walk on G⁺, |Max+| = 11,452, dim Z = 275).

Sanity anchors all pass: Γ(e) = 2400 = n(n−2) ✓, ψ_Z(e) = 275 = dim Z ✓,
Γ_δ(e) = 0 ✓ (the proved tr(Φ_δ) = 0).

## KILLED: Γ_δ is not quantized

With u := ‖δ‖²/24 = 0.434326 (the p=5 unit), Γ_δ/u at p=7 reads

    −44.81, −15.03, −14.79, −10.07, −3.34, −0.64, 0,
      1.18,   4.23,   4.42,   7.66,  14.79, 30.66

— **not integers.**  The p=5 integrality (exact multiples of 64/65) was an
artifact of Max+ being a single G-orbit at p=5; at p=7 it cannot be, since
|G⁺|/|Max+| ∉ ℤ.  So the "find a closed form for Γ_δ via integrality"
shortcut is off.  Do not reopen it.  (No other unit was fished for — that
would be exactly the pre-asymptotic fitting error this discipline exists
to prevent.)

## SURVIVES: the class-function compression, and it is stronger than expected

| p | \|G⁺\| | distinct Γ values |
|---|---|---|
| 5 | 31,200 (full enumeration) | **14** |
| 7 | 117,600 (3,254 sampled) | **15** (lower bound — sampled) |

The group grows 3.8× and the number of distinct Γ values goes 14 → 15.
PSL(2,49) has ≈ q+4 = 53 conjugacy classes, so Γ is *far* more compressed
than one-value-per-class: many classes share a Γ value.  If that
persists, Γ is carried by ~15 numbers largely independent of p.  Caveat:
the p=7 count comes from a 2.8% sample, so 15 is a floor, not a census.

## What is theorem (not empirical) and remains the frame

* `Γ(g) = tr(Φ·π(g))`  ⟹  `Γ = Σ_c λ_c χ_c`: the λ_c ARE Γ's Fourier
  coefficients, λ_c = ⟨Γ, χ_c⟩.
* `Γ = λ̄·ψ_Z + Γ_δ` with ψ_Z known in closed form (15.589 B), so the
  explicit half of Γ is solved at every p and all open content is Γ_δ.
* `λ_min(Φ) ≥ 6  ⟺  λ̄ + min_c ⟨Γ_δ, χ_c⟩ ≥ 6`.

## Revised next step

The closed form must come from structure, not integrality.  The concrete
lead: determine which conjugacy classes share a Γ value, and whether the
partition is by the classical PGL(2,q) element types (split / non-split /
unipotent / central, by trace).  If Γ is constant on those four families
plus a few strata, it is carried by O(1) numbers and the character pairing
is a finite closed-form computation for all p.  That is testable at p=5
by full enumeration (classes are computable there) before any claim.

---

# Step 2 (p=5, full class enumeration): what makes Γ collapse

`scripts/gamma_conjugacy_classes.py`, p=5: all 31,200 elements of G⁺,
conjugacy classes computed by orbit under conjugation (not assumed).

**36 conjugacy classes → 14 distinct Γ values.**  The mechanism:

## Γ(−g) = Γ(g) — a theorem, and the first factor of 2

s(y, (−g)·y) = −s(y, g·y) and Γ squares it, so Γ is invariant under the
central sign twist.  Confirmed in the data: every Γ value is carried by a
class pair identical in (order, #fixed points) and differing only in sign
pattern (#neg = 0 vs 12), e.g. (12,fix2), (10,fix1), (4,fix2), (2,fix6).
The extreme case is Γ(−I) = Γ(I) = 624 = n(n−2).
So **Γ descends to G⁺/{±1}: 36 → 18 classes**, then four further
coincidences bring it to 14.

## Γ is NOT a function of (order, fixed-point count)

Two distinct ±-class-pairs both have order 10 and exactly 1 fixed point,
yet Γ = −17.2308 and Γ = +8.6154.  So the coarse invariants do not
determine Γ; it needs the genuine class datum (for semisimple elements of
PGL(2,q), the eigenvalue ratio up to inversion).  This is consistent with
the principal/discrete-series split of the character table, and it means
a closed form must be indexed by that ratio, not by order or fixed points.

## Where the route stands

Established as theorem, independent of the failed quantization:
* Γ(g) = tr(Φ·π(g)) ⟹ λ_c = ⟨Γ, χ_c⟩,
* Γ = λ̄·ψ_Z + Γ_δ with ψ_Z known in closed form,
* Γ(−g) = Γ(g), so everything lives on G⁺/{±1},
* λ_min(Φ) ≥ 6 ⟺ λ̄ + min_c ⟨Γ_δ, χ_c⟩ ≥ 6.

Established empirically: the compression is severe and stable —
36 classes → 14 values at p=5; ≥15 values at p=7 out of ≈53 classes.

Remaining for a closed form: **Γ(K) as a function of the eigenvalue-ratio
class parameter and q.**  With that, the pairing against the explicit
PSL(2,q) character table is a finite computation valid for all p, and R1
(hence leftovers 1 and 3) closes.  Without it, this route is a better
frame but not yet a proof.  Any proposed form must be fitted at p=5,7 and
then PREDICT p=11 out of sample before it is claimed.

Clarification vs the opening claim: `λ_min(Φ)≥6` is leftover 1 as a
whole (QVAR on `W_e` **and** the principal floor).  It is not Type I
leftover 3.  15.595 is a variance sufficient condition, not that unit.

---

# Step 3 (p=5 fit, p=7 gate): Γ on PSL is a class-parameter function; the O(1)-in-p split formula is dead

`scripts/gamma_class_parameter.py`.  Full Aut(C): p=5 |G⁺|=31,200; p=7
|G⁺|=235,200.  No-Frobenius half is signed PSL (`2|PSL|`).  Möbius
recovered from the permutation; class parameter `τ = tr²/det ∈ F_q`.
`Γ` on GPU (CuPy batched gather).  No flag flipped.

## What is theorem, or holds at both primes

Write `d=n/2`, `n=q+1`, `q=p²`.  Restrict to PSL (no Frobenius).  Then:

1. **Elliptic (anisotropic, 0 fixed points):** `Γ ≡ 0` and `χ_{V_+} ≡ 0`.
   All 7,200 elements at p=5; all 56,448 at p=7.  This is the
   character-table prediction: Z contains only `W_e` and principal series
   (15.589 B); principal series vanish off Borel classes; the even Weil
   of `PSL(2,p²)` is 0 on anisotropic classes (integer-valued Weil,
   `√q=p`).  So `Γ = ∑ λ_c χ_c` vanishes there automatically and **does
   not constrain any `λ_c`**.

2. **Identity:** `Γ(e) = n(n-2)` (already proved; `tr Φ`).  The extra
   signed element with `π=id`, `d≡-1` has the same `Γ` (`Γ(-g)=Γ(g)`).

3. **Involution (`τ=0`, ratio `r=-1`):** `Γ = 2(n-2)`.  Exact at p=5
   (`48`) and p=7 (`96`).  Not promoted to a general-p theorem here.

4. **Unipotents (parabolic, 1 fixed point):** exactly two `Γ` values, and
   `χ_{V_+}` takes the Weil unipotent values `(1±p)/2`.  Classical
   principal-series characters are `1` on non-identity unipotents, so

       Γ(u_±) = λ_exc · (1±p)/2 + (tr Φ − λ_exc d)/n.

   This is an identity, not a bound: it **reproduces** the two measured
   unipotent `Γ` at p=5 (`752/13`, `-128/13`) and p=7 (`34752/409`,
   `4512/409`) from the known `λ_exc`.  Unipotents do not separate the
   principal scalars.

5. **Split:** `τ` determines `Γ` (one value per `τ`).  Number of split
   `τ` is `(q-1)/4` (6 at p=5, 12 at p=7).  Further collapse: 6→3 values
   at p=5, 12→4 at p=7 (including the involution).  `χ_{V_+}=±1`.

## KILLED: an O(1)-in-p formula for all split `τ`

p=5-only fit, **predicted** for split no-Frob:

| class | predicted `Γ` | p=7 measured | |
|---|---|---|---|
| involution `τ=0` | `2(n-2)=96` | `96` | holds |
| elliptic | `0` | `0` | holds |
| `τ ∈ F_p^*` | `-4(n-2)/n = -3.84` | `-2208/409` and `-6240/409` (two values) | **fail** |
| `τ ∉ F_p` | `-4(n-2)p/n = -26.88` | `-3360/409` (one value) | **fail** |

p=5 had only **one** non-involution subfield class (`F_5^*/{r∼r^{-1}, r≠±1}`
is a singleton), so “`τ ∈ F_p` is one number” was an artifact.  Do not
reopen that formula.  Non-subfield split **did** collapse to a single
`Γ` at both primes; the value is not `-4(n-2)p/n`.

## Where leftover 1 actually sits

Identity, unipotent, and elliptic are either tautological or equivalent
to `λ_exc` plus the trace.  Every remaining Fourier coefficient of `Φ`
is on **split classes**:

    Γ(r) = λ_exc χ_W(r) + ∑_{α ∈ A_e} λ_α (α(r)+α(r^{-1})),

`A_e` the `(q-9)/8` principal-series parameters in 15.589 A, `r` the
eigenvalue ratio.  `λ_min(Φ)≥6` is this finite system together with
`λ_exc≥6`.  Next named step: identify `A_e` from the square map on the
Weil character (Adams 6.4), invert the split Fourier transform at p=5
(2 parameters) and p=7 (5 parameters), then **predict** the p=11
principal spectrum before claiming.

No flag flipped.  `leftover1` stays `global_qvar AND r1_l2`, both False.

---

# Step 4: A_e is the 4|k principal series; Fourier inversion is exact

`scripts/gamma_ae_fourier.py`.  Canonical PSL lifts (`d[0]=+1`, no
Frobenius).  `χ_Z = (χ_W(g)² + χ_W(g²))/2 − #fix(π)` (15.589 A).
PSL principal series `ρ(α_k)`: even `k ∈ (0,(q−1)/2)`, excluding the
quadratic `k=(q−1)/2`.  Inner products over `|PSL|`.  p=7 Γ on the V100.

## Identification of A_e

⟨χ_Z, ρ_k⟩ ∈ {0,1} exactly.  The 1's are:

| p | PS indices | A_e (`⟨χ_Z,ρ⟩=1`) | 4\|k in PS |
|---|---|---|---|
| 5 | 2,4,6,8,10 | **4, 8** | 4, 8 |
| 7 | 2,4,…,22 | **4, 8, 12, 16, 20** | 4, 8, 12, 16, 20 |

**Count identity (all odd p, q=p²≡1 mod 8).**  Even `k` in `(0,(q−1)/2)`
with `4|k` and `k ≠ (q−1)/2`: `(q−1)/2` is itself `0 mod 4` and is the
excluded quadratic, so there are `(q−1)/8 − 1 = (q−9)/8` such `k`.
This is exactly `|A_e|`.  Equivalently: `α` is trivial on the unique
subgroup `μ_4 ⊂ F_q^*`.

Do not promote the inner-product identification to a general-p character
theorem yet (verified at two primes + the count).  Do not reopen "A_e is
α trivial on F_p^*": that is 3 series at p=7, not 5.

## Fourier inversion (exact at both primes)

`χ_{W_e⊂Z} = χ_Z − ∑_{k∈A_e} ρ_k`  (not `tr(U|_{V_+})`; that pairing
missed `λ_exc`).  Then `λ_c = ⟨Γ, χ_c⟩`:

| p | `λ_exc=⟨Γ,χ_W⟩` | principal `λ_k=⟨Γ,ρ_k⟩` |
|---|---|---|
| 5 | `176/13` | `k=4: 80/13`, `k=8: 144/13` |
| 7 | `4320/409` | `4: 3360/409`, `8: 4032/409`, `12: 3648/409`, `16: 3072/409`, `20: 3360/409` |

Reconstruction `Γ = λ_exc χ_W + ∑ λ_k ρ_k` has **0 mismatches** on all
7800 PSL elements at p=5 and all 58800 at p=7 (every family).

p=7: `λ_4=λ_20` (multiplicity `2n` in the Φ spectrum).  `λ_8 ≠ λ_16`.
The pairing `k ↔ (q−1)/2−k` is `ρ(α)` vs `ρ(χ_2 α)`, not automatic
equality of scalars.

## p=11 spectrum (stored Φ, not a Γ census)

`phiZ_p11.npy`, dim Z=1769.  Clusters:

- `λ_exc=8.664378` (mult 61=`n/2`)
- 6 principal values, mult 122=`n`
- 4 principal values, mult 244=`2n` (coincident pairs)
- 6+4·2=14=`(121−9)/8` constituents.  All in `[8.054, 8.637]`.

Numerically leftover 1 holds at p=11 with margin (already in 15.593 G).
Not a p-law.

## Where leftover 1 sits, sharpened

```
λ_min(Φ)≥6  ⇔  λ_exc≥6  and  λ_k≥6 for every k=4,8,…,(q−1)/2−4.
```

The binding principal scalar is at p=5: `80/13≈6.154`.  At p=7 the
minimum is `3072/409≈7.511`; at p=11, `8.054`.  A closed form
`λ(k,q)≥6` would close the principal floor; QVAR remains `λ_exc≥6`.
p=7 cosine fit `λ(k)=(3456 − 192 cos 2θ − 480 cos 3θ)/409`,
`θ=2πk/(q−1)`, is **one prime** — do not predict p=11 from it.

No flag flipped.

---

# Step 5: leftover 1 is QVAR for every even character of F_q^*

`scripts/lambda_as_char_moment.py`.  Same constant as 15.589 E, now on
the whole of F=Z^U.

On F_q write z=y|_{F_q}, D={z=−1}, N(a)=|D∩(D−a)|,
Z_α=∑_{a≠0} α(a) N(a).  Even characters: α(−1)=1, i.e. even k.

**Certified p=5 and p=7, exact:** for every even k∈(0,(q−1)/2)
(all even nontrivial non-quadratic α, up to inverse),

    32 E|Z_α|² / [q(q−1)]

equals an eigenvalue of Φ, and these values are **exactly** spec(Φ).
The QVAR character is the unique even α with α²=χ (k=(q−1)/4); it
gives λ_exc (the **top** at p=5,7,11).  The other even k give the
principal scalars.  Count: (q−5)/4 = dim F (15.278).

Same factor for every such α because the Aut_∞-circulants β(t)=α(t)
all have ‖β‖²=q−1; 15.589 E computes the constant once.

Quadratic α=χ gives 2(n−2), not a Φ-eigenvalue (B_χ∉Z).  For
nontrivial even α, R(a)=∑ z_x z_{x+a} satisfies Z_R=−4 Z_N.

Leftover 1 is therefore the single inequality

    E|Z_α|² ≥ 3q(q−1)/16

for every even character α of F_q^* other than 1 and χ.  QVAR is one
character in that list, not the binding one (p=5 binds at k=2: 80/13;
p=7 at k=8: 3072/409).  That is 15.277 G, named.

Not proved.  Do not replace leftover 1 by QVAR alone.  The 6+4/n
trace is false (mixes Wick-6 with λ̄; “only k=4 gets 4-point” dies at
p=7).

Parseval (character orthogonality, Max+-free): even α satisfy
∑_{α(−1)=1} |Z_α|² = (q−1) ∑_{a≠0} N(a)².
At p=5 this checks: E∑ N²=1200, (q−1)∑N²=28800, and the even Fourier
mass (trivial 26100 + χ 900 + the rest 1800) matches.  Threshold
3q(q−1)/16=112.5; min E|Z_α|²=1500/13≈115.38; gap 75/26.  Mean of
the leftover characters is 180, so the average clears and the min
just clears.  Average > threshold is not a proof of the min
(equal-density shortcut, already killed).

p=11 even-character scan (37.4M, 86-way fork) with **15590 encoding
is not a Step 5 kill**: that array is labeled by
`paley_conference_prime_power` (`e=a+bp`, ω²=2).  HIP rerun on nuka
9070 XT with that field (`even_char_hip_nuka.py minmax`): **29/29**
even k match spec(Φ) to ~1e-6.  GPU 9.75s / wall 12.5s.

| p=11 | value |
|---|---|
| QVAR k=30=(q−1)/4 | 8.664378 = λ_exc (mult 61) |
| binding min k=8,28,32,52 | 8.054447 = Φ_min (mult 244) |
| threshold 3q(q-1)/16 | 2722.5; min E\|Z_α\|²=3654.7 |

Step 5 is now a three-prime identity (p=5,7 MuLab; p=11 minmax
labeling).  Leftover 1 is still the **min** over those even α, not
the mean and not QVAR alone.  No flag flipped.

This is the same object as **15.279 Φ|_F** (even characters of
F_q^×/{±1}, pairs {ψ,χψ}, dim F=(q−5)/4).  Identified constants:
λ=2 E|M_ψ|²/(p⁴(p²−1))=32 E|Z_α|²/[q(q−1)].  Theorem D:
λ=8+ˆR_rest/q², so leftover 1 is ˆR_rest≥−2q² for every even
ψ∉{1,χ}, equivalently Gauss 4-distinct pairing of m₄ ≥0 (15.279 L),
equivalently F̂(ψ)≥0 on squares (15.279 M).  That is the QVAR goal
language, now required for **all** even characters, not just α²=χ.

Measured ˆR_rest at the binding character (budget −2q²):

| p | λ_min | ˆR_rest | −2q² | margin |
|---|---|---|---|---|
| 5 | 80/13 | −15000/13 | −1250 | 1250/13=(λ−6)q² |
| 7 | 3072/409 | −1174.08 | −4802 | 3628 |
| 11 | 8.054447 | +791 | −29282 | λ_min>8 so remainder **positive** |

|ẑ|² on Ω (15.279 support): energy off Ω is 0; E[|ẑ|²]=2q.  Not
two-valued (p=5: five values in Q(√5), including 0; p=7: many
values, including 0).  Pointwise |M|²=0 is attained; ensemble
average is the floor.  No flag flipped.

---

# Step 6: Aut-orbit values of L; two p-laws killed at p=11

`scripts/aut_orbit_L_and_lambda_fit.py`.  L(r)=E ∑_{δ≠0} N(δ)N(rδ)
on squares.  Aut-invariant (15.279 Q), so constant on ⟨Frob,inv⟩-orbits
of T=F_q^×/{±1}.  Full Max+ mix (both y_∞) at p=5,7; y_∞=+1 half
matches Wick+Plancherel for L(1) and is the p=11 `eps1` convention.

## Named orbit values (full mix)

p=5, 4 orbits (formula (p+3)²/16=4):

| rep | size | order | Paley (χ(r−1),χ(r+1)) | L |
|---|---|---|---|---|
| 1 | 1 | 1 | (0,+) | 1200 |
| i ∈ F_p | 1 | 4 | (++) | 14550/13 |
| cube | 2 | 3 | (−+) | 14550/13 |
| binding | 2 | 12 | (−−) | 14250/13 |

p=7, 6 orbits:

| order | Paley | L |
|---|---|---|
| 1 | (0,+) | 8232 |
| 3 and 4 (i and ω₃) | (++) | 3239880/409 |
| 12 (−−) and 24 (+−) | | 3236352/409 |
| 8 (−−) binding | | 3227532/409 |

y_∞=+1 only: L(1)=450 (p=5), 4116 (p=7), 76230 (p=11 sample).
These three match Wick Q(1)=8q² + Plancherel + A=p(2−p)+4N, which
is **y_∞=+1 only**.  y_∞=−1 uses A=−p(p+2)+4N (nD=p(p+1)/2).
Binding orbit stays the unique smallest leftover at both primes.

## KILLED: L(i)=L(ω₃) as a p-law

At p=5 and p=7, the 4th-root and 3rd-root orbits have **the same** L
(full mix and y_∞=+1 half).  That is **not** Aut-invariance: they are
distinct ⟨Frob,inv⟩-orbits.  On every Type+ 1D lift the difference is
a nonzero constant (p=5: +100 on all 3 lifts; p=7: −490 on all 4),
so the equality is an ensemble-mixing identity, not a 1D identity.

p=11, 200k sample of `maxplus_p11_eps1.npy` (y_∞=+1):

| r | L mean | se |
|---|---|---|
| 1 | 76230.22 | 1.67 |
| i | 74346.90 | 1.02 |
| ω₃ | 74365.84 | 0.76 |

L(i) ≠ L(ω₃) at ~14σ on this sample.  Do not promote the p=5,7
coincidence.  Leftover Aut-dofs stay n_orb−2.

## KILLED: p-independent cosine model of λ(k)

Fit λ(k)−8 = a₀ + ∑_{m=1}^M a_m cos(2π m k/(q−1)) on p=5+7 jointly,
predict p=11 (HIP 29 even k).  Max |λ−pred| ≳ 0.95 already at M=1
and does not drop below ~2 at M=6.  Per-prime: p=5 is exact at
M=2, including QVAR,

    λ(k) = 8 + 8/13 − (64/13) cos(π k / 6)     (even k, p=5 only).

p=7 is not a cosine in θ=2πk/(q−1) of degree ≤7 across all even k
(the PLAN p=7 formula is PSL indexing, principal series only, not
QVAR).  Polynomial in x=cos(4πk/(q−1)) is exact at p=5 (deg 1) and
p=7 (deg 5) and still 0.35 maxerr at p=11 deg 1.  Do not predict
p=11 from a p=5/p=7 cosine.

Naive bound L≥L_min on leftover squares gives
E|Z|² ≥ 2(L(1)−L_min)=2700/13≈207.7 at p=5, above threshold 112.5,
but the true min is 1500/13≈115.4: the binding character's oscillation
against the leftover δ=L−L_min eats the slack.  Uniform L_min is not
a proof.

No flag flipped.  Floor is still F̂(ψ)≥0 for every even ψ∉{1,χ}.
Next constraint that could cut Aut-dofs is the Boolean cubic on Ω
(15.279 T) as a linear relation among leftover orbit values, not
another Paley-type or cyclotomic collapse.

---

# Step 7: Boolean cubic on Ω does not cut leftover Aut-dofs of Q

`scripts/boolean_cubic_orbit_relations.py`.  15.279 T, all Max+ at
p=5,7 (MuLab; vectorized DFT).

The cubic is the Fourier of \(z_x^2=1\):

    ∑_η ẑ(η)ẑ(ξ−η) = q² 1_{ξ=0}.

Restricted to ξ∈Ω with ẑ supported on {0}∪Ω this is

    2 ẑ(0) ẑ(ξ) + ∑_{t∈R} ẑ(tξ)ẑ((1−t)ξ) = 0,

ẑ(0)=p y_∞.  Residual of `2ẑ(0)ẑ+∑B` is 2e-13 (p=5) / 2e-12 (p=7).
The published form `2pẑ+∑B` holds only on y_∞=+1 (fails on y_∞=−1
by |4pẑ|).  Magnitude |∑B|=2p|ẑ| on both signs.

Squared and averaged: ∑_{t,s∈R} Γ_{t,s} = 4p² E u = Q(1)=8q², and
`4p² Q(r)=E[|∑B|² u(r·)]` to 1e-10.  Both are the cubic tautology
(15.279 T: “Boolean rewrite returns M”).  They do **not** give a
linear relation on leftover Aut-orbit Q beyond Q(±1)=8q² and the
row-sum ∑_{r□} Q(r)=2q²(q−1) (checked exact).

Diagonal of the bilinear Gram is Q: Γ_{t,t}=Q(ρ(t)), ρ(t)=(1−t)/t.
ρ(R) **misses** leftover orbits (p=5: the i-orbit is not in ρ(R)).
Naive |ẑ|²-Wick of Γ is 5× too big (∑=25000 vs 5000 at p=5).

Off-diagonal Γ is **not** linear in Aut-orbit Q (lstsq maxerr 757
at p=5, 1226 at p=7).  At p=7 it is not even a class function of
(orb ρt, ρs, s/t, t=s, t=1−s): 10 split buckets, spreads ~10³.
(At p=5 those keys do not split, still not linear in Q.)

Γ ⪰ 0 is automatic (Gram of the B_t).  Putting a Q-model on the
off-diagonal and requiring PSD is the SOS-4 / linear 4-point
relaxation already killed.

**Killed:** Boolean cubic on Ω as a source of new linear constraints
on leftover Aut-orbit Q.  Aut leftover dofs stay n_orb−2.  Floor
F̂(ψ)≥0 still OPEN.  No flag flipped.

Named Q on Aut-orbits (full mix; Wick off ±1 is 4q²):

| p=5 orbit | Q | Q−4q² |
|---|---|---|
| ±1 | 5000 | +2500 |
| i and ω₃ | 30000/13 | −2500/13 |
| binding (order 12) | 20000/13 | −12500/13 |

p=7 binding (order 8) is the unique min leftover Q≈8077.7 vs Wick 9604.

Next that is not a rewrite of z²=1 or Aut-invariance: Boolean 4-point
of V_+ as a function of Aut-orbits of 4-tuples (cross-ratio, not a
linear model in {κ,star,φ} — that already fails at p=5,7).

---

# Step 8: Boolean 4-point of V_+ is Aut-constant and not a function of (κ,CR,star)

`scripts/m4_aut_orbit_vplus.py`.  Max+ ⊂ V_+ ∩ {±1}^n, p=5 and p=7
exact (MuLab).  4-sets of P¹, permutation Aut(C).

m₄(S)=E[∏_{i∈S} y_i] is constant on every Aut-orbit: 42/42 at p=5,
128/128 at p=7 (dead_labels=0).  Distinct finite values:

- p=5 (denom 65=N/4): ±{1,3,9,13,21}/65  (ten values)
- p=7 (denom 2863=N/4=7·409): ±{45,61,65,67,99,109,139,191,279,327}/2863
  (twenty values)

Linear {1,κ,φ,star} maxerr 0.114 (p=5) / 0.027 (p=7).  15.597
particular solution vs true m₄: same size error.  Confirmed: true
Boolean 4-point is **not** the particular m₄.

(κ, has_∞) splits all 8 keys.  (κ, has_∞, CR) still splits 6/21
(p=5) and 27/44 (p=7).  Adding star and φ does not kill the split
(4 leftover at p=5, 24 at p=7).  The scratch line “affine m₄ is a
function of (cross-ratio, κ_C, star)” is **false** as a complete
invariant: Aut-orbits are strictly finer.  (CR labels used
`canon_cr_fn` minmax encoding; Aut-orbits themselves are 15590/C
and do not depend on that.)

⟨m₄, κ_A⟩ lives entirely on **finite** 4-sets (∞-containing orbits
contribute 0 for both QVAR and the binding even character):

| p | ⟨m₄,κ_QVAR⟩ | ⟨m₄,κ_bind⟩ |
|---|---|---|
| 5 | 14.135 | 0.180 |
| 7 | 41.917 | 6.941 |

Both positive (census floor).  Binding at p=5 is tight.  No
single-orbit closed form in p presented itself (values are odd
over N/4, not a Jacobi list we could match).  Character-sum of m₄
on Aut-orbits of 4-sets is the same open as 15.48 for g_min.

**Killed:** m₄ as a function of {κ,φ,star}; m₄ as a function of
(κ,CR,star).  Aut-constancy is Aut-invariance, not a floor.
Do not add an identity file.  Floor F̂(ψ)≥0 still OPEN.  No flag
flipped.

---

# Step 9: q-dependent split-Γ formula — cosine and Jacobi fail the p=11 gate

`scripts/split_gamma_dilation_ansatz.py`.  Square dilation z↦t z
on P¹ (∞↦∞) is the gauge of a split PSL element; s=1 Aut lift iff
t is a square.  Γ(t)=E[(y·πy)²]−2n on that lift.

## What holds

Identity dilation t=1: Γ=n(n−2) (tr Φ).  Involution t=−1:
Γ=2(n−2) **exactly** at p=5 (48) and p=7 (96).  p=11 150k sample
237.3 vs 240 (sampling).  **Not pointwise:** y·Uy takes many values
(p=5: {−14,−2,2,6,10,14,26}, not only ±2p).  Ensemble only.  Do not
promote to a general-p unit without a Max+-free argument; 4-point
of this one involution is a special pairing, not the floor.

Split Aut Γ takes **2** leftover values at p=5 (−240/13, −48/13)
and **3** at p=7 (−6240/409, −3360/409, −2208/409).  These are the
PLAN Step 3 class-function values.  Binding is the unique min
(p=5 order-12 (−−); p=7 order-8 (−−)).  At p=5 and p=7, Γ is a
function of (χ(t),χ(t−1),χ(t+1),in F_p, order) with 0 splits —
a two-prime artifact (cf. “τ∈F_p is one number” at p=5).

## KILLED at the p=11 gate

1. **λ−8 = a₀(q)+a₁(q)cos(4πk/(q−1))+…** with a_j=A+B/q, n_harm≤3,
   fitted on all even k at p=5,7.  p=11 max|λ−pred|≈3.7–4.0 and
   predicted min 4.8–5.8, **below 6**.  Do not reopen p-independent
   or A+B/q cosine.

2. **λ−8 linear in Re J(χ,α_k)/q, |J|/q, 1/q.** |J|=p constantly
   (Weil).  Fit on p=5,7 already maxerr 2.87.  Phase of J does not
   linearly give λ (p=5 QVAR k=6 has Re J=2.24 but the largest λ−8;
   k=4 has Re J=p and medium λ).

3. **(p−1)/2 unique split Γ values** (2,3 at p=5,7).  p=11 sample
   has many distinct split Γ (23 at 0.1 rounding; involution se is
   ~3 so this oversplits, but it is not 5).  Paley type of t±1 plus
   order does **not** determine Γ at p=11 (order 60 appears in many
   buckets).  The 0-split 5-tuple at p=5,7 is the same small-prime
   collapse as subfield-one-value.

Dilation-Γ and even-character λ are Mellin duals.  A short
q-dependent formula on one side is a short formula on the other.
Bounded-term cosine/Jacobi on λ failed the gate; Aut leftover
dofs of Γ(t) remain.  Floor F̂(ψ)≥0 still OPEN.  No flag flipped.
No identity file.

---

# Step 10: involution 4-distinct mass is 2(n−2); Kloosterman is not Γ(t)

`scripts/involution_kloosterman.py`.

Involution dilation t=−1, s=∑ d_i y_i y_{πi}, π(x)=−x.  Index split of
E[s²] (empirical, p=5 and p=7, exact):

| |{i,πi,j,πj}| | #(i,j) p=5 | mass | 2-point pred |
|---|---|---|---|
| 1 (fixed pts) | 2 | 2 | 2 |
| 2 | 50=2q | 2q | 2q |
| 3 | 96 | 0 | 0 |
| 4 (two 2-cycles) | 528 | **2(n−2)** | unknown |

So Γ(−1)=2(n−2) iff the 4-distinct mass is 2(n−2).  2-point/Boolean
on collisions already match and cannot supply it.

Those 4-sets are involution rectangles {a,−a,b,−b}.  m₄ is **not**
constant on them (p=5: six values, all the positive Aut-orbit m₄'s;
p=7: ten).  The mean is **2/(q−3)** = 1/11 at p=5, 1/23 at p=7.
Boolean on the cycle monomials m_P=z_a z_{−a}=±1 gives this mean
iff E[(∑_P m_P)²]=q−1.  That is equivalent to the involution identity,
not a 2-point proof.  (DFT: ∑_x z_x z_{−x}=(1/q)∑_ξ ẑ(ξ)², square not
|ẑ|²; the |ẑ|² Plancherel would have forced the pairing constant,
which is false.)

Kloosterman Kl(1,t) and χ-Bessel ∑ χ(x) e(Tr(t x+x^{−1})) do not
linearly give Γ(t) (maxerr hundreds including identity; on leftover,
|Kl| takes more values than Γ and is not monotone with the binding
orbit at p=7).  Not a formula.

No identity file.  Floor still OPEN.  No flag flipped.
