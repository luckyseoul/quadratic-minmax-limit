# Fable query — W2 as a p-law (counting identity on the split-involution conic)

## Setup (self-contained)

p odd prime, q = p², N = q−1 = 2^a·m with a = v_2(N) ≥ 2 and m odd.
Max− = the (−p)-eigenvectors of the Paley conference matrix of order n = q+1
on P¹(F_q).  U = the {∞,0} pair-slice of Max−.  dir(U) = F2-span of the
pairwise differences of U.

W ≅ F2[X]/(X^N+1);  W_0 ≅ R = F2[X]/h,  h = (X^N+1)/(X+1) = (X+1)^{2^a−1}·g^{2^a},
g = (X^m+1)/(X+1).  D = the multiplier action of a generator of F_q^*,
I(z) = 1/z, Frob(z) = z^p; these generate Aut({0,∞}) and preserve U, so
I_U := dir(U) ∩ W_0 is an Aut-invariant ideal of R.

**15.612 (PROVED):** the maximal proper Aut-invariant ideals of R are exactly
(X+1)R and (f_O)R for each ⟨I,Frob⟩-orbit O of irreducible factors of g,
f_O = ∏_{f∈O} f.  Hence Walsh (15.406 E, ⇒ leftover 2) ⟺ **W1 ∧ W2**:

* **W1:** I_U ⊄ (X+1)W_0 — some U-difference has (X+1)-valuation 0.
* **W2:** for every orbit O, some U-difference is nonzero mod f_O.

W2 is vacuous at p=3, implied by W1 at p=5,7 (g = Φ_3 irreducible), and first
live at p=11 (m=15; orbits {Φ_3}, {Φ_5}, {the two Φ_15-quartics}).

## The question

**Prove W2 for every odd prime p** (currently certified only at p = 5, 11, 17, 31).

## The specific attack we want evaluated

Take the **switched split-involution class**
`C_p = {±[[α,β],[γ,−α]] : α²+βγ = 1}` ⊂ PGL(2,q) — a conic, |C_p| = p(p+1)/2 —
all of whose elements are Max− under switching.  Census:

| p | \|C_p\| | in U | W2 hits | rate |
|---|---|---|---|---|
| 17 | 153 | 49 | 17 | 11.1% |
| 31 | 496 | 146 | 76 | 15.3% |

The hit rate is a growing positive fraction, so the natural completion is a
**counting identity** rather than a clever named witness (15.627 B states the
gap in exactly this form: "some conjugate is always coprime").

Concretely, we want to know:

1. Can the W2 condition — "the U-difference attached to π ∈ C_p is nonzero
   mod f_O" — be written as the **non-vanishing of a polynomial of bounded
   degree** in the conic coordinates (α,β,γ), with degree bounded
   independently of p (or by deg f_O, itself bounded by the orbit structure
   of ⟨I,Frob⟩ on the factors of g)?
2. If so, does Weil / Lang–Weil on the conic give
   `#{π ∈ C_p : W2 holds} = c·p² + O(p^{3/2})` with `c > 0` explicit,
   hence W2 ≠ ∅ for all p above an explicit bound (small p by census)?
3. If the degree is **not** p-bounded, is there a substitute — e.g. an
   averaging/second-moment argument over C_p, or a Chebotarev/Frobenius
   equidistribution statement over the factors of g — that still forces a
   nonzero residue for at least one π ∈ C_p?
4. Independent of the above: is there a **structural** reason some conjugate
   must be coprime to f_O — e.g. the ⟨I,Frob⟩-orbit structure on the
   factors of g forcing a non-degenerate pairing?

## Already dead — do not propose these

* **A single named W2 witness as a p-law.** t = −2 works at p=17 and fails at
  p=31 (15.626).  x/(x−1) hits at p=31 (15.627 B).  PGL(2,q)·z is Φ_3-dead
  and the χ_p-pullback misses Φ_3 (15.620).  The named-witness route has
  repeatedly failed out of sample.
* **Krylov / gcd shortcuts:** the claim that z+Dz ∈ ker g is false;
  f(D)(z+Dz) ≠ 0 for every irreducible factor of g at p = 5,7,11 (15.616,
  corrected in 15.617 — membership is `w ∈ (f)R` iff f divides the
  γ-content, *not* iff f(D)w = 0).
* **W1 side-quests.** W1 is proved for p ≡ 5 (mod 8), p ≡ 17 (mod 24),
  p ≡ 73, 97 (mod 120), and all (2/p)₄ = −1.  Its residual class
  (2/p)₄ = +1 has both the **linear box** in (a,b,i,k) killed over 61 primes
  (15.627 A) and the **entire interval family** killed
  (2026-08-23: ε(−(p−1)/4) ≡ 0 and ε(−(p−1)/8) ≡ 0 on the whole residual
  class; ε(−(p−1)/16) exists only for p ≡ 1 (mod 16) and is mixed there).
  W1 is a separate problem — please stay on W2.

## Ground rules

Answer must be a proof sketch precise enough to code as a `True`/`False`
predicate, or an explicit statement that the approach fails and why.  No
soft-closing: if the argument only covers large p, say so and give the bound.
Any numerical criterion must be fitted below p = 2000 and then **predict**
p ∈ (2000, 10⁴) — a pre-asymptotic fit on ~8 primes has already burned this
project once.
