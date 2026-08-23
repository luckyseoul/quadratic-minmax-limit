# Fable query — is W2 a p-law? (open)

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
live at p=11 (m=15; orbits {Φ_3}, {Φ_5}, {the two Φ_15-quartics}). W1 is a
separate problem (tracked elsewhere) — this query is about W2 only.

## The question

**Is W2 true for every odd prime p?** If so, what's the right proof
mechanism? If not — or if it's not clear — what does the failure or
uncertainty look like?

We are NOT asking you to evaluate one particular proof strategy. We want
your own read on the right way to attack this: a counting argument, a
structural/algebraic argument about the ⟨I,Frob⟩-orbit action, an inductive
or recursive construction of a witness, a reduction to something already
known, a reason to doubt it's a p-law at all, or something we haven't
thought of.

## Data available, for context (not a hint toward any particular method)

W2 is certified (by direct computation) at p = 5, 7, 11, 17, 31.  At p=11
the orbits are {Φ_3}, {Φ_5}, {the two Φ_15-quartics} and W2 holds for all of
them.  One thing we noticed while exploring: the "switched split-involution"
class `{±[[α,β],[γ,−α]] : α²+βγ=1}` ⊂ PGL(2,q) (all Max− under switching)
has an increasing hit rate for W2-satisfying elements — 17/153 ≈ 11% at
p=17, 76/496 ≈ 15% at p=31 — but we do not know whether this is the right
object to look at, whether the trend continues, or whether it's even
necessary that a witness come from this particular class. Take it or leave
it.

## Already tried and failed — so you don't repeat them, not to redirect you

* **Named W2 witnesses as a uniform p-law.** t = −2 works at p=17, fails at
  p=31. x/(x−1) hits at p=31 but isn't shown to generalize.
  PGL(2,q)·z is Φ_3-dead; the χ_p-pullback misses Φ_3.
* **A Krylov/gcd shortcut** claiming z+Dz ∈ ker g — this is FALSE;
  f(D)(z+Dz) ≠ 0 for every irreducible factor of g at p=5,7,11. (The
  correct membership test is: w ∈ (f)R iff f divides the γ-content of w,
  not iff f(D)w = 0.)

## Ground rules

Answer must be precise enough to eventually code as a `True`/`False`
predicate, or an honest statement that you don't see how to close it and
why. No soft-closing — if an argument only covers large p or a subfamily,
say so explicitly and give the boundary. Any numerical pattern must be
checked out of the fitting range before being trusted — we've been burned
once already by a hypothesis that held for 3 primes and broke on the 4th.
