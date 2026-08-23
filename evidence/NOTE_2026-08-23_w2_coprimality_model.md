# W2: unit-content heuristic and the corrected factor-orbit target

Date: 2026-08-23.  Target: W2 (15.612), the odd-factor half of Walsh, which
with W1 closes leftover 2's 15.406 E.  **No flag flipped; W2 stays open.**

## Reframing (structural, from the PIR)

R = F2[X]/h is a principal ideal ring, so I_U = (γ) for a unique monic
γ | h, Aut-invariant: γ = (X+1)^i·∏_O f_O^{e_O}.  Walsh ⟺ γ=1;
W1 ⟺ i=0; W2 ⟺ all e_O=0; deg γ = codim I_U.

For an individual content c, the bad event for the maximal Aut-invariant
ideal `(f_O)` is

    c ∈ (f_O)  ⟺  f_O | c
                 ⟺  every irreducible factor in O divides c.

Thus a content coprime to g is a sufficient simultaneous witness for all
orbits, but it is not necessary.  W2 only asks, separately for each O, for
some c not divisible by the whole product f_O; the witness may depend on O.
In particular `gcd(c,f_O)=1` is too strong whenever O has more than one
irreducible factor.

## The model and the census

Random-element heuristic: Pr[content coprime to g] = ∏_{f|g}(1−2^{−deg f}).
Census of the switched split-involution class (rates conditioned on in-U —
the class-level 11%/15% figures conflate three filters):

| p | m | model | measured in-U rate | direction predicted/observed |
|---|---|---|---|---|
| 17 | 9 | 0.738 | 17/49 = 0.347 | — |
| 31 | 15 | 0.618 | 76/146 = 0.521 | ↑ / ↑ |
| 41 | 105 | 0.458 | 63/237 = **0.266** | **↓ / ↓** |
| 47 | 69 | 0.749 | 205/318 = **0.645** | **↑ / ↑** |

p=41 and p=47 were run AFTER the model was fixed (this session,
`scripts/w2_class_scan_general.py`, 60-way ProcessPool over `_switched`).
The measured rate zigzags 0.35 → 0.52 → 0.27 → 0.65 and the model called
BOTH reversals in advance: down at p=41 because 7 | m brings two cubic
factors, up at p=47 because m=69 has only three orbit-factors.  A naive
trend-reading would have been wrong twice.  Measured/model ratio stays in
[0.47, 0.86]: a persistent deficit (differences are not perfectly random)
but bounded away from collapse, with the variation structure entirely the
factorization of m.  Also uniform across all four primes: EVERY class
element is Max− under switching (n_eigen = n_class), and the in-U fraction
drifts slowly (0.320, 0.294, 0.275, 0.282).

## The strong unit-content gap (historical route)

w_t mod f is NOT low-degree in the conic parameter t — the difference bits
carry χ(γj+δ) switches, so the residue is a mixed character sum and the
"few roots" closure does not apply.  What remains is a second-moment
estimate: for each irreducible f | g,

    |B_f| := #{t : f | w_t}  ≲  2^{−deg f}·|T|·(1+o(1)),

via Σ_t |w_t(ζ_f)|² (Weil-class input).  A union bound could leave a
unit-content witness, but that would prove a strictly stronger statement
than W2.  The corrected per-orbit target below makes this route unnecessary.

## Superseded next steps from the unit-content route

1. γ(p) exactly at p ≤ 41 (one gcd chain over a spanning set of
   U-differences).  If γ = (X+1)^c uniformly, W2 is subsumed by the named
   conjecture γ-form.
2. Per-factor breakdown of the p=41 misses (which f kills each) against
   the per-factor model 2^{−deg f}; mismatch localizes the structure a
   proof must use.
3. The second-moment estimate itself.

---

## Logical correction: the old union-bound "kill" was an artefact

The earlier version of this note reversed the ideal-membership condition:
it called an orbit bad when the content shared *some* irreducible factor
with f_O.  That is the event `gcd(c,f_O) != 1`, not the W2 event `f_O | c`.
Consequently its per-orbit table, summed bad rates, and claimed
union-bound obstruction did not measure W2 and are retracted.  Those data
remain relevant only to the stronger unit-content problem.

No union across O is required at all.  For each O independently it is
enough to prove

    B_O = {t ∈ T : f_O | c_t}  is a proper subset of T.

This is a much weaker and cleaner target.

## Complete normalized pole family

Take the split involutions

    π_t(x) = x/(t x - 1),   t ∈ F_p^×,

with the standard switching used by 15.622.  This is not the earlier
one-parameter family `x/(x-τ)` dismissed in `HANDOFF.md`: after projective
normalization its varying lower-left entry is essential.

For the named halfspace point z, y_0=+1 automatically.  If λ is the
ω-coordinate of σ^{-1}, then y_∞=-1 exactly when

    u = λ/t ∈ {(p+1)/2, ..., p-1}.

Hence the admissible set T has exactly (p-1)/2 elements, with an explicit
upper-half-interval parameter u.  `scripts/w2_ramanujan_mask_spectrum.py`
factors g completely, constructs the factor orbits under X↦X^p and
X↦X^{-1}, and tests the correct product-divisibility event.

| p | |T| | unit c_t | one c_t clears every f_O | max_O |B_O|/|T| | #O |
|---:|---:|---:|---:|---:|---:|
| 5 | 2 | 1 | 1 | 1/2 | 1 |
| 7 | 3 | 2 | 2 | 1/3 | 1 |
| 11 | 5 | 2 | 3 | 2/5 | 3 |
| 13 | 6 | 4 | 4 | 2/6 | 3 |
| 17 | 8 | 5 | 5 | 3/8 | 2 |
| 19 | 9 | 5 | 7 | 1/9 | 5 |
| 23 | 11 | 11 | 11 | 0/11 | 3 |
| 29 | 14 | 4 | 7 | 7/14 | 7 |
| 31 | 15 | 8 | 8 | 5/15 | 3 |
| 37 | 18 | 11 | 11 | 6/18 | 7 |
| 41 | 20 | 9 | 9 | 9/20 | 8 |
| 43 | 21 | 10 | 12 | 9/21 | 7 |
| 47 | 23 | 18 | 18 | 5/23 | 3 |
| 53 | 26 | 12 | 12 | 13/26 | 11 |
| 59 | 29 | 13 | 13 | 14/29 | 9 |
| 61 | 30 | 17 | 17 | 8/30 | 15 |
| 67 | 33 | 23 | 23 | 8/33 | 14 |
| 71 | 35 | 15 | 19 | 8/35 | 15 |
| 73 | 36 | 23 | 23 | 10/36 | 7 |
| 79 | 39 | 24 | 30 | 8/39 | 9 |
| 83 | 41 | 18 | 25 | 15/41 | 13 |
| 89 | 44 | 17 | 17 | 19/44 | 11 |
| 97 | 48 | 24 | 24 | 18/48 | 5 |
| 101 | 50 | 15 | 25 | 21/50 | 25 |
| 103 | 51 | 39 | 39 | 11/51 | 15 |
| 107 | 53 | 37 | 37 | 16/53 | 7 |
| 109 | 54 | 18 | 29 | 16/54 | 15 |
| 113 | 56 | 18 | 27 | 27/56 | 15 |
| 127 | 63 | 29 | 43 | 17/63 | 7 |
| 137 | 68 | 36 | 36 | 31/68 | 13 |

Thus this family collectively certifies the W2 condition at every tested
prime p=5,...,137 (30 primes, excluding p=3 where W2 is vacuous).  The
stronger empirical statement

    |B_O| ≤ |T|/2

holds throughout this continuous range; equality occurs at p=29 and p=53.
It is not a theorem and in fact fails in the sparse large-p check below.

| p | |T| | unit c_t | one c_t clears every f_O | max_O |B_O|/|T| | #O |
|---:|---:|---:|---:|---:|---:|
| 191 | 95 | 57 | 67 | 21/95 | 7 |
| 193 | 96 | 65 | 65 | 30/96 | 5 |
| 223 | 111 | 66 | 84 | 25/111 | 11 |
| 257 | 128 | 57 | 57 | **65/128** | 10 |

At p=257 the Φ_3 singleton orbit has 65 bad parameters, killing the
half-bound.  Collective W2 remains true (63 Φ_3-good parameters, and every
other singleton factor has at most three bad parameters).  The viable
general target is therefore only `B_O != T`, exactly what W2 needs.

The p≤73 rows used the full cyclic coordinate solve.  The p≥79 extension
uses `scripts/w2_pole_fourier_fast.py`.  If ζ is a root of an irreducible
factor f and w=c(D)γ, then on either multiplicative point-orbit

    Σ_j ζ^j w(g^j) = c(ζ) Σ_j ζ^j γ(g^j).

For each f, cyclicity of γ makes the right generator projection nonzero on
at least one of the square/nonsquare point-orbits.  On that orbit the left
sum vanishes exactly when f|c.  This gives the same complete irreducible
mask without constructing or inverting the N×N γ-orbit coordinate matrix.
The fast and coordinate backends agree bit-for-bit at p=29 and p=89.

### Exact normalization from the Bose affine relative difference set

The norm identity used by the factor-free backend is a theorem, not an
empirical normalization.  Put `G=F_{p^2}^*`, `N=F_p^*`, and let `L` be the
nonzero F_p-linear functional whose kernel is the line `F_p b` used in the
definition of the named cyclic generator.  Then

    R = {x in G : L(x)=1}

is the classical affine `(p+1,p-1,p,1)` relative difference set in `G`
relative to `N`; see Dukes--Ling,
[Relative difference sets partitioned by cosets](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v24i3p64/pdf/),
Section 2.  Its integral group-ring identity is

    R R^(-1) = p e + (G-N).

Let `H=ker(L)\{0}=bN`.  On `G`, the named generator is `gamma=R+H` over
F2.  Also `RN=G-H`, so both cross-products `RH^(-1)` and `HR^(-1)` reduce
to `G+N`, while `HH^(-1)=(p-1)N=0`.  Reducing the displayed relative
difference-set identity modulo two therefore gives the exact identity

    gamma gamma^(-1) = e + G + N       in F2[G].

Project to the square subgroup `S=G^2`.  At every nonprincipal odd-order
character of `S`, the `S`-sum vanishes.  The `N`-sum also vanishes: either
the restricted character is nonprincipal, or it is principal and its sum
is `|N|=p-1=0` in F2.  Hence the square/nonsquare two-component Fourier
norm of gamma is exactly `1` at every root of g.  Consequently, for
`w=c(D)gamma`, the folded autocorrelation computed by
`scripts/w2_translated_antipodal_norm_scan.py` is exactly
`c(X)c(X^-1)` modulo g.  Coprimality of that norm with g is therefore a
factorization-free certificate that c is a unit modulo g.

### Antipodal poles: the formal identity and a valid translated family

For the original basepoint, pairing the formal pole parameters `u` and
`-u` produces unit content for every pair at every prime `5<=p<=137` in
the complete scan.  More strongly, its folded autocorrelation equals the
gamma norm exactly, so its content satisfies

    c(X)c(X^-1) = 1 mod g.

This does **not** prove W2: when `u` is in the normalized upper interval,
`-u` is not, so one endpoint is outside U.  Translating the halfspace
basepoint changes the content, and the tempting fixed translation `s=-1`
already fails at p=5,29,53.  Exact coordinate searches at
p=5,7,11,13,17,29 also show that the formal pair is not in the span of all
U-points obtained from the translated pole family.  Thus there is no
hidden bridge from the formal identity to I_U inside this family.

There is nevertheless a valid translated construction.  Write `s=-a` and
`m=(p-1)/2`.  Both antipodal pole endpoints belong to U exactly when

    1 <= a < u <= m,       a+u > m.

The xor of their two basepoint differences is then a genuine difference
of two U-points.  On the boundary edge `u=m`, the factorization-free scan
finds at least one unit-content pair for **every prime 7<=p<=401**.  Through
p=251 the first unit always occurs with `a<=9`.  In the range 257<=p<=401
the first unit has `a<=14` except at p=373, where the first unit is
`(a,u)=(29,186)`.  These small observed bounds are not claimed as p-laws.
At p=5, whose translated boundary edge is the sole failure, the original
normalized pole family already contains a direct unit-content U-difference.
Hence W2 is computationally certified at all 77 primes `5<=p<=401`, while
the uniform proof remains open.

The edge scan is substantially stronger than a factorized orbit census:
it verifies U-membership of both endpoints and proves `gcd(c,g)=1` from
the norm alone, without factoring g or solving the full cyclic coordinate
system.  The complete evidence is generated by
`scripts/w2_translated_antipodal_norm_scan.py`; the mesh sweep is evidence,
not a theorem.

Literature cross-check: the Gleason–Prange theorem and Ding–Liu–Tonchev's
classification of binary cyclic codes invariant under PSL(2,n)
([arXiv:1704.01199](https://arxiv.org/abs/1704.01199)) are nearby but do not
close this target.  Here I_U is forced invariant only by the {0,∞}-pair
stabilizer used in 15.612, and its cyclic module has length
N=(p²−1)/2 rather than the prime projective-line code length classified in
that paper.  Importing full PSL invariance would assume the missing step.

Several tempting finite shortcuts are already falsified.  Reflection of the
upper interval does not pair bad with good.  For example the Φ_3 orbit has
both members bad in reflected pairs at p=11,29,31,41,43,53,59.  The simple
XOR of all admissible contents also fails some orbit-products.  The first
two normalized endpoints u=(p+1)/2,(p+3)/2 clear every orbit through p=103,
but fail together at p=107: λ=1 gives t=2 and t=2/3=72, and both contents
are divisible by Φ_3 (and by no other factor of g).  Any proof must use the
whole residue/character structure rather than a fixed two-map disjunction.

`scripts/w2_phi3_endpoint_scan.py` removes factorization entirely and scans
only the universal Φ_3 gate.  Over all 166 primes 5≤p≤997, with the first
eight endpoints available, it finds no eight-endpoint failure but does kill
the next finite repairs:

* p=263 has bad-bad-bad-good at offsets 0,1,2,3;
* p=499,599,811 have bad-bad-bad-bad-good at offsets 0,1,2,3,4.

Thus the first three endpoints are not a p-law either, and the growing bad
prefix makes a small fixed disjunction a poor proof target.  The scan's
maximum first-good offset is 4; this is census evidence, not a uniform bound.

### Collective boundary norms: a stronger factor-free reduction

The unit-content boundary scan leaves useful information on the floor.  For
each valid translated boundary pair let

    N_a = c_a(X)c_a(X^-1) mod g.

If the polynomial gcd of `g` and all the `N_a` is 1, then for every
irreducible `f|g` at least one genuine U-difference has `f` not dividing its
content.  In particular no Aut orbit-product can divide every content, so
W2 follows for that prime.  This is stronger than W2, but weaker than asking
one `c_a` to be a unit: different rows may remove different factors.

The running-gcd backend now implements this test directly, with no
factorization.  On all 302 primes `5<=p<=2003`, the boundary norms generate
the unit ideal except at the already isolated `p=5`; the direct normalized
pole witness handles that prime.  Thus W2 is computationally certified
through `p=2003`.  The most delayed clearance is `p=1721`, where Φ_3 divides
the first fourteen norms and `a=15` clears it.  In 51 primes the common gcd
becomes 1 before any scanned row is itself a unit, concretely demonstrating
that the collective certificate is the right object.

There is substantial but nonuniform cyclotomic structure.  Of the 78 primes
whose common gcd after `a=1,2` remains nontrivial, 74 residuals are repetition
polynomials `R_d=1+X+...+X^(d-1)`.  The four exceptions are not noise: at
`p=953,1613,1709,1721` they factor into complete order layers / reciprocal
pairs.  Every one of the 78 residuals was checked to be invariant under both
`X->X^p` and `X->X^-1`, exactly the Aut factor-orbit structure relevant to
15.612.  The earlier tempting statement “two rows leave one R_d” is therefore
false; the corrected target is a union of complete Aut factor-orbits.

A dedicated GF(4) boundary scanner evaluates Φ_3 without FFT or
factorization.  Across all 429 primes `5<=p<=3001`, its only translated
boundary failure is `p=5`, and its maximum first-good row is `a=15` at
`p=1721`.  This kills the apparent `a<=10` bound and gives no reason to expect
any small fixed prefix to be a theorem.

Compact results and hashes are in
`evidence/w2_translated_common_norm_5_2003.json`; full traces are archived at
`/mnt/storage/e1work/maxplus_p13/w2_translated_attack_2026-08-23`.

### Consecutive translations reduce Phi_3 to four affine lines

The boundary family has a much sharper finite-difference description.  Put
`m=(p-1)/2` and let

    B_a(r) = 1_{(r+a mod p)>m}.

Over F2 its consecutive difference is supported on exactly two levels:

    B_(a+1)(r) + B_a(r)
      = 1_{r=m-a} + 1_{r=-1-a}.

The switched-pole signs do not depend on `a`, while the two copies of the
untransformed base vector cancel in the antipodal pair.  Consequently
`w_(a+1)+w_a` is the xor of the pullbacks of these two affine lines under
the two antipodal Mobius maps.  If `F_t(r)` denotes one such line residue,
the substitution `z=t*y` in `x=y/(t*y-1)` also gives the exact antipodal
relation `F_(-t)(r)=F_t(-r)`: multiplication by `-1` has trivial Phi_3
phase on both point-orbits.  Thus only one pole's two component tables need
to be stored.  The support identity was checked bit-for-bit against directly
reconstructed rows at representative small, medium, and delayed primes,
including `p=29,59,1721`.

At a root `zeta` of `Phi_3`, the named cyclic generator has exactly one
nonzero square/nonsquare multiplicative Fourier component.  This follows
from its exact two-component norm `1`: every nonzero GF(4) element has norm
`1`, so precisely one of the two component norms contributes.  Hence one
direct value of `w_0`, together with the two component tables for one pole,
recovers the complete sequence

    c_0(zeta), c_1(zeta), ..., c_(p-1)(zeta)

in `O(p^2)` field operations.  The previous method rebuilt a vector of
length `p^2` separately for every translation.  The implementation is
`scripts/w2_translated_phi3_sequence.py`; it asserts closure of the
recurrence, vanishing of the inactive component, and the exact reflection
`c_(m-a)=c_a`.  Its output agrees with the direct boundary scanner through
`p=101`, including every first-good representative.

Parameterizing one affine line turns each table entry into a cubic
multiplicative-character sum of the shape

    sum_{x in F_p} eta((x+A)/(x+B)),

with the pole and line level absorbed into `A,B`.  Equivalently, this is a
cubic cyclotomic sum on a Baer subline of `P^1(F_(p^2))`.  This is the first
reduction in the attack that replaces the translated `p^2`-point Boolean
vectors by a one-dimensional finite-field object.  The literature and
repository searches found general work on Baer sublines and character
sums, but no theorem that gives the required parity/nonvanishing on this
particular admissible interval.

The complete `5<=p<=5003` sequence census covers 668 primes.  Its only
all-valid-bad case is the already exceptional `p=5`; every `p>5` has a
nonzero valid Phi_3 residue.  For every tested `p == 5 (mod 12)`, the
sequence is binary and obeys the striking exact count

    #{a in F_p : c_a(zeta)=0} = (2p-7)/3,
    #{a in F_p : c_a(zeta)!=0} = (p+7)/3.

All 170 tested primes in this congruence class obey the formula.  This is an
observed law, not yet a theorem.  More importantly, even a proof
of the global count would not finish the interval problem.  At `p=101`, 65
of all 101 translations are zero and 39 of the 49 valid translations are
zero.  Thus the false inequality "total bad < valid interval length" cannot
be repaired into the desired conclusion.  The reflected sequence also has
linear complexity `p-1` in the tested `p == 5 (mod 12)` cases, so the exact
count is not coming from a bounded-order recurrence.  Apart from
`a -> m-a`, searches found no affine symmetry that forces a valid bad point
to pair with a good point.

The durable proof target is therefore an *incomplete* cubic/Baer-subline
sum: show that the zero set of this character-parity sequence cannot cover
the valid interval `1<=a<m/2`.  A square-root discrepancy estimate of the
right strength would suffice in the `p == 5 (mod 12)` class, and the census
suggests such cancellation, but no proved estimate has yet been matched to
the characteristic-two parity being measured here.

### Revised gap

For `Phi_3`, prove that the incomplete cubic/Baer-subline parity sequence
just derived is nonzero somewhere on the valid half-edge.  For a general
Aut factor-orbit `O`, derive the analogous translated boundary
`c_a mod f_O` as a finite-field multiplicative Fourier/Jacobi sum and show
it is not identically zero there.  Since `O` is handled separately, no
factor-count or union-bound loss occurs.  The formal `a=0` antipodal pair is
a unit, but the growing Φ_3 prefixes and maximal linear complexity rule out
propagation by a bounded-row recurrence.  A unit-content witness remains a
convenient computational certificate, but is no longer the proof target.

W2 stays OPEN.  No proposition number or closure flag is claimed.
