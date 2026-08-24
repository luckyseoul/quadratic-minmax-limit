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

### Frobenius and projective-line parity solve the Phi_3 gate at p=5 mod 12

There is one clean congruence-class phenomenon hidden by the mixed-prime
endpoint census.  Suppose `p == 5 (mod 12)` and choose the first normalized
endpoint

    u=(p+1)/2,        t=lambda/u=2 lambda in F_p.

The named halfspace is invariant under `p`-Frobenius.  In the repository's
quadratic basis `omega^2=d` with `d` a nonsquare, writing `x=a+b omega`
shows that the level of `omega^-1 x` is `a/d`, unchanged when
`x^p=a-b omega`.  Since `t` lies in `F_p`, both the pole involution and its
quadratic switch commute with Frobenius.  The endpoint difference `w` is
therefore Frobenius invariant.

Write the nonsquare multiplicative orbit as `omega_0 g^j`, where `g` is the
square generator, and let `C_j` be the number of support points with exponent
class `j mod 3`.  Frobenius acts on the exponent by

    j -> (p-1)/2 + p j = 2-j mod 3.

Consequently `C_0=C_2`.  Class 1 is preserved, but the nonsquare orbit has no
`F_p`-points, so `C_1` is even.  At a root `Z` of `Phi_3`, the endpoint
residue is thus

    C_0 + C_1 Z + C_2 Z^2 = (C_0 mod 2) Z.

The named cyclic generator has the exact nonsquare class counts

    ((p+1)/6, (7p-11)/6, (p+1)/6),

and hence residue `Z`, since `(p+1)/6` is odd.  The normalized endpoint
content is therefore exactly `1` if and only if `C_0` is odd.  This reduces
the complete `Phi_3` endpoint question in this congruence class to one
integer parity, with the Frobenius equalities proved before computation.

Projective-line parity closes the remaining lemma.  Since `6|(p+1)`, scalar
multiplication by `F_p^*` preserves every sextic class in `F_(p^2)^*`.
Class 0 therefore consists of `(p+1)/6` projective `F_p`-directions, an odd
number.  Consider any finite nonsquare direction

    x=c(1+k omega),   c in F_p^*,   n=1-d k^2 nonsquare.

The base halfspace level is `c/d`, so its line parity is `(p-1)/2=0 mod 2`.
For the endpoint pole, after clearing the square denominator `d^2`, the
quadratic switch is controlled by

    Delta(c)=4n c^2-4d c+d^2.

Its discriminant `16d^3 k^2` is nonsquare.  The standard quadratic-character
sum is therefore `sum_c chi(Delta(c))=-chi(4n)=1`, so exactly `(p-1)/2`
switches occur, again even.  The pole-image halfspace level is

    R(c)=(2n c^2-dc)/Delta(c).

For a target level `r`, its fiber is

    2n(2r-1)c^2+d(1-4r)c+r d^2=0.

At the upper-half boundary `u=1/2=(p+1)/2`, this becomes a nonconstant
linear equation and has exactly one solution.  Away from `u`, a quadratic
fiber has odd cardinality only when its discriminant vanishes.  After the
irrelevant square factor `d^2`, that discriminant is

    Q(r)=1+8d k^2 r(2r-1).

It obeys `Q(u-r)=Q(r)`.  Its roots have no fixed point under this involution,
because `Q(u/2)=n!=0`; neither `0` nor `u` is a root.  The involution pairs
roots within the lower half or within the upper half.  Hence the upper half
contains an even number of branch values, while the unique linear fiber at
`u` remains.  Every finite nonsquare projective direction therefore has odd
endpoint-support parity.

The omitted direction `F_p^* omega` is a nonsquare cube: `omega^(p-1)=-1`
and `(p+1)/3` is even.  It lies in exponent class 1, not class 0.  Thus all
`(p+1)/6` class-0 directions are finite and contribute odd parity.  Their
number is odd, proving `C_0` odd.  Consequently the first normalized
endpoint has `Phi_3` content exactly `1` for every prime
`p == 5 (mod 12)`.

The same line parity clears a complete factor layer, not only `Phi_3`.
Put `M=(p+1)/2` and write the endpoint's nonsquare Fourier sequence as

    W_j=w(omega_0 g^j),       0<=j<N=(p-1)M.

The subgroup `F_p^*` is generated by `g^M`.  Folding modulo `X^M+1` therefore
xors along projective directions:

    Wbar_r = sum_(h=0)^(p-2) W_(r+hM).

Every finite nonsquare direction has parity one, while the sole
infinite-slope direction `F_p^* omega` has `(p-1)/2=0 mod 2` support points.
If its quotient index is `r_0`, then exactly

    Wbar(X) = R_M(X)+X^(r_0),
    R_M(X)=1+X+...+X^(M-1).

For every nontrivial `M`-th root `alpha`, this gives
`Wbar(alpha)=alpha^(r_0)!=0`.  The multiplicative Fourier identity
`Wbar(alpha)=c(alpha) Gamma(alpha)` then forces `c(alpha)!=0`; no separate
nonvanishing theorem for the generator component is needed.  Finally `M` is
odd and divides the odd part of `N`, so `R_M|g`.  Thus the same first
normalized endpoint is coprime to the entire projective-torus factor layer

    R_((p+1)/2) = (X^((p+1)/2)+1)/(X+1)

for every prime `p == 5 (mod 12)`.  The earlier `Phi_3` theorem is its
degree-two special case.

This also identifies the exact unresolved layer.  If
`s=oddpart(p-1)`, then the odd part of `N` is `Ms` and

    g=R_(Ms),             g/R_M=R_s(X^M).

Thus only the `p-1` layer `R_s(X^M)` remains for this endpoint.  In
particular, whenever `p-1` is a power of two, `s=1` and `g=R_M`; the first
normalized endpoint is then a unit-content W2 witness, not merely a witness
against one factor layer.  Exact Krylov reconstruction at `p=5,17,29,41`
independently gives `gcd(c,R_M)=1`, checking the Fourier convention against
the actual content polynomial.

`scripts/w2_phi3_endpoint_frobenius.py` checks the coordinate formula,
Frobenius action, projective direction counts, line parities, quotient fold,
and `R_M|g` directly.
Through `p=5003`, all 170 primes in this congruence class have zero failures.
The reproducible in-repository result is
`evidence/w2_phi3_endpoint_frobenius_5_2003.json`; the 5003 mesh extension is
archived beside the translated sequences under
`/mnt/storage/e1work/maxplus_p13/w2_translated_attack_2026-08-23` as
`w2_projective_torus_endpoint_5_5003.json`.

### Reciprocal-orbit reduction and the two-endpoint repair

The residual `p-1` layer does not admit the strongest possible continuation:
the first normalized endpoint is not always a W2 witness.  The obstruction can
nevertheless be isolated without factoring the ambient repetition polynomial.

Keep `p == 5 (mod 12)`, `M=(p+1)/2`, `s=oddpart(p-1)`, and `H=Ms`.  Frobenius
acts on the square and nonsquare exponent sequences by

    j -> p j,                 j -> (p-1)/2 + p j,

respectively.  Since the endpoint difference is Frobenius invariant, vanishing
of either two-component multiplicative Fourier transform is invariant under
`alpha -> alpha^p`.  Therefore the bad irreducible factors of the content are
unions of `p`-orbits.  The full Aut orbit adds inversion, so endpoint W2 can
fail only when a bad `p`-orbit and its reciprocal `p`-orbit are both bad.

The exact generator norm from the Bose relative-difference-set identity makes
this reduction computationally sharp.  If `c` is the endpoint content, the
xor of the square and nonsquare folded autocorrelations is exactly

    c(X)c(X^-1) mod g.

Consequently

    d = gcd(g, c(X)c(X^-1))

contains exactly the reciprocal closure of the bad factors.  Factors outside
`d` cannot participate in an endpoint W2 failure.  Thus one first computes
`d` without factorization, factors only `d`, evaluates the endpoint Fourier
residue on those factors, and forms their orbits under `X->X^p` and
`X->X^-1`.  The implementation is
`scripts/w2_endpoint_reciprocal_orbits.py`.

Across all 77 primes `p == 5 (mod 12)`, `5<=p<=2003`:

* the first endpoint has unit content at 66 primes;
* it is already a W2 witness at 72 primes;
* its five genuine W2 failures are `p=401,461,1061,1289,1877`;
* the second endpoint `u=(p+3)/2` clears every factor in the exceptional
  reciprocal closure at all five failures.

The first two failures are the self-reciprocal `Phi_5` factor.  At `p=1061`
the obstruction is the reciprocal degree-four pair in `Phi_15`; at
`p=1289,1877` it is the reciprocal cubic pair `Phi_7`.  Their norm-gcd degrees
are respectively `4,4,8,6,6`.  These examples disprove the universal
one-endpoint conjecture, but replace the huge residual layer by three explicit
low-order phenomena in the tested range.  The two-endpoint statement is still
census evidence, not a theorem.

The literature check did not duplicate this reduction.  Wu--Yue--Fan
[classify self-reciprocal irreducible factors of `X^n-lambda`](https://arxiv.org/abs/2001.04766),
which names when singleton reciprocal obstructions can occur but gives no
nonvanishing result for these endpoint residues.  Chai--Li's
[Soto--Andrade character-sum framework](https://www2.math.upenn.edu/~chai/papers_pdf/chaili_v34_part2.pdf)
studies norm-one-torus character-sum families and their monodromy; it likewise
does not decide the binary vanishing of this halfspace/pole sum.  The remaining
proof target is now precise: show that the two consecutive normalized endpoint
residues cannot vanish on a complete orbit generated by `p` and inversion in
the residual `R_s(X^M)` layer.

### Oriented content reconstruction removes factorization entirely

The reciprocal norm is not the strongest consequence of the Bose generator
identity.  Let `W_sq,W_ns` and `Gamma_sq,Gamma_ns` be the two multiplicative
Fourier components of `w=c(D)gamma`.  At every root of `g`,

    W_sq(alpha) = c(alpha) Gamma_sq(alpha),
    W_ns(alpha) = c(alpha) Gamma_ns(alpha),

while the exact generator norm says

    Gamma_sq(alpha) Gamma_sq(alpha^-1)
      + Gamma_ns(alpha) Gamma_ns(alpha^-1) = 1.

Multiplying and adding therefore reconstructs the oriented content itself:

    c(X) = W_sq(X) Gamma_sq(X^-1)
         + W_ns(X) Gamma_ns(X^-1)                 mod g.

The right side is two folded cross-correlations.  It requires neither a
Krylov coordinate solve nor factorization.  Direct Krylov reconstruction at
`p=5,17,29,41` agrees coefficient-for-coefficient, fixing the correlation
orientation as well as the gcd result.

For a Frobenius-invariant endpoint put

    A(c) = gcd(g, c(X), c(X^-1)).

Bad factors are already unions of `X->X^p` orbits, so `A(c)` is exactly the
product of complete Aut orbits missed by that endpoint.  Hence the two
consecutive endpoints certify W2 exactly when

    gcd(A(c_0), A(c_1)) = 1.

This is implemented by `scripts/w2_two_endpoint_oriented_content.py`.  On the
same 77 primes through `p=2003`, endpoint zero individually clears 72 primes,
endpoint one clears 46, and the oriented gcd certifies all 77 pairs.  The
first endpoint failures remain exactly `401,461,1061,1289,1877`; there are no
two-endpoint failures.  Unlike the earlier norm/factor scan, this conclusion
is a direct polynomial certificate and does not inspect any irreducible
factor.  The split exact outputs are
`evidence/w2_two_endpoint_oriented_content_5_1000.json` and
`evidence/w2_two_endpoint_oriented_content_1001_2003.json`.

The residual Aut-bad layers in the `p<=2003` census are very small.  Endpoint
zero has only orders `5,7,15` in its five failures there.  After removing the
projective layer, endpoint one has residual misses only at `p=401` (the
reciprocal `Phi_15` quartics), `p=953` (the reciprocal `Phi_7` cubics), and
`p=1613` (the reciprocal order-31 quintics; its additional `Phi_3` factor is
projective).  No residual layer is shared by the two endpoints at the same
prime.  This is census structure, not yet a proof that higher orders cannot
appear.

### The second endpoint's complete projective quotient

The projective-line parity argument extends exactly to the second endpoint
`u=(p+3)/2=3/2`.  Write `tau=1/u=2/3`.  On a finite nonsquare direction
`x=c(1+k omega)`, `n=1-dk^2`, the image level is

    R(c) = (tau n c^2-dc)
           /(tau^2 n c^2-2 tau d c+d^2).

The fiber above level `r` is

    n tau(r tau-1)c^2 + d(1-2r tau)c + r d^2 = 0,

and its discriminant, apart from the square `d^2`, is

    Q(r)=1+4 tau(n-1)r(1-tau r),       Q(u-r)=Q(r).

The involution `r->u-r` preserves the lower and upper half-intervals except
for the single crossing pair `{1,u-1}`.  The linear fiber at `r=u` contributes
one.  Consequently the line support is even precisely when

    Q(1)=0  iff  d k^2=9/8.

These two exceptional slopes exist exactly for `p=5 (mod 24)`; they are
themselves nonsquare directions.  For `p=17 (mod 24)` there are none.  The
base and switch terms have even parity, and the pure `F_p^* omega` direction
is also even: its switch count is even, while its image level depends on
`c^2` and therefore has even fibers.

Let `r_inf` be the pure-omega quotient index.  Frobenius sends a quotient
index to `-r-1`, fixes `r_inf`, and exchanges the two exceptional slopes.
Thus their indices are `r_inf+a,r_inf-a`.  The second endpoint quotient is

    R_M + X^r_inf                                      (p=17 mod 24),
    R_M + X^r_inf + X^(r_inf+a) + X^(r_inf-a)          (p=5 mod 24).

At a nontrivial `M`-th root the second formula is a monomial times
`1+X^a+X^(2a)=Phi_3(X^a)`.  Therefore every possible second-endpoint zero in
the projective layer is explicit, and endpoint zero's already-proved
nonvanishing on all of `R_M` clears it automatically.  The executable check
is `scripts/w2_second_endpoint_projective.py`.  Its independent scan of all
77 eligible primes through `p=2003` has no pure-direction, finite-slope,
symmetry, quotient, or bad-polynomial mismatch; the output is
`evidence/w2_second_endpoint_projective_5_2003.json`.  The genuinely
unresolved two-endpoint theorem is now confined to the scalar `p-1` layer.

### GPU spectral census: order 17 appears

The factor-free calculation admits a useful GPU form.  Fold each square and
nonsquare orbit sequence modulo the odd order `H` *before* taking a Fourier
transform.  This is the same calculation in `F_2[X]/(X^H+1)`, but removes the
full 2-part of `(p^2-1)/2` from the transform length.  All requested endpoints,
both point-orbits, and the two Bose-generator components can then be processed
as one spectral batch.  The oriented products recover `c(X)`; the endpoint
autocorrelations independently recover `c(X)c(X^-1)`.  Thus one run tests four
different statements:

1. exact collective W2 via the gcds `A(c)`;
2. the stronger collective norm/unit-ideal condition;
3. the projective factor layer;
4. the residual scalar factor layer.

This is implemented in `scripts/w2_endpoint_norm_gpu.py`.  It reuses the
Wieferich engine's persistent-allocation principle, transfers endpoint batches
rather than individual transforms, writes evidence atomically, and uses a tiny
NTL bridge for compiled `GF(2)` polynomial gcds.  At `p=1997` on the RX 9070 XT,
these changes reduced the four-endpoint job from `53.8s` to `11.7s`.  V100
double precision and RX 9070 XT single/double precision give identical exact
certificates on cross-checked cases; every run also asserts the exact Bose
generator norm, so unsafe FFT rounding fails loudly.

The independent GPU rerun exactly matches every endpoint degree and pair
witness in the earlier 77-prime census through `p=2003`.  The new census adds
all 35 eligible primes `2004<=p<=3000`.  The first two endpoints again certify
W2 at every prime, so the exact pair now has no failures on all 112 eligible
primes through `p=3000`.  Four endpoint norms generate the unit ideal throughout
the new range as well.

More importantly, the wider scan falsifies the tempting bounded-order pattern.
Endpoint zero has five new failures:

    p=2081: reciprocal order-15 quartics, degree 8;
    p=2381: Phi_5 times a self-reciprocal order-17 factor, degree 12;
    p=2441: Phi_5, degree 4;
    p=2549: reciprocal order-7 cubics, degree 6;
    p=2621: Phi_5, degree 4.

At `p=2381`, the degree-eight factor is

    X^8+X^7+X^6+X^4+X^2+X+1,

whose roots have order 17.  This is the first endpoint-zero order outside
`{5,7,15}`.  At `p=2549,2621`, endpoint one misses only the projective `Phi_3`
while endpoint zero misses only the displayed scalar factor, making the layer
complementarity explicit.  Repository searches found no prior W2 record of
these primes.  Exact outputs are
`evidence/w2_four_angle_gpu_2004_2500.json` and
`evidence/w2_four_angle_gpu_2501_3000.json`.

### Sparse atomic low-order census: the endpoint pair first fails at p=5237

The dense spectral calculation spends most of its time constructing the full
polynomials even when the question concerns a small root order.  For an odd
order `d`, reduction modulo `X^d+1` needs only the parity of each endpoint
orbit in exponent classes modulo `d`.  The Bose-generator components admit
the same reduction.  Two cyclic cross-correlations of these tiny residue
vectors recover the oriented content modulo `X^d+1`, after which

    gcd(Phi_d, c(X), c(X^-1))

is the exact primitive-order-`d` Aut obstruction.  No FFT, factorization, or
probabilistic test is involved.

`scripts/w2_low_order_atomic_gpu.py` implements this reduction with a
CUDA/HIP grid-stride kernel.  Each block accumulates endpoint residue parities
with shared-memory `atomicXor`, then atomically merges only the small residue
table.  The multiplicative square/nonsquare orbits and inverse logarithm are
generated in the compiled `src/gf2x_ntl.cpp` bridge, removing the former
Python orbit loop.  The same source runs unchanged on the V100 and RX 9070 XT.

`scripts/w2_low_order_atomic_validate.py` compares sparse records directly
with dense exact endpoint gcds.  Through `p=3000`, every primitive odd order
`3<=d<=63` matches coefficient-for-coefficient: 112 primes, 1,246
endpoint/order comparisons, and all 60 exceptional records.  At `p=2381`,
the sparse RX 9070 XT run takes 0.182 seconds end-to-end (0.011 seconds in the
atomic kernel), versus 17.013 seconds for the earlier dense exact run.

The extended scan tests every odd `d<=255` dividing the ambient odd order at
all 197 eligible primes `3001<=p<=10000`.  It records 114 endpoint/order
exceptions, 29 in the unresolved scalar layer.  New higher-order scalar
examples include order 85 at `p=5441`, order 255 at `p=6869`, order 63 at
`p=8513`, and orders 51 and 255 at `p=9929`.  Repository and history searches
found no earlier W2 record of these cases.

Most importantly, the first two normalized endpoints share an obstruction:

    p=5237, endpoint 0: Phi_7 Phi_11,   aut_bad = 0x157d5, degree 16;
    p=5237, endpoint 1: Phi_3 Phi_7,    aut_bad = 0x17d,   degree 8;
    common pair gcd:    Phi_7,          aut_bad = 0x7f,    degree 6.

Thus the exact two-endpoint statement that held through `p=3000` is false.
An independent full-polynomial V100 run confirms both endpoint gcds and the
common `Phi_7`; this is not an inference from the bounded-order scan.  The
next endpoint repairs the failure completely: endpoint offset 2 has
`aut_bad=1` and is itself a unit-content witness.  The exact prefix bad degrees
are therefore `16,6,0` for one, two, and three endpoints.

The order-63 hit at `p=8513` was also independently checked with the dense
RX 9070 XT calculation.  Endpoint zero has exact scalar Aut-bad polynomial
`0x13f9` of degree 12, while endpoint one is a unit-content witness and the
pair succeeds.  The broad sparse census proves only the reported primitive
orders through 255; it does not exclude additional higher-order factors or
additional pair failures in `3001<=p<=10000`.  Accordingly these findings do
not change W2's open status.

Raw sparse outputs are
`evidence/w2_low_order_atomic_gpu_5_1500_o63.json`,
`evidence/w2_low_order_atomic_gpu_1501_3000_o63.json`,
`evidence/w2_low_order_atomic_gpu_3001_8000_o255.json`, and
`evidence/w2_low_order_atomic_gpu_8001_10000_o255.json`.  The two independent
dense confirmations are `evidence/w2_endpoint_norm_gpu_p5237_triple.json` and
`evidence/w2_endpoint_norm_gpu_p8513.json`.

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

### General four-line differences isolate a three-step scalar target

The four-line reduction is not special to `Phi_3`.  For every odd order `d`,
fold each pullback line by exponent class modulo `d`, reconstruct the oriented
content with the two Bose-generator correlations, and form the running gcd

    gcd(Phi_d, Delta c_1, Delta c_1^*,
                 Delta c_2, Delta c_2^*, ...),

where `Delta c_a=c_(a+1)+c_a`.  If this gcd is one, no primitive-order-`d`
Aut orbit can kill every genuine boundary row.  Unlike the original contents,
the differences contain no interval halfspace at all: each is exactly four
Baer-subline incidence vectors.

`scripts/w2_boundary_line_difference_gpu.py` computes all line-level residue
tables in one `O(p^2)` CUDA/HIP atomic pass and then runs the tiny exact gcd
chain.  Its `Phi_3` result agrees with the independent full translation
sequence at `p=17`, including the first nonzero difference `c_3+c_2`.

Across all 309 primes `p == 5 (mod 12)`, `5<=p<=10000`, every odd order
`d<=255` dividing the ambient odd order is cleared by valid four-line
differences, except the degenerate `p=5` boundary with no consecutive valid
rows (the known direct normalized witness handles it).  This comprises 2,663
primitive orders.  The existing projective-line theorem already handles the
1,393 projective orders.  Of the 1,270 genuinely scalar orders, 1,249 clear
after `Delta c_1`, 20 after `Delta c_2`, and exactly one after `Delta c_3`:
`p=2141,d=15`, whose residual chain is

    Phi_15, Phi_15, 1.

Thus the new finite theorem candidate in this congruence class is precise:
prove that the first three explicit four-line differences have no common
factor in the scalar layer `R_s(X^M)`.  This is materially different from the
false fixed-endpoint conjectures: it is an algebraic incidence statement with
the interval term cancelled before the gcd.  The census still tests only
orders through 255 and is not a proof; higher scalar orders could falsify the
three-step statement.

Evidence is
`evidence/w2_boundary_line_difference_gpu_5_5000_o255.json`,
`evidence/w2_boundary_line_difference_gpu_5001_10000_o255.json`, and the
exceptional trace
`evidence/w2_boundary_line_difference_gpu_p2141_o15_trace.json`.

The earlier durable proof target was an *incomplete* cubic/Baer-subline
sum: show that the zero set of this character-parity sequence cannot cover
the valid interval `1<=a<m/2`.  A square-root discrepancy estimate of the
right strength would suffice in the `p == 5 (mod 12)` class, and the census
suggests such cancellation, but no proved estimate has yet been matched to
the characteristic-two parity being measured here.

### Revised gap

For `p == 5 (mod 12)`, prove or falsify beyond order 255 that the first three
four-line differences have scalar common Aut-gcd one.  The projective layer
is already proved, so this would finish this congruence class.  If the
three-step statement fails at higher order, the fallback remains the full
incomplete cubic/Baer-subline problem: for each Aut factor-orbit `O`, show
that the translated boundary sequence is not identically zero on the valid
half-edge.  Since `O` is handled separately, no factor-count or union-bound
loss occurs.  Other prime congruence classes still require their corresponding
projective/scalar analysis.  A unit-content witness remains a convenient
computational certificate, but is no longer the proof target.

W2 stays OPEN.  No proposition number or closure flag is claimed.
