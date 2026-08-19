# Strategy for the remaining E(1) leftovers (post-15.588)

State: 15.588 shipped (28 tests green, no flag flipped). Two structural
results now stand under everything: (I) Max+ = polynomial line-sum profiles,
stratified by the number k of active square directions, with the repo's
"unclassified full family" identified as the k >= 4 strata; (II) the floor
lambda_min(Phi) >= 6 is `<T, chi>/<chi, chi> >= 6` per irreducible constituent
of Z, for ONE class function T(g) = E_y <y, gy>^2 - 2n on Aut(C), whose
constituent decomposition is already proved for general p (15.278 F).

The wall ("Max+ moments at general p") is therefore two concrete objects now:
per-stratum counts/moments (Part I) and the class function T (Part II). They
meet: T is a Max+ average, so it decomposes over the strata,
  T = sum_k (N_k / N) T_k,
and each T_k is an average of <y, gy>^2 over an EXPLICIT polynomial family.

## Step 1 — finish p = 11 (compute, hours)

k=6 via dilation-gauged outers x GPU inner (see scripts/maxplus_profile_enum/
README): ~111 + ~1342 orbit reps at seconds each. Then:
- N(11), D(11) = N/44 for the first time past p = 7.
- Phi spectrum at p = 11 (spectrum.py; check D*lambda integrality, the
  floor margin 8 - lambda_min, bottom multiplicity n?, top n/2?).
- |kappa|=1 four-point mu at p = 11 vs L = 9/242 (moments.py).
This does not flip a flag; it pins the p-laws to name.

## Step 2 — close the k=4 stratum as a theorem (math, the model case)

The k=4 family mod translations is pure parabolas rho_j = lam w0_j s^2 + v_j
with a two-adjacent-carry condition — a quadratic-character counting problem.
Data so far: per (subset, lam) the v-count is 15 (p=7, harmonic), 4/2
(p=11, generic/harmonic), 2/0 (p=13, only 7 of 35 subsets nonempty).
The counts are collapsing with p and are cross-ratio-functions, not
p-polynomials; get p = 17, 19, 23 (k4_scan.py, cheap), classify the passing
(lam, v) by quadratic invariants, then prove the selection rule by Gauss-sum
manipulation. This is the template: the SAME carry analysis at general k is
what a full classification needs.

## Step 3 — T stratum by stratum (the join, and the real prize)

T_1 (the 1d family) is a closed-form object: y = v(t(x)) with v a +-1
function on F_p, and <y, gy> is a correlation of two 1d profiles along the
g-twisted pairing — expressible by Jacobi/Kloosterman-type sums. T_3 is the
affine family — same machinery one level up. If the k >= 4 strata are
o(N) as p grows (the k=4 collapse at p=13 and the 1d family's 2^p growth
both point that way, but k=5 at p=11 is 1.3M vs k4's 58k, so this is NOT
settled — the p=11 k=6 count is the test), then
  T = T_{1d} + O(poly(p)/2^p corrections),
and the floor becomes: main term (explicit character sums) >= 6 + error
bound. The ~25% slack in both leftovers 1 and 3 for p >= 7 is exactly the
room such an asymptotic argument needs, with p <= P0 finite-checked by the
enumerator. That is the intended shape of the eventual flag flip:
- leftover 1: lambda_pi = <T, chi_pi>/<chi_pi, chi_pi> >= 6 for each of the
  finitely many constituent types (principal series / Weil / Steinberg),
  main term 8 with error O(p^{-1/2}) from Weil bounds on the stratum sums.
- leftover 3: mu on |kappa|=1 four-sets = same stratum decomposition of a
  4-point average; needs |mu| <= (p-2)/2p^2, margin ~25%.

## Step 4 — which strata survive as p grows (the decisive unknown)

Get exact stratum totals at p = 11 (k=6 pending) and p = 13 (k up to 7;
k=7 needs degree-5 profiles — the same toolkit, one more gauge level, and
C(7,k) subsets). Two scenarios:
- If max-k strata keep exploding (k=5 at p=11 suggests the k ~ m families
  dominate), the asymptotic route needs cancellation in T_k, not smallness
  of N_k — harder, but T_k is still an explicit polynomial-family average.
- If high-k strata collapse the way k=4 does by p = 13, the 1d+affine
  domination argument closes the wall outright.

## Discipline

No census as a p-law. Every count above is exact enumeration with Cy = py
verified at residual 0. A flag flips only when the general-p inequality is
proved; the enumerator's role is to name the law and to supply the finite
checks below the crossover prime.
