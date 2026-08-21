# Signed ν on the |κ|=3 locus — exact dataset at p=5,7,11; reduction validated end-to-end

> **SCOPE CORRECTION (same day):** where this note says signed ν values are
> *required* and magnitude-only bounds insufficient, that is proved only for
> the **affine-parametrization** hypothesis set.  The w-line convolution in
> `NOTE_2026-08-21_nu_convolution_reduction.md` is a magnitude-only estimate
> that closes for p≥11.  Read that note for the current reduction.

Date: 2026-08-21.  Follows `NOTE_2026-08-21_leftover3_contraction_closure.md`
(prop 15.590).  Data: `evidence/nu_of_w_dataset.json` (45 exact values),
`evidence/p11_orbit_mu_nu_values.json`.  Scripts in
`/mnt/storage/e1work/leftover3_mu/` (`nu_p11.py`, `nu_of_w.py`,
`affine_p11_validate.py`, `proj_rank_p7.py`).

## Frame normalization

By 3-transitivity every four-set is equivalent to S = {∞,0,1,w}.  Then

    κ(S) = 1 + χ(w) + χ(1−w),      φ(S) = Σ_r χ(r(r−1)(r−w)) = −a_q(E_w)

with E_w the Legendre curve y²=x(x−1)(x−w) over F_q.  So the |κ|=3 locus
is exactly {w : χ(w)=χ(1−w)=1} (both-squares), and ν is a function
ν(w) = ½(m₄⁺−m₄⁻)({∞,0,1,w}) on it.  Locus sizes: 5, 11, 29 at p=5,7,11.

## The dataset (canonical values, no orbit-sign conventions)

ν(w) computed for EVERY w on the locus: p=5,7 from enumerated Max±
(V±-completion), p=11 from the stored pair Gram
(`G_pairmoment_p11.npy`, Sp) with Max− via the 15.254 swap transport (Sm).
Validations, all exact: Gram↔C consistency 0/500; ν(w)=0 on all nine
|κ|=1 orbits at p=11; max|μ| = 17827/1560713 = census.

Fiber structure at p=11 (5 fibers = the 5 signed orbits): generic
quadratic (12), quadratic (6), F_p-rational generic {3,4,5,7,8,9} (6),
harmonic {2,6,10} (3), equianharmonic w²−w+1=0 (2).  φ=−a_q confirms:
rational fiber has a_p²=16; harmonic and equianharmonic are both
supersingular (a_p=0, p≡3 mod 4) with IDENTICAL φ=22 — yet ν differs
(−500128 vs −417120 over N): **ν sees strictly more than the elliptic
trace**; the distinguishing invariant is the live question.

## Empirical laws (3 primes, 45 points — observations, not theorems)

1. ν(w) < 0 on the entire locus at p=5,7,11 (canonical, frame-fixed).
2. μ(w) > 0 on the entire locus (all κ=+3; no κ=−3 orbit exists at
   p=5,7,11 in this frame family).
3. max|ν|·p⁴ ≈ 96.2, 109.0, 98.6 at p=5,7,11 — |ν| ~ c/p⁴ scaling,
   two orders below the trivial bound.  NOT fitted to a formula; recorded
   only as scaling.

## Reduction validated end-to-end at p=11

The data-free degree-4 affine map (kgen convention, rank 15/19, 4 free
δ-columns) + the 4 free signed ν values reproduces ALL 14 orbit μ values
at p=11 with Fraction equality (`affine_p11_validate.py`: 14/14 OK).
Same at p=5, 7 earlier.  So the chain

    {signed ν on |κ|=3 fibers}  →  affine map  →  μ on |κ|=1  →  leftover 3

is exact at every prime with data.  Remaining: closed form (or sign-exact
estimate) for ν(w) at general p.

## Degree-escalation verdict, sharpened

The joint degree-4+6 kernel at p=7 (dim 4) projects onto the (μ₄,δ₄)
block with rank 2 — equal to the degree-4 kernel dim (`proj_rank_p7.py`).
Degree 6 pinned 1-of-1 δ-combination at p=5 but 0-of-2 at p=7: its
marginal value on the 4-point sector is already zero at p=7.

## Next attacks on ν(w)

- The harmonic/equianharmonic split at equal φ says the needed invariant
  is finer than a_q: candidates: the class of E_w's twist pair under the
  15.267 V₄/pairing-determinant structure, or Jacobi-sum arguments
  (15.496 territory), or the CM discriminant.
- 15.268's pairing argument gives ν=0 on |κ|=1; its |κ|=3 analogue (all
  three pairing determinants +1, all pole-squares) may yield exact
  ν-relations instead of vanishing — the first thing to derive.
- Any conjectured ν(w) form must reproduce all 45 values and then be
  tested out-of-sample at p=13 via `p13p17enum.md`'s Route A (orbit-rep
  moments), NOT fitted at 3 primes and shipped (METHOD note).
