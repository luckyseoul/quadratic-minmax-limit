# Leftover 3: eigen-contraction closure of μ/δ — determination at p=5, level-4 kill

Date: 2026-08-21.  Companion code: `src/e1_gmin_m4_prop15590.py`,
tests `tests/test_prop15590.py`.  Lab artifacts (arrays, exploratory
scripts, PSD scans): `/mnt/storage/e1work/leftover3_mu/` on soulkiller.
**No leftover flag flipped.**  Current correction: the bounds `|μ|≤L` and
15.191 K's stronger `|μ|≤2/n` are sufficient only for the `|κ|=1`
half of the multi-level Type-I reduction.  The independent `|κ|=3` signed
`(μ,ν)` inequality and its nonparticular δ remainder also remain open.

## New Max-free identities (exact, zero violations at p=5 and p=7)

With μ=½(m₄⁺+m₄⁻), δ=½(m₄⁺−m₄⁻) and only Cy=±py, y_i²=1,
E±[y_iy_j]=±C_ij/p (15.189):

    (★μ)   Σ_{l∉{i,j,k}} C_kl μ({i,j,k,l}) = C_ij            (all triples)
    (★δ)   Σ_{l∉{i,j,k}} C_kl δ({i,j,k,l}) = −(2/p) C_ik C_jk
    (out)  Σ_{b∉{i,j,a},b≠l} C_lb μ({i,j,a,b}) = p δ({i,j,a,l})
           Σ_{b∉{i,j,a},b≠l} C_lb δ({i,j,a,b}) = p μ({i,j,a,l})
                       − (1/p)(C_li C_ja + C_lj C_ia + C_la C_ij)

plus the degree-6 analogues coupling (μ₆,δ₆) down to (μ,δ) (prop 15.590 C).
The l∈{i,j} boundary terms *cancel between the ensembles* in (★μ) — the
inside contraction closes on μ alone.  Also: 1ᵀy = (1±p)y_∞ pointwise
(so E±[s²], E±[e₄] are exact), and δ is **twisted-dead** on every |κ|=1
signed orbit at p=5,7 — m₄⁺=m₄⁻ on |κ|=1 by sign-consistency alone
(independent mechanism from 15.268's pairing argument).

## Determination theorem at p=5

Signed-orbit counts of four-sets under ⟨PΓL(2,q), anti⟩: p=5 has 4 orbits,
p=7 has 7.  The **complete** equivariant linear system (all relations above,
built from marked-set orbit representatives — provably complete by
equivariance) yields, by exact Fraction elimination:

| system | p=5 | p=7 | p=11 | p=13 |
|---|---|---|---|---|
| μ-orbits of four-sets | 4 | 7 | 14 | 19 |
| … of which \|κ\|=1 | 2 | 4 | 9 | 12 |
| live δ-orbits (= \|κ\|=3 orbits) | 2 | 3 | 5 | 7 |
| degree-4 kernel dim K₄ | 1 | 2 | **4** | **6** |
| joint degree-4+6 kernel dim | **0** | 4 | — | — |

The p=11 row is **data-free**: the degree-4 system is homogeneous in the
orbit variables, so K₄(11) needs only orbits + sign cocycles, never Max±.
It was computed independently of any enumeration (`kernel_p11.py`), and the
same script reproduces K₄(5)=1 and K₄(7)=2 — the values obtained by the
separate data-driven route — as a cross-check.  ("3-set orbits = 1" at every
prime is the correctness canary: PGL(2,q) is 3-transitive.)

Two exact structural facts, verified **set-wise** at p=5,7,11:

* **The kernel is carried entirely by δ**: every free column is a
  δ-coordinate, and K₄(p) = (#live δ-orbits) − 1.  So the degree-4
  contraction system already determines μ as an affine function of the
  δ-orbit values and pins exactly one linear combination of them.
  Leftover 3 at degree 4 is therefore *not* an unstructured search — it is
  the problem of bounding a K₄(p)-dimensional δ-vector (dim 4 at p=11).
* **{δ-dead orbits} = {|κ|=1 μ-orbits}** as sets (2/2, 4/4, 9/9).  This
  re-derives m₄⁺=m₄⁻ on |κ|=1 at p=11, where the census route is out of
  reach, from sign-consistency alone.

At p=5 the joint system has rank = #unknowns = 38 (4 μ + 2 δ + 13 μ₆ +
19 δ₆ live orbit values): **μ, δ, μ₆, δ₆ at p=5 are the unique equivariant
solution — no census needed.**  The solved values reproduce the enumerated
moments exactly (per-set array equality).

The degree-4 kernels touch the |κ|=1 μ-coordinates at both primes, and
implied-functional tests (e₄-sum, ⟨κ⟩, ⟨φ⟩, s²-pair sums) all vanish on
the kernel — the linear theory really is complete at each degree; the
kernel is not a missing-equation artifact.

## Route kill: level-4 moment/SoS cannot close leftover 3

The degree-4 relaxation {complete linear theory + signed equivariance +
M± = E±[vec(yyᵀ)vec(yyᵀ)ᵀ] ⪰ 0, compressed to Sym²(V±)} admits feasible
points violating BOTH targets:

| p | PSD-feasible max|μ_int| | 2/n·N | L·N | true max |
|---|---|---|---|---|
| 5 | 31.2 (τ-window [−240/7, 66]) | 20 | 15.6 | 12 |
| 7 | 2466 | 916.16 | 1168.57 | 872 |

So no argument whose only inputs are 4-point moment consistency, symmetry,
and positivity can prove |μ|≤2/n or |μ|≤L.  This subsumes and explains the
earlier failures (|μ|≤|f4| false at p=7; Wick/hull/Ext slack; crude Σ|per|):
**do not reopen level-4 majorants for leftover 3.**  Degree 6 closes the
gap at p=5 but NOT at p=7 (kernel dim 4 — the four free directions touch
every 4-point coordinate).

## Isolated remaining estimate (replaces nothing, sharpens the route)

Let K_D(p) = equivariant kernel of the complete degree-≤D contraction
system.  Data: K_4(5)=1, K_4(7)=2, K_4(11)=4, K_6(5)=0, K_6(7)=4.

K_4 = #liveδ − 1 was verified at a FOURTH prime by prediction: the p=13
run gave 7 live δ-orbits and K_4(13)=6, exactly (#liveδ − 1).  Both
structural facts (kernel entirely in δ; {δ-dead} = {|κ|=1}) hold set-wise
at all four primes (`scripts/kernel_dim_degree4.py`, KP=<prime>).

**Verdict on degree escalation:** K_4 grows with p (1,2,4,6 ≈ linearly,
tracking the |κ|=3 orbit count), and degree 6 closes only p=5.  Chasing a
kernel-zero degree D(p) is a census in disguise; treat it as dead unless
someone proves the closing degree is bounded.

## The affine reduction, calibrated — signs of ν are load-bearing

Since the kernel is carried by δ, the degree-4 system determines μ on
|κ|=1 as an EXACT affine function of the δ-orbit values (= ν on the
|κ|=3 orbits; 15.268's ν vanishes on |κ|=1, is free on |κ|=3):

    μ|κ=1  =  a  +  C·ν|κ=3        (data-free a, C; exact at p=5,7:
                                     reconstruction matches enumerated
                                     moments per orbit, Fraction equality)

Calibration (`affine_calibrate.py`): the triangle inequality
|μ| ≤ |a| + ‖C‖₁·max|ν| closes NOTHING — at p=7 one |κ|=1 orbit has
|a| = 1402.3 > both thresholds (2/n: 916.2, L: 1168.6) and only signed
cancellation brings the true value to 488.  Worse, keeping the TRUE
magnitudes of ν and flipping signs breaks both bounds at both primes
(p=5: 77.3 vs 40/31.2; p=7: 2316.6 vs 916.2/1168.6).  So within this
(exact) reduction, **no magnitude-only estimate on ν can close leftover
3 — the signed values are required.**

## Isolated remaining estimate (the deliverable)

> **Leftover 3 reduces to: determine, with correct signs and error
> smaller than the calibrated slack, the values of ν = ½(m₄⁺−m₄⁻) on
> the |κ|=3 signed orbits of four-sets (2, 3, 5, 7 of them at
> p = 5, 7, 11, 13).  The data-free affine map then forces μ on |κ|=1,
> and |μ| ≤ L (Type-I close) or ≤ 2/n follows by evaluation.**

The |κ|=3 orbits are exactly the special cross-ratio classes
(15.267/15.268 territory: T-products all equal, harmonic/equianharmonic
structure); a closed form for ν there is a character-sum problem on a
slowly growing family, not a census of Max±.  This replaces "prove
|μ|≤2/n directly" as the sharpest known form of leftover 3.

1. Find the character-algebra formula for dim K_D(p) (15.589's PSL(2,q)
   machinery is the toolkit; the constraint spaces are explicit
   PΓL-modules).  Determine whether some D(p) has K_D(p)=0 for all p
   (D(5)=6; conjecture to test: D(p)=p+1, i.e. K_8(7)=0 — unverified).
2. If yes: μ has a symmetry-forced closed form for every p; leftover 3
   reduces to bounding that form (and leftover 1's principal-room
   E[(y·z)⁴]=Σ T⁺(S)² becomes computable from the same determination —
   shared blocker removed).
3. If no: the PSD windows must be re-examined at the closing degree; the
   level-6 odd-block (V₊⊕Sym³(V₊)) PSD scan at p=7 was NOT run (open).

Done: K_4(11)=4, K_4(13)=6 (data-free, above).  Next cheap experiments:
the |κ|=3 orbit-count formula (pure Burnside over PΓL(2,q)); ν values on
those orbits at p=5,7 as cross-ratio data to seed the character-sum hunt;
K_4(17) if the count formula needs another point.

## Tooling notes

- Max± enumeration by V±-basis completion (dim n/2 pivoted basis, 2^{n/2}
  completions): p=5 instant, p=7 in ~13 s — ~10³× faster than prior paths
  at p=7.  Wall at p=11 (2^61) unchanged, consistent with fable.md.
- Six-set moment sums via triple-Gram (X₃ᵀX₃, fp32-exact): 15.9M sets in
  ~40 s at p=7.
- Signed orbits by vectorized label propagation with sign cocycles and
  dead-orbit (sign-conflict) detection; anti-automorphisms (UᵀCU=−C)
  included for μ, twisted by s for δ.
