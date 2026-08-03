# Dual-gap / Hypothesis H — attack status (2026-08-03)

## Preferred target (Props 15.159–15.160)

\[
G=\frac{d}{32}(16I-\Phi|_Z)\succeq I
\quad\Longleftrightarrow\quad
\lambda_{\max}(\Phi)\le 16\frac{d-2}{d}
\quad\Longleftrightarrow\quad 16N.
\]

**Hypothesis H:** \(\mathrm{ray}_{\max}\le H(p)=(p+2)^2/d\).

| Fact | Status |
|------|--------|
| \(\mathrm{thr\_ray}-H=(3p+7)(p-5)/(2d)\) | **Proved** (15.160) |
| H ⇒ G≽I ⇒ 16N for all primes \(p\ge5\) | **Proved** (15.160) |
| H ≤ 5 (⇒ 16N if ray≤H) | **Proved** (15.63/15.160) |
| ray=H=thr_ray at p=5; ray<H<thray at p=7 | **Certified** |
| G eigs {1,2,4} at p=5 | **Certified** (structure clue) |
| H for general \(p\ge5\) | **OPEN** |
| residual δ²≤room_hyp/24 general | **OPEN** |
| L = lim α_n | **OPEN** |

## Attacks tried this arc (not closing)

1. Rational closed form for λ_max / orth — thrash; no simple poly fit for general p.
2. Max+ IP association scheme — dead (15.158).
3. G-Schurian BM on Max+: **124** pair orbitals at p=5 ≫ 3 Φ levels — Schurian commutant too large; Φ is special element with extra degeneracy.
4. Delsarte/Welch/Gershgorin on φ-frame — too weak for λ_max≤16.
5. Majorization with mult≥d without orth bound — too weak.
6. 2×sphere bound \(16d/(d+2)\): holds numerically p=5,7; **no proof** that Max+ Rayleigh ≤ 2× spherical/Gaussian degree-2 energy (would give 16N; equality pattern at p=3 with λ≡16).
7. G·hs character-sum proxy — dead (δ²_Ghs ≫ room, 15.135).
8. Weil on full Max+ — needs Max+ parametrization beyond O_hs (15.135.E).

## What would close 16N / bi-tight

Any one of:

- Prove H for all primes \(p\ge5\)
- Prove G≽I (dual gap) directly by BM/SDP dual/rep ID of G
- Prove orth≤room (or ED4≤wick_hi / κ²≤96n / δ²≤room_hyp/24)
- Prove Max+ degree-2 Rayleigh ≤ 2 × Gaussian (cov 2I on V+)

Then: 16N → bi-tight empty (15.61) for p≥5. Still need E(1)/deep ND for Main / L.

## Non-goals respected

No soft-close from sandwich. No class_key thrash. No Prop mill that only renames δ.
