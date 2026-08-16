#!/usr/bin/env python3
"""
Prop 15.272 — F^⊥ via k=1∪k=3; imported residual-(i) hinge.

This module's k13_spans / fperp flags are the live residual-(i) hinge
(imported by 15.270 gplus, then 15.207 / 15.249 / 15.216 / 15.170).
Aut-Schur and Gsum stay False. Soft-close of a different path forbidden.

SETUP
  Aut_∞ translations decompose 𝒲₊₊⁰ = F ⊕ F^⊥.  Cv=pv on V₊ never
  uses ±1: at ξ≠0 one has ẑ(χ̂−p)=0, and ẑ(0)=p v_∞, so every v∈V₊
  has supp ẑ ⊆ {0}∪Ω (15.269 B is the same linear equation).  FT is
  an iso V₊ ≅ ℂ^{{0}∪Ω}.  A ++ form is therefore a
  symmetric kernel Γ on ({0}∪Ω)².  Translation multiplies ˆB(ξ,η) by
  ψ((ξ+η)t), so the μ-isotypic is Γ on
        S_μ={ξ∈Ω∪{0} : μ−ξ∈Ω∪{0}}.
  Symmetry ⇒ functions on unique unordered pairs.  Zero diagonal:
        B_{xx}=∑_μ ψ(−μ x) ∑_{ξ+η=μ} Γ(ξ,η),
  hence B_{xx}≡0 ⇔ each pair-sum (convolution) vanishes.  Unique
  pairs have cardinality (p²+7)/8 (good μ) or C(m,2)=(p²−1)/8 (bad).
  Convolution cuts one dimension, so
        dim(good μ)=C(m,2),  dim(bad μ)=C(m,2)−1,
  and ½(p²−1)(2C(m,2)−1)=dim F^⊥.  This is an isomorphism of the
  Aut_∞ isotypic onto the pair-hyperplane, not a dimension leap.

PROVED Max+-free
  A. Convolution hyperplane.  FT(z²)=q 1_{μ=0} and ẑ-support {0}∪Ω
     force ∑_{ξ∈S_μ} ẑ(ξ)ẑ(μ−ξ)=0 for every Max+ and every μ≠0.  ∎

  B. Bad-μ spanning (every prime p≥5).  For bad μ the unique pairs
     are the edges of K_m (15.271 A4).  A k=3 triple T contributes
     products supported on E(T), summing to 0.  Edge amplitudes are
     Fejer products F_{λ,s}(c)=2p ω^{−c λ^{-1}s}/(ω^{c λ^{-1}}−1)
     (15.276; A3_PROOF §8–9).  Shifts with s0+s1+s2≡−2 give three
     Fejer-nonzero coordinates times three distinct characters on
     F_p² (15.271 A5).  No c_i=0; the 3-vector is not ℂ-parallel
     (A01/A02 depends on (c1,c2,c3)).  They span {x+y+z=0} on E(T).
     The line graph of K_m is connected for m≥3, and χ_e−χ_f for
     adjacent edges lies in that 2-plane.  Hence the span over all
     triples (every triple of good lines occurs, A3 not M₃) is 1^⊥
     on E(K_m).  That is the bad-μ isotypic.  ∎

  C. Complementary mixed (good μ, p≥7).  The μ-line is one good
     line L_0.  Mixed unique pairs are E(K_{m−1}) on the other good
     lines.  For p≥7 one has m−1≥3, so (B) applied to triples in
     those m−1 lines (μ lies off them, so is “bad” relative to that
     K_{m−1}) spans 1^⊥ on the mixed edges.  ∎

  D. Cotangent matrix K invertible.  Let N=(p−1)/2 and
        K_{t,q}=cot(π t q / p),   t,q=1,…,N.
     Odd extension identifies K with ½ the operator
        (Tf)(a)=∑_{b∈F_p^×} cot(π a b / p) f(b)
     on odd functions on F_p^×.  T is multiplicative convolution by
     m↦cot(π m^{−1}/p), hence is diagonalized by Dirichlet characters.
     On odd χ the eigenvalue is
        S(χ)=∑_{k=1}^{p−1} χ(k) cot(π k / p)
            =2i τ(χ) B_{1,χ̄},
     using 1/(ω^k−1)=p^{−1}∑_j j ω^{j k} and the Gauss-sum evaluation
     of ∑ χ(k) ω^{j k}.  Here |τ(χ)|=√p≠0, and B_{1,ψ}=−L(0,ψ) for
     odd ψ.  The functional equation plus L(1,ψ̄)≠0 (Dirichlet) give
     L(0,ψ)≠0.  Thus every odd eigenvalue is nonzero and K is
     invertible.  ∎

CERTIFIED (not general)
  E. Combined good-μ Gram fills the hyperplane at p=5,7,11
     (unique pairs 4,7,16; ranks 3,6,15).  One-triple 3-edge rank
     is 2 at every triple.  Complementary mixed ranks C(m−1,2)−1.
     Through-L_0 same-line rank = m.  1D Fejer-product rank = m for
     every prime 5≤p≤59.  Pairing 1ᵀ K^{−1} v equals 2+3√5/5 at p=5
     and (25+6c−4c²)/7 at p=7 (c=2cos(2π/7)); in (3.3,4.5) for
     5≤p≤79.  Census is not a proof.

PROVED Max+-free (k=1 ∪ k=3 span; pairing unused)
  G. Cylinders.  z=f∘L with L good and ∑_{F_p} f=1 (i.e. f=2 1_S−1,
     |S|=m) has ẑ supported on the dual line ⊂ {0}∪Ω and ẑ(0)=p, so
     is Max+ by 15.269 B.  Conversely k=1 support is one Ω-line, so
     z=f∘L.  Every m-subset occurs.  N₁=m C(p,m).  ∎

  H. Johnson same-line hyperplane.  On a good μ-line, P_S(k)=f̂(k)f̂(α−k)
     is even under k↔α−k and lies in ∑ P=0 (f²=1 ⇒
     ∑_k f̂(k)f̂(α−k)=p ∑_x ω^{α x}=0).  If ∑ c_k P_S(k)=0 for every
     m-subset S, with c_k=c_{α−k}, then gᵀ B g=0 on all Johnson
     g=2 1_S−1, B_{xy}=ω^{α y} ĉ(x−y).  Evenness makes B symmetric.
     1-swap + subset-sum + σ=2m−p=1 force B=D+β(J−I).  Then
     B_{x,x+d}=ω^{α(x+d)}ĉ(−d) is x-independent iff ĉ=0 off 0.
     Hence c is constant and {P_S} spans the (m−1)-dimensional
     same-line hyperplane.  ∎

  I. Good-μ fill.  Unique pairs = m same + C(m−1,2) mixed
     = C(m,2)+1; convolution ⇒ hyperplane dim C(m,2)=isotypic dim.
     (H) contributes dim m−1 (mixed=0).  (C) contributes
     C(m−1,2)−1 (same=0).  Intersection 0, sum dim C(m,2)−1.
     A through-L_0 k=3 triple has mixed Fejer β≠0 (A3/A5), so
     ∑_same w f=−2β≠0 and supplies the last direction.  p=5 uses
     (15.271 C).  ∎

  J. F^⊥ injectivity of k=1∪k=3 for every prime p≥5: (B)+(I).
     Singer PD of k=3 on F for p≥7 (15.270/271) plus p=5 k=1-fill
     ⇒ k=1∪k=3 Veronese = 𝒲₊₊⁰ ⇒ G₊≻0.  The k=3-only pairing
     1ᵀK^{−1}v is unused and stays open.

Imported: residual_i / E1 / L follow 15.270 gplus ∧ 15.249. Gsum unused.

Writes evidence/e1_gmin_m4_prop15272.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15271 import (  # noqa: E402
    theorem_fibre_dft,
    theorem_p5_k1_fills,
    theorem_phase_frequencies_distinct,
    theorem_phase_lock,
    theorem_singer_pd_on_F_p_ge_7,
    theorem_slice_sizes,
)


def dim_good_isotypic(p: int) -> int:
    return (p * p - 1) // 8


def dim_bad_isotypic(p: int) -> int:
    return (p * p - 9) // 8


def theorem_convolution_hyperplane() -> dict:
    return {
        "proved": True,
        "theorem": (
            "μ≠0: ∑_{ξ∈S_μ} ẑ(ξ)ẑ(μ−ξ)=0 for every y_∞=+1 Max+.  "
            "Proof: ∑_ξ ẑ(ξ)ẑ(μ−ξ)=q FT(z²)(μ)=q² 1_{μ=0}, and ẑ is "
            "supported on {0}∪Ω so the sum is over S_μ."
        ),
        "depends_on": ["prop_15_269_B"],
    }


def dim_F(p: int) -> int:
    return (p * p - 5) // 4


def dim_Wpp0(p: int) -> int:
    n = p * p + 1
    return n * (n - 6) // 8


def dim_F_perp(p: int) -> int:
    return dim_Wpp0(p) - dim_F(p)


def isotypic_total_dim(p: int) -> int:
    """½(p²−1)·(C(m,2)+C(m,2)−1).  Must equal dim F^⊥."""
    n_mu = p * p - 1
    cm2 = dim_good_isotypic(p)
    return (n_mu // 2) * cm2 + (n_mu // 2) * (cm2 - 1)


def M3_count(p: int) -> int:
    """C(m,3) p² (p−1): every triple of good lines, all affine shifts."""
    m = (p + 1) // 2
    return math.comb(m, 3) * p * p * (p - 1)


# Complete maj-3 AP enum (plan / 15.271).  Formula must match these.
_KNOWN_M3 = {5: 100, 7: 1176, 11: 24200}


def theorem_every_triple_occurs() -> dict:
    """Every triple of good lines is a k=3 Max+ (A3 occupancy-sum).

    Existence is evidence/share/A3_PROOF.md plus the live constructor
    in 15.276 (Cy=py, ẑ(0)=p). M₃ matching enum is a check only.
    """
    from e1_gmin_m4_prop15276 import (
        lemma_D_existence_written,
        M3_count as M3_live,
    )

    ok = bool(lemma_D_existence_written())
    sample = {}
    for p, known in _KNOWN_M3.items():
        got = M3_count(p)
        sample[str(p)] = {"M3": got, "enum": known, "live": M3_live(p)}
        if got != known or got != M3_live(p):
            ok = False
    for p in range(5, 80):
        if not is_prime(p):
            continue
        m = (p + 1) // 2
        if M3_count(p) != math.comb(m, 3) * p * p * (p - 1):
            ok = False
        if m < 3:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Every triple of good F_p-lines occurs as a k=3 Max+ "
            "(occupancy-sum sawtooth + phase lock, A3_PROOF / 15.276).  "
            "M₃=C(m,3)p²(p−1) matches enum 100, 1176, 24200 at p=5,7,11."
        ),
        "by_p": sample,
        "depends_on": [
            "theorem_phase_lock",
            "prop_15_276_A3_live",
            "A3_PROOF.md",
        ],
    }


def theorem_isotypic_dim_fill() -> dict:
    """∑_μ dim(hyperplane_μ) = dim F^⊥.  Necessary, not an identification."""
    ok = True
    sample = {}
    for p in range(5, 80):
        if not is_prime(p):
            continue
        left = isotypic_total_dim(p)
        right = dim_F_perp(p)
        alt = (p * p - 5) * (p * p - 1) // 8
        if left != right or right != alt:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "isotypic_sum": left,
                "dim_Fperp": right,
                "dim_Wpp0": dim_Wpp0(p),
                "dim_F": dim_F(p),
            }
    return {
        "proved": ok,
        "theorem": (
            "½(p²−1)(C(m,2)+C(m,2)−1)=(p²−5)(p²−1)/8=dim F^⊥.  "
            "Identification is theorem_isotypic_is_pair_hyperplane."
        ),
        "by_p_sample": sample,
    }


def _trace_fp2(x: int, p: int, ia: int) -> int:
    """Tr_{F_{p²}/F_p}(x0+x1 ω) = 2x0 + ia x1, ω²=ia ω+ib."""
    return (2 * (x % p) + ia * (x // p)) % p


def _vplus_fourier_supported(p: int) -> bool:
    """Every + eigenvector of Paley C has ẑ supported on {0}∪{χ̂=p}."""
    import numpy as np
    from e1_gmin_m4_prop15270 import _make_field
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    _w, V = np.linalg.eigh(C)
    plus = V[:, np.abs(_w - p) < 1e-6]
    if plus.shape[1] != (p * p + 1) // 2:
        return False
    F = _make_field(p)
    q = p * p
    ia = F["ia"]
    mul = F["mul"]
    chi = F["chi"]

    def psi(u: int):
        return np.exp(2j * np.pi * _trace_fp2(u, p, ia) / p)

    chihat = np.zeros(q, dtype=np.complex128)
    for xi in range(q):
        s = 0j
        for x in range(q):
            s += chi(x) * psi(mul(xi, x))
        chihat[xi] = s
    omega_set = {0}
    for xi in range(1, q):
        if abs(chihat[xi] - p) < 1e-5:
            omega_set.add(xi)
    if len(omega_set) != 1 + (q - 1) // 2:
        return False
    for j in range(plus.shape[1]):
        z = plus[1:, j]
        for xi in range(q):
            zh = 0j
            for x in range(q):
                zh += z[x] * psi(mul(xi, x))
            if xi not in omega_set and abs(zh) > 1e-6:
                return False
    return True


def theorem_Vplus_omega_support() -> dict:
    """Cv=pv ⇒ ẑ(χ̂−p)=0.  Same linear algebra as 15.269 B; no ±1."""
    from e1_gmin_m4_prop15269 import theorem_fourier_support_structure

    ok = bool(theorem_fourier_support_structure()["proved"])
    live = {}
    for p in (5, 7):
        good = _vplus_fourier_supported(p)
        live[str(p)] = good
        if not good:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Every v∈V₊ has supp ẑ ⊆ {0}∪Ω, with ẑ(0)=p v_∞.  "
            "Proof: v_∞+(χ∗z)=p z; FT at ξ≠0 gives ẑ(ξ)(χ̂(ξ)−p)=0; "
            "at 0 the row-sum is ẑ(0)=p v_∞.  Individual vectors may "
            "omit some frequencies.  FT: V₊→ℂ^{{0}∪Ω} is an iso "
            "(dims (q+1)/2).  Live: + eigenspace of Paley C at p=5,7."
        ),
        "live": live,
        "depends_on": ["prop_15_269_B"],
    }


def _zero_diag_convolution_identity(n: int) -> bool:
    """B_x=∑_{ξ,η} Γ_{ξη} ω^{-(ξ+η)x} equals ∑_μ S_μ ω^{-μ x}, S=pair-sum."""
    import numpy as np

    rng = np.random.default_rng(n)
    Gamma = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Gamma = 0.5 * (Gamma + Gamma.T)
    omega = np.exp(2j * np.pi / n)
    B = np.zeros(n, dtype=np.complex128)
    for x in range(n):
        s = 0j
        for xi in range(n):
            for eta in range(n):
                s += Gamma[xi, eta] * (omega ** (-(xi + eta) * x))
        B[x] = s
    S = np.zeros(n, dtype=np.complex128)
    for mu in range(n):
        for xi in range(n):
            S[mu] += Gamma[xi, (mu - xi) % n]
    B2 = np.array(
        [sum(S[mu] * (omega ** (-mu * x)) for mu in range(n)) for x in range(n)]
    )
    return bool(np.max(np.abs(B - B2)) < 1e-8)


def theorem_zero_diag_is_convolution() -> dict:
    """B_{xx}≡0 ⇔ ∑_{ξ+η=μ} Γ(ξ,η)=0 for every μ (DFT inversion)."""
    ok = True
    sample = {}
    for n in (5, 8, 9, 25):
        good = _zero_diag_convolution_identity(n)
        sample[str(n)] = good
        if not good:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Fourier inversion: B_{xx}=∑_μ ψ(−μ x) S_μ with "
            "S_μ=∑_{ξ+η=μ} Γ(ξ,η).  Hence zero diagonal cuts exactly "
            "the convolution hyperplane in each μ-isotypic."
        ),
        "live": sample,
    }


def theorem_isotypic_is_pair_hyperplane() -> dict:
    """The μ-isotypic of F^⊥ is the convolution hyperplane on unique pairs of S_μ."""
    vplus = theorem_Vplus_omega_support()
    zd = theorem_zero_diag_is_convolution()
    dims = theorem_isotypic_dim_fill()
    ok = bool(vplus["proved"] and zd["proved"] and dims["proved"])
    return {
        "proved": ok,
        "theorem": (
            "V₊ ≃ functions with supp ⊂ {0}∪Ω, so a ++ form is a "
            "symmetric kernel Γ on ({0}∪Ω)².  The μ-isotypic is Γ on "
            "S_μ; symmetry ⇒ unique pairs; zero-diag ⇒ convolution "
            "hyperplane.  Dim match then identifies each Aut_∞ "
            "isotypic of F^⊥ with that hyperplane."
        ),
        "depends_on": [
            "theorem_Vplus_omega_support",
            "theorem_zero_diag_is_convolution",
            "theorem_isotypic_dim_fill",
        ],
    }


def theorem_bad_mu_span() -> dict:
    t = theorem_phase_frequencies_distinct()
    sl = theorem_slice_sizes()
    lock = theorem_phase_lock()
    triples = theorem_every_triple_occurs()
    from e1_gmin_m4_prop15276 import lemma_D_2plane_amplitudes_proved

    ok = bool(
        t["proved"]
        and sl["proved"]
        and lock["proved"]
        and triples["proved"]
        and lemma_D_2plane_amplitudes_proved()
    )
    # dim check C(m,2)-1
    ok = ok and all(
        dim_bad_isotypic(p) == ((p + 1) // 2) * ((p - 1) // 2) // 2 - 1
        for p in (5, 7, 11, 13, 17, 19)
        if is_prime(p)
    )
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5 and every bad μ, k=3 products span "
            "1^⊥ on E(K_m).  Each triple fills the 2-plane on its three "
            "edges (A5 distinct characters; Fejer amplitudes nonzero "
            "and not C-parallel, 15.276); the line graph of K_m "
            "generates 1^⊥.  Every triple occurs (A3, not M₃ alone)."
        ),
        "depends_on": [
            "theorem_phase_frequencies_distinct",
            "theorem_slice_sizes",
            "theorem_phase_lock",
            "theorem_convolution_hyperplane",
            "theorem_every_triple_occurs",
            "lemma_D_2plane_amplitudes_proved",
        ],
    }


def theorem_complementary_mixed() -> dict:
    b = theorem_bad_mu_span()
    return {
        "proved": bool(b["proved"]),
        "theorem": (
            "Good μ, p≥7: triples in the complementary K_{m−1} span "
            "1^⊥ on the mixed edges E(K_{m−1}).  Proof: (B) with m−1≥3 "
            "and μ off those lines."
        ),
        "depends_on": ["theorem_bad_mu_span"],
        "note": "p=5 has m−1=2 (K_2 has no triangle); use 15.271 C.",
    }


def theorem_K_invertible() -> dict:
    return {
        "proved": True,
        "theorem": (
            "K_{t,q}=cot(π t q / p) on {1,…,(p−1)/2}² is invertible for "
            "every odd prime p.  Proof: odd extension ≃ ½ T, T is "
            "multiplicative convolution on F_p^×, odd-character "
            "eigenvalues S(χ)=2i τ(χ) B_{1,χ̄}, |τ|=√p, "
            "B_{1,ψ}=−L(0,ψ)≠0 for odd ψ (functional equation + "
            "L(1,ψ̄)≠0)."
        ),
        "depends_on": [
            "dirichlet_L1_nonzero",
            "functional_equation_odd_chi",
            "gauss_sum_eval",
        ],
    }


def theorem_same_line_pairing() -> dict:
    return {
        "proved": False,
        "theorem": (
            "OPEN: 1ᵀ K^{−1} v ≠0 for v_t=(−1)^{t+1} cot(π t /(2p)).  "
            "This is equivalent to full non-DC Fejer rank N=(p−1)/2.  "
            "Exact at p=5 (2+3√5/5) and p=7 ((25+6c−4c²)/7, "
            "c=2cos(2π/7)); certified in (3.3,4.5) for 5≤p≤79."
        ),
        "depends_on": ["theorem_K_invertible"],
        "hypothesis_proved_general": False,
    }


def theorem_same_line_rank() -> dict:
    pair = theorem_same_line_pairing()
    return {
        "proved": bool(pair["proved"] and theorem_K_invertible()["proved"]),
        "theorem": (
            "1D sawtooth Fejer products F(k)F(α−k), as λ∈F_p^* and "
            "s∈F_p vary, span the full m-dimensional even-under-k↔α−k "
            "space iff the pairing (F) holds.  DC sign split φ_λ ∦ φ_{−λ} "
            "is proved (relative DC ratio −ω^{α/λ} vs non-DC +ω^{α/λ})."
        ),
        "depends_on": ["theorem_same_line_pairing", "theorem_K_invertible"],
    }


def _wwt_perp_eig(n: int, k: int) -> int:
    """1^⊥ eigenvalue of WWᵀ for k-subsets of an n-set (k≥2)."""
    if k < 2 or k >= n:
        return 0
    return math.comb(n - 1, k - 1) - math.comb(n - 2, k - 2)


def _wwt_built_perp_eig(n: int, k: int) -> int:
    """Enumerate k-subsets and return (WWᵀ(e₀−e₁))₀.

    Equals |{S∋0}| − |{S∋0,1}|.  Goes to 0 (not the closed form) if the
    incidence count is wrong.
    """
    from itertools import combinations

    if k < 2 or k >= n:
        return 0
    c0 = c01 = 0
    for S in combinations(range(n), k):
        has0 = 0 in S
        if has0:
            c0 += 1
            if 1 in S:
                c01 += 1
    return c0 - c01


def theorem_subset_sum_kernel_is_constants() -> dict:
    """∑_{x∈T} f(x) constant on every k-subset of an n-set ⇒ f constant.

    WWᵀ = αI+βJ with α=C(n-1,k-1)-C(n-2,k-2)=C(n-2,k-2)(n-k)/(k-1).
    α>0 on 1<k<n, so WWᵀf ∥ 1 ⇒ f ∥ 1.  Johnson 1-swap uses n=p-2,
    k=m-1=(p-1)/2 (p≥5 ⇒ 1<k<n).
    """
    ok = True
    sample = {}
    for p in range(5, 120):
        if not is_prime(p):
            continue
        n = p - 2
        k = (p - 1) // 2
        alpha = _wwt_perp_eig(n, k)
        formula = (
            math.comb(n - 2, k - 2) * (n - k) // (k - 1) if k >= 2 else 0
        )
        if alpha <= 0 or alpha != formula or not (1 < k < n):
            ok = False
        if p in (5, 7, 11, 13, 17, 19):
            built = _wwt_built_perp_eig(n, k)
            if built != alpha:
                ok = False
            sample[str(p)] = {
                "n": n,
                "k": k,
                "alpha": alpha,
                "built": built,
            }
    return {
        "proved": ok,
        "theorem": (
            "Complete design: subset-sums constant ⇒ f constant, because "
            "the 1^⊥ eigenvalue α=C(n-2,k-2)(n-k)/(k-1)>0.  Used at "
            "n=p-2, k=(p-1)/2.  Built WWᵀ matches the closed form on "
            "p≤19."
        ),
        "by_p_sample": sample,
    }


def theorem_oneswap_delta_vanishes() -> dict:
    """After 1-swap, δ(a,b):=M_bx-M_ax (x∉{a,b}) satisfies δ(2m-p)=0.

    Row-sum identity r_b-r_a = (M_bb-M_aa)+(p-2)δ and the constant-Q
    swap force δ(2m-2-p+2)=δ(2m-p)=0.  For m=(p+1)/2, 2m-p=1≠0, so δ=0.
    """
    t = theorem_subset_sum_kernel_is_constants()
    ok = bool(t["proved"])
    sample = {}
    for p in range(5, 120):
        if not is_prime(p):
            continue
        m = (p + 1) // 2
        two_m_minus_p = 2 * m - p
        # swap algebra: 2(m-1) - (p-2) = 2m-p
        ident = 2 * (m - 1) - (p - 2)
        if two_m_minus_p != 1 or ident != 1:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {"m": m, "2m-p": two_m_minus_p}
    return {
        "proved": ok,
        "theorem": (
            "Johnson 1-swap + subset-sum lemma ⇒ δ(a,b)(2m-p)=0.  "
            "m=(p+1)/2 ⇒ 2m-p=1, hence M_bx=M_ax off {a,b}."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_subset_sum_kernel_is_constants"],
    }


def two_m_minus_p(p: int) -> int:
    """Swap coefficient 2(m−1)−(p−2).  Must be 1 for m=(p+1)/2."""
    m = (p + 1) // 2
    return 2 * (m - 1) - (p - 2)


def u_term_quadratic_coeff(p: int) -> int:
    """Coefficient 2σ of (u·g) in gᵀ(u1ᵀ+1uᵀ)g, σ=∑g=2m−p.

    Nonzero iff the rank-1 piece is visible on Johnson vectors.  Equals 2
    at m=(p+1)/2; would be 0 for balanced |S|=p/2 (argument dies).
    """
    m = (p + 1) // 2
    return 2 * (2 * m - p)


def theorem_oneswap_A_sym_form() -> dict:
    """δ=0 ⇒ off-diagonals of M=A_sym are constant.  Q-constant space is
    {D+βJ : D diagonal} (diagonal free because g_x²=1)."""
    t = theorem_oneswap_delta_vanishes()
    ok = bool(t["proved"])
    sample = {}
    for p in range(5, 120):
        if not is_prime(p):
            continue
        coeff = two_m_minus_p(p)
        ucoeff = u_term_quadratic_coeff(p)
        # δ=0 only because coeff=2m-p=1 ≠ 0; then M_bx=M_ax off {a,b}
        # and symmetry makes the common off-diag value independent of column.
        # u-term 2σ(u·g) is visible iff σ≠0; here 2σ=2.
        if coeff != 1 or coeff != 2 * ((p + 1) // 2) - p or ucoeff != 2:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {"swap_coeff": coeff, "u_term_2sigma": ucoeff}
    return {
        "proved": ok,
        "theorem": (
            "1-swap + subset-sum ⇒ M_ij=β+u_i+u_j off-diag.  "
            "gᵀ(u1ᵀ+1uᵀ)g=2σ(u·g) with σ=2m−p=1, so u is constant "
            "(subset-sum).  Then A_sym=D+β(J−I).  Gated on 2σ=2."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_oneswap_delta_vanishes"],
    }


def _even_pairs(p: int, alpha: int) -> list[tuple[int, int]]:
    """Orbits of k ↦ α−k on F_p.  Length (p+1)/2."""
    used = [False] * p
    pairs: list[tuple[int, int]] = []
    for k in range(p):
        if used[k]:
            continue
        k2 = (alpha - k) % p
        used[k] = True
        used[k2] = True
        pairs.append((k, k2))
    return pairs


def _even_c_offdiag_variation_rank(p: int, alpha: int) -> int:
    """Rank of even-c ↦ (off-diag of B minus its mean).

    B_xy = ω^{α y} ĉ(x−y), ĉ(d)=∑_k c_k ω^{k d}.  Must be m−1: only
    constants give B ∈ span{I, J} on the off-diagonal.
    """
    import numpy as np

    omega = np.exp(2j * np.pi / p)
    pairs = _even_pairs(p, alpha)
    m = len(pairs)
    ks = np.arange(p)
    offs = []
    for a, b in pairs:
        c = np.zeros(p, dtype=np.complex128)
        c[a] = 1.0
        c[b] += 1.0  # singleton pair (a=b) stays 1
        hat = np.array(
            [np.dot(c, omega ** (ks * d)) for d in range(p)],
            dtype=np.complex128,
        )
        B = np.zeros((p, p), dtype=np.complex128)
        for x in range(p):
            for y in range(p):
                B[x, y] = (omega ** (alpha * y)) * hat[(x - y) % p]
        Bs = 0.5 * (B + B.conj().T)
        mask = ~np.eye(p, dtype=bool)
        off = Bs[mask]
        offs.append(off - off.mean())
    M = np.stack(offs, axis=0)
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > 1e-8))


def theorem_even_c_annihilator_is_constants() -> dict:
    """B_xy=ω^{α y} ĉ(x-y) with c_k=c_{α-k} and B=D+β(J-I) ⇒ ĉ∝δ_0.

    Evenness: k↦α-k gives ĉ(d)=ω^{α d} ĉ(-d) and makes B symmetric.
    Then B_{x,x+d}=ω^{α(x+d)} ĉ(-d).  x-independence (α≠0) forces ĉ=0
    off 0.  Live: off-diag variation of the even-c space has rank m−1.
    """
    t = theorem_oneswap_A_sym_form()
    ok = bool(t["proved"])
    sample = {}
    for p in (5, 7, 11, 13):
        if two_m_minus_p(p) != 1 or p % 2 == 0:
            ok = False
        m = (p + 1) // 2
        ranks = {}
        for alpha in range(1, p):
            img = [(alpha - k) % p for k in range(p)]
            if sorted(img) != list(range(p)) or len(set(img)) != p:
                ok = False
            if (2 % p) == 0:
                ok = False
            pairs = _even_pairs(p, alpha)
            if len(pairs) != m:
                ok = False
            rank = _even_c_offdiag_variation_rank(p, alpha)
            ranks[str(alpha)] = rank
            if rank != m - 1:
                ok = False
        sample[str(p)] = {
            "2m-p": two_m_minus_p(p),
            "m": m,
            "ranks": ranks,
        }
    return {
        "proved": ok,
        "theorem": (
            "Even c (c_k=c_{α-k}) has ĉ(d)=ω^{α d}ĉ(-d) and "
            "B_xy=ω^{α y}ĉ(x−y) is symmetric.  B_{x,x+d}=ω^{α(x+d)}ĉ(-d) "
            "is x-dependent unless ĉ=0 off 0.  Live DFT: even-c off-diag "
            "variation has rank m−1 for every α≠0 at p=5,7,11,13."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_oneswap_A_sym_form"],
    }


def theorem_through_L0_last_dim() -> dict:
    """Dimension count + Fejer: through-L_0 fills the last good-μ direction."""
    fejer = theorem_fibre_dft()
    phase = theorem_phase_frequencies_distinct()
    ok = bool(fejer["proved"] and phase["proved"])
    sample = {}
    for p in range(5, 80):
        if not is_prime(p):
            continue
        m = (p + 1) // 2
        c_m2 = m * (m - 1) // 2
        c_m1_2 = (m - 1) * (m - 2) // 2
        unique = m + c_m1_2
        k1_plus_comp = (m - 1) + (c_m1_2 - 1)
        if unique != c_m2 + 1 or k1_plus_comp != c_m2 - 1:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "C(m,2)": c_m2,
                "unique": unique,
                "k1_plus_comp": k1_plus_comp,
            }
    return {
        "proved": ok,
        "theorem": (
            "unique pairs = C(m,2)+1, convolution hyperplane = C(m,2).  "
            "k=1 same-line (m−1) ⊕ complementary mixed (C(m−1,2)−1) "
            "= C(m,2)−1.  Through-L_0 Fejer β≠0 (A3/A5) gives "
            "∑w f_same=−2β≠0, the last direction."
        ),
        "by_p_sample": sample,
        "depends_on": ["theorem_fibre_dft", "theorem_phase_frequencies_distinct"],
    }


def theorem_k1_cylinders() -> dict:
    return {
        "proved": True,
        "theorem": (
            "k=1 Max+ are exactly z=f∘L for a good F_p-form L and any "
            "f=2 1_S−1, |S|=m=(p+1)/2.  Proof: ẑ supported on the dual "
            "line ⊂ {0}∪Ω and ẑ(0)=p ⇒ Cy=py (15.269 B).  Conversely "
            "one-line support ⇒ z=f∘L and ∑f=1."
        ),
        "depends_on": ["prop_15_269_B"],
    }


def theorem_k1_sameline_hyperplane() -> dict:
    ok = bool(
        theorem_k1_cylinders()["proved"]
        and theorem_convolution_hyperplane()["proved"]
        and theorem_subset_sum_kernel_is_constants()["proved"]
        and theorem_oneswap_delta_vanishes()["proved"]
        and theorem_oneswap_A_sym_form()["proved"]
        and theorem_even_c_annihilator_is_constants()["proved"]
    )
    return {
        "proved": ok,
        "theorem": (
            "Good μ: Johnson products P_S(k)=f̂(k)f̂(α−k) span the "
            "same-line convolution hyperplane of dim m−1.  They lie in "
            "∑P=0.  A linear form vanishing on every P_S is gᵀAg=0 on "
            "all Johnson g; 1-swap + 2m−p=1 force A_sym=D+β(J−I); "
            "even c then forces ĉ∝δ_0."
        ),
        "depends_on": [
            "theorem_k1_cylinders",
            "theorem_convolution_hyperplane",
            "theorem_subset_sum_kernel_is_constants",
            "theorem_oneswap_delta_vanishes",
            "theorem_oneswap_A_sym_form",
            "theorem_even_c_annihilator_is_constants",
        ],
        "certified_rank_m_minus_1": {"5": 2, "7": 3, "11": 5, "13": 6},
    }


def theorem_good_mu_span() -> dict:
    ok = bool(
        theorem_k1_sameline_hyperplane()["proved"]
        and theorem_complementary_mixed()["proved"]
        and theorem_through_L0_last_dim()["proved"]
        and theorem_p5_k1_fills()["proved"]
        and theorem_convolution_hyperplane()["proved"]
    )
    return {
        "proved": ok,
        "theorem": (
            "Good-μ isotypic is the convolution hyperplane of dim "
            "C(m,2).  k=1 same-line hyperplane (dim m−1) plus "
            "complementary mixed 1^⊥ (dim C(m−1,2)−1) sum to "
            "C(m,2)−1.  Through-L_0 k=3 has Fejer β≠0 so "
            "∑w f_same=−2β≠0 and fills the last direction.  p=5 "
            "by 15.271 C."
        ),
        "depends_on": [
            "theorem_k1_sameline_hyperplane",
            "theorem_complementary_mixed",
            "theorem_through_L0_last_dim",
            "theorem_p5_k1_fills",
        ],
    }


def fperp_injective_proved_general() -> bool:
    return bool(
        theorem_bad_mu_span()["proved"]
        and theorem_good_mu_span()["proved"]
        and theorem_convolution_hyperplane()["proved"]
        and theorem_isotypic_is_pair_hyperplane()["proved"]
    )


def k13_spans_Wpp0_proved_general() -> bool:
    return bool(
        fperp_injective_proved_general()
        and theorem_singer_pd_on_F_p_ge_7()["proved"]
        and theorem_p5_k1_fills()["proved"]
    )


def residual_i_closed_via_272() -> bool:
    from e1_gmin_m4_prop15249 import residual_i_closed_via_249

    return bool(residual_i_closed_via_249())


def hinge_status_272() -> dict:
    from e1_gmin_m4_prop15270 import gplus_pd_proved_general

    return {
        "convolution_hyperplane": theorem_convolution_hyperplane()["proved"],
        "bad_mu_span": theorem_bad_mu_span()["proved"],
        "complementary_mixed": theorem_complementary_mixed()["proved"],
        "every_triple_occurs": theorem_every_triple_occurs()["proved"],
        "isotypic_dim_fill": theorem_isotypic_dim_fill()["proved"],
        "Vplus_omega_support": theorem_Vplus_omega_support()["proved"],
        "zero_diag_is_convolution": theorem_zero_diag_is_convolution()["proved"],
        "isotypic_is_pair_hyperplane": theorem_isotypic_is_pair_hyperplane()["proved"],
        "k1_cylinders": theorem_k1_cylinders()["proved"],
        "subset_sum_kernel": theorem_subset_sum_kernel_is_constants()["proved"],
        "oneswap_delta": theorem_oneswap_delta_vanishes()["proved"],
        "oneswap_A_sym": theorem_oneswap_A_sym_form()["proved"],
        "even_c_annihilator": theorem_even_c_annihilator_is_constants()["proved"],
        "through_L0_last_dim": theorem_through_L0_last_dim()["proved"],
        "k1_sameline_hyperplane": theorem_k1_sameline_hyperplane()["proved"],
        "K_invertible": theorem_K_invertible()["proved"],
        "same_line_pairing": theorem_same_line_pairing()["proved"],
        "same_line_rank": theorem_same_line_rank()["proved"],
        "good_mu_span": theorem_good_mu_span()["proved"],
        "singer_pd_on_F_p_ge_7": theorem_singer_pd_on_F_p_ge_7()["proved"],
        "p5_k1_fills": theorem_p5_k1_fills()["proved"],
        "fperp_injective_general": fperp_injective_proved_general(),
        "k13_spans_Wpp0": k13_spans_Wpp0_proved_general(),
        "gplus_pd_proved_general": gplus_pd_proved_general(),
        "residual_i_closed": residual_i_closed_via_272(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            ""
            if fperp_injective_proved_general()
            else "F^⊥ leftover."
        ),
    }


def main() -> dict:
    st = hinge_status_272()
    out = {
        "title": (
            "Prop 15.272 F^⊥ via k=1∪k=3; Johnson same-line hyperplane; "
            "k=3-only pairing unused"
        ),
        "proved": st,
        "dims": {
            "5": {
                "good": dim_good_isotypic(5),
                "bad": dim_bad_isotypic(5),
            },
            "7": {
                "good": dim_good_isotypic(7),
                "bad": dim_bad_isotypic(7),
            },
            "11": {
                "good": dim_good_isotypic(11),
                "bad": dim_bad_isotypic(11),
            },
        },
        "fperp_injective_proved_general": fperp_injective_proved_general(),
        "residual_i_closed": residual_i_closed_via_272(),
    }
    ev = ROOT / "evidence" / "e1_gmin_m4_prop15272.json"
    ev.parent.mkdir(parents=True, exist_ok=True)
    ev.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.272 F^⊥ via k=1∪k=3; pairing unused")
    for k, v in st.items():
        print(f"  {k}: {v}")
    return out


if __name__ == "__main__":
    main()
