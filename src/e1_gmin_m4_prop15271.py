#!/usr/bin/env python3
"""
Prop 15.271 — k=3 Veronese: ker ⊆ Aut_∞-circulants F; Singer kills
F for p≥7; residual (i) still OPEN (F^⊥ injectivity not proved general).

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP
  G₊≻0 on 𝒲₊₊⁰ ⇔ span{yyᵀ−S : y∈Max+} = 𝒲₊₊⁰ (15.212).
  k=3 Max+ are the explicit 3-line sawtooths z=maj(s₁,s₂,s₃) (A3_PROOF),
  M₃=C((p+1)/2,3) p²(p−1), Cy=py exact.
  F := Aut_∞-circulants in 𝒲₊₊⁰
     = {B_{∞,x}=0, B_{x,y}=f(x−y), f even on Ω, mean 0},
     dim F=(p²−5)/4.  This is the Singer / even-on-Ω space of 15.270.

PROVED Max+-free
  A. Majority identity.  For sᵢ=±1,
        maj(s₁,s₂,s₃)=(s₁+s₂+s₃−s₁s₂s₃)/2
     = A + s₃ B with A=½(s₁+s₂), B=½(1−s₁s₂).  ∎

  A2. Phase lock.  If α∈ker(L₁,L₂,L₃) (the Plücker relation) and
     Nᵢ=1+(αᵢ c+sᵢ) mod p, then N₁+N₂+N₃∈{p+1,2p+1}≡1 (mod p)
     while α₁L₁+α₂L₂+α₃L₃≡0 forces
        s₀+s₁+s₂ ≡ −2 (mod p).  ∎

  A3. Fibre-sum DFT.  For ξ on the dual of a good form Lᵢ, ẑ(ξ) equals
     the 1-variable DFT of the occupancy fibre-sum Gᵢ=2Nᵢ−p.  A3 says
     Gᵢ is a sawtooth, so Fejer: ẑ(ξ)≠0 at every nonzero frequency on
     that line.  ∎

  A4. Additive slices.  For μ≠0,
        |S_μ| = (p²+3)/4 if F_p·μ⊂Ω∪{0} (good),
        |S_μ| = (p²−1)/4 if F_p·μ∩Ω=∅ (bad).
     Bad S_μ is in bijection with ordered pairs of distinct good lines
     (no same-line pairs, so no odd-odd same-line sector).  Good S_μ
     is the whole good line through μ (p points) plus mixed pairs among
     the other m−1 lines.  ∎

  A5. Phase frequencies on a triangle.  Shift theorem plus (A2) gives
     edge-product phase exponents (c₁,c₂), (c₁−c₃,−c₃), (−c₃,c₂−c₃)
     in F_p².  Any equality of two of these forces some cᵢ=0, which
     Fejer forbids.  So the three characters on F_p² are distinct.  ∎

  B. Singer / Fejer / Weil on F (15.270).  The locked-triple circulant
     on the m=(p+1)/2 good lines is PD for every prime p≥7:
     μ=0 block eigenvalues 3(m−1)(m−2)/2 and (m−2)(m−3)/2, both >0
     iff m≥4; Fejer invertibility of even stretches; Weil |T|≤2√p+1
     + Gershgorin p≥37 + DFT 7≤p≤31.  Hence G_{k=3}≻0 on F for p≥7.
     At p=5, m=3, the μ=0 off-diagonal eigenvalue (m−2)(m−3)/2=0.  ∎

  C. p=5 fill.  Only k=1 and k=3 exist.  Veronese of k=1∪k=3 has rank
     65=dim 𝒲₊₊⁰ (15.212).  k=3 alone is 61/65: the 4-dim kernel is
     the singular Singer subspace of F (Type C square-line indicators
     ⊕ Type D χ₅-pullbacks ⊥ f_hit).  Detected by k=1 cylinders.  ∎

CERTIFIED (not general)
  D. ker(G_{k=3}) ⊆ F, equivalently k=3 Veronese injective on F^⊥:
        p=5:  rank on F = 1/5,  on F^⊥ = 60/60,  total 61/65
        p=7:  rank on F = 11/11, on F^⊥ = 264/264, total 275/275
        p=11: rank on F = 29/29, on F^⊥ = 1740/1740, total 1769/1769
     GPU SVD of the explicit A3 family (Cy−py=0).  Per-triple rank is
     constant (135 at p=7, 365 at p=11).  Do **not** treat D as a
     general proof: Fejer-even products miss odd-odd Fourier sectors
     unless a full spanning lemma is proved (referee BLOCK 2026-08-15).

OPEN (blocks residual i / predicate flip)
  E. Prove F^⊥-injectivity for every prime p≥5 (every nontrivial
     additive isotypic of Aut_∞ on 𝒲₊₊⁰ pairs nondegenerately with
     {ẑ⊗ẑ : z k=3}, including odd-odd sectors and the ∞-row).  Then
     ker ⊆ F, and (B) gives ker=0 for p≥7, so G₊≻0; p=5 uses (C).
     Then 15.207 ker=sc + 15.249 cost_D empties dual-eq.

Until E: this module's residual_i_via_271 stays False (k=3-only).
The k=1∪k=3 close lives in 15.272 and does not flip this flag.

Writes evidence/e1_gmin_m4_prop15271.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15270 import (  # noqa: E402
    certify_small_dft,
    gplus_pd_proved_general,
    theorem_fejer_stretch,
    theorem_gershgorin_large_p,
    theorem_mu0_block,
    theorem_weil_T,
)


def dim_F(p: int) -> int:
    return (p * p - 5) // 4


def dim_Wpp0(p: int) -> int:
    n = p * p + 1
    return n * (n - 6) // 8


def theorem_majority_identity() -> dict:
    return {
        "proved": True,
        "theorem": (
            "maj(s1,s2,s3)=(s1+s2+s3-s1 s2 s3)/2 = A+s3 B "
            "for si=±1, A=½(s1+s2), B=½(1-s1 s2)."
        ),
    }


def theorem_phase_lock() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Plücker α∈ker(L1,L2,L3) and 2-valued occupancy sum force "
            "s0+s1+s2≡-2 (mod p).  Proof: α·L≡0 ⇒ ∑(α_i L_i + s_i) is "
            "the constant ∑s_i; N_i=1+(α_i L_i+s_i) mod p so "
            "∑N≡3+∑s (mod p); ∑N∈{p+1,2p+1}≡1."
        ),
    }


def theorem_fibre_dft() -> dict:
    return {
        "proved": True,
        "theorem": (
            "If ξ lies on the dual of a good form L, then ẑ(ξ) is the "
            "DFT of the fibre-sum G=2N-p on that form.  A3: G is a "
            "sawtooth, so Fejer ẑ(ξ)≠0 at every nonzero line-frequency."
        ),
    }


def theorem_slice_sizes() -> dict:
    ok = True
    sample = {}
    for p, g, b in ((5, 7, 6), (7, 13, 12), (11, 31, 30)):
        good = (p * p + 3) // 4
        bad = (p * p - 1) // 4
        sample[str(p)] = {"good": good, "bad": bad}
        if good != g or bad != b:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "μ≠0: |S_μ|=(p²+3)/4 on a good line, (p²−1)/4 on a bad line. "
            "Bad S_μ ↔ ordered pairs of distinct good lines (no same-line "
            "pairs).  Proof: N_Ω(μ)=(p²-3-2σ χ(μ))/4 by expanding "
            "1_Ω=(1+σ χ)/2 and ∑χ(x(μ-x))=-1."
        ),
        "by_p": sample,
    }


def theorem_phase_frequencies_distinct() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On a triple, edge-product phase vectors in F_p² are "
            "(c1,c2), (c1-c3,-c3), (-c3,c2-c3).  Any two equal ⇒ some "
            "ci=0, contradicting Fejer.  Hence the three characters "
            "a,b↦ω^{linear} are pairwise distinct."
        ),
    }


def theorem_singer_pd_on_F_p_ge_7() -> dict:
    ok = bool(
        theorem_fejer_stretch()["proved"]
        and theorem_mu0_block()["proved"]
        and theorem_weil_T()["proved"]
        and theorem_gershgorin_large_p()["proved"]
        and certify_small_dft()["certified"]
    )
    return {
        "proved": ok,
        "theorem": (
            "G_{k=3}≻0 on the Aut_∞-circulant space F for every prime p≥7, "
            "by 15.270 Fejer + μ=0 + Weil + Gershgorin/DFT.  "
            "At p=5 the μ=0 block is singular ((m-2)(m-3)/2=0)."
        ),
        "depends_on": ["prop_15_270"],
    }


def theorem_p5_k1_fills() -> dict:
    return {
        "proved": True,
        "theorem": (
            "At p=5 only k=1,3 exist.  Veronese(k=1∪k=3) has rank 65 "
            "(15.212).  k=3 kernel is 4-dimensional inside F "
            "(Type C ⊕ Type D), filled by k=1."
        ),
        "depends_on": ["prop_15_212"],
    }


def certify_fperp_injective() -> dict:
    """Census only. Not a general proof."""
    return {
        "certified": True,
        "proved_general": False,
        "by_p": {
            "5": {
                "dim_F": 5,
                "dim_Wpp0": 65,
                "rank_on_F": 1,
                "rank_on_Fperp": 60,
                "rank_total": 61,
                "ker_in_F": True,
            },
            "7": {
                "dim_F": 11,
                "dim_Wpp0": 275,
                "rank_on_F": 11,
                "rank_on_Fperp": 264,
                "rank_total": 275,
                "ker_in_F": True,
            },
            "11": {
                "dim_F": 29,
                "dim_Wpp0": 1769,
                "rank_on_F": 29,
                "rank_on_Fperp": 1740,
                "rank_total": 1769,
                "ker_in_F": True,
            },
        },
        "note": (
            "GPU SVD of explicit A3 k=3 family.  F^⊥ injectivity general "
            "is the leftover (odd-odd / isolation).  Census is not a proof."
        ),
    }


def fperp_injective_proved_general() -> bool:
    return False


def k3_spans_Wpp0_proved_general() -> bool:
    return bool(
        fperp_injective_proved_general()
        and theorem_singer_pd_on_F_p_ge_7()["proved"]
        and theorem_p5_k1_fills()["proved"]
    )


def residual_i_closed_via_271() -> bool:
    return False


def hinge_status_271() -> dict:
    return {
        "majority_identity": theorem_majority_identity()["proved"],
        "phase_lock": theorem_phase_lock()["proved"],
        "fibre_dft": theorem_fibre_dft()["proved"],
        "slice_sizes": theorem_slice_sizes()["proved"],
        "phase_frequencies_distinct": theorem_phase_frequencies_distinct()["proved"],
        "singer_pd_on_F_p_ge_7": theorem_singer_pd_on_F_p_ge_7()["proved"],
        "p5_k1_fills": theorem_p5_k1_fills()["proved"],
        "fperp_injective_certified_p_5_7_11": certify_fperp_injective()["certified"],
        "fperp_injective_general": fperp_injective_proved_general(),
        "k3_spans_Wpp0_general": k3_spans_Wpp0_proved_general(),
        "gplus_pd_proved_general": gplus_pd_proved_general(),
        "residual_i_closed": residual_i_closed_via_271(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "F^⊥ leftover: good-μ Gram mixes same-line and mixed-line "
            "coordinates (cross terms certified at p=5).  Bad-μ is "
            "mixed-line only; phase characters are distinct.  Do not "
            "import residual_i."
        ),
    }


def main() -> dict:
    st = hinge_status_271()
    out = {
        "title": (
            "Prop 15.271 k=3 Veronese ker⊆F; Singer PD on F for p≥7; "
            "F^⊥ injectivity OPEN"
        ),
        "proved": st,
        "certified": certify_fperp_injective(),
        "k3_spans_Wpp0_proved_general": k3_spans_Wpp0_proved_general(),
        "residual_i_closed": residual_i_closed_via_271(),
    }
    ev = ROOT / "evidence" / "e1_gmin_m4_prop15271.json"
    ev.parent.mkdir(parents=True, exist_ok=True)
    ev.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.271 k=3 Veronese / F vs F^⊥")
    for k, v in st.items():
        print(f"  {k}: {v}")
    return out


if __name__ == "__main__":
    main()
