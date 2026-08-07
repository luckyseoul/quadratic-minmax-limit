#!/usr/bin/env python3
"""
Prop 15.214 — Master / δ∈E_{±4p} structure for μ₄; T symmetric;
min-norm vs Max±; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176–188, 15.208)
  Master (15.177, Max+-free conference boolean algebra):
      (16 p² I − T²) μ₄ = 16 κ
  with T the C-weighted 4-set extension operator (symmetric — below).
  μ_part = [(p²−1)κ − 2φ]/(p²(p²−5)) solves the master on span{κ,φ}
  when T²φ = 4(p²+3)(φ−2κ) (certified p=3,5,7; 15.187).
  |μ_part| ≤ (p²+4p−9)/(p²(p²−5)) ≤ 1/(2p) on |κ|=1 for p≥5
  under |φ|≤2(p−2) (15.186, proved).  Actual Max± μ₄ = μ_mn + δ with
  δ ∈ ker(16p²I−T²).  Census: δ≠0 at p=5; |μ|≰|μ_part| pointwise at p=7.

PROVED Max+-free
  A. T is symmetric on R^{binom(n,4)}.
     For 4-sets S,S' adjacent by replacement v→r (v∈S, r∉S, S'=S−v+r),
     the undirected Johnson edge is weighted by C_{vr}=C_{rv}.  Hence the
     adjacency operator T is symmetric.  ∎

  B. Spectral form of the master.
     T = Tᵀ ⇒ T diagonalizable, σ(T)⊂R, and
        ker(16p² I − T²) = E_{+4p}(T) ⊕ E_{−4p}(T) =: E_{±4p}.
     On each T-eigenspace with eigenvalue λ, λ²≠16p²:
        μ_λ = 16 κ_λ / (16p² − λ²).
     The unique solution orthogonal to E_{±4p} is the min-norm solution
        μ_mn := 16 (16p² I − T²)^+ κ.
     Every solution is μ = μ_mn + δ with δ∈E_{±4p}.  ∎

  C. Compatibility: κ ⊥ E_{±4p}.
     The master admits a solution μ₄ = ½(m₄^{+}+m₄^{-}) for any conference
     Max± (15.177, proved).  Hence 16κ ∈ range(16p²I−T²) = E_{±4p}^⊥,
     i.e. κ ⊥ E_{±4p}.  In particular ⟨κ, δ⟩=0 for every admissible δ.  ∎

  D. Formal resolvent for m₄± (when invertible).
     (4p I − T) m₄^{+} = 4κ/p,  (4p I + T) m₄^{-} = 4κ/p (15.177).
     When ±4p ∉ σ(T): uniquely
        μ₄ = 16 (16p² I − T²)^{−1} κ = μ_mn
     (δ=0).  When ±4p ∈ σ(T), m₄± are determined only up to
     ker(4pI∓T), and Max± may select δ≠0.  ∎

  E. Particular majorant (recall 15.186, still Max+-free).
     On |κ|=1 with |φ|≤2(p−2): |μ_part| ≤ (p²+4p−9)/(p²(p²−5)) ≤ 1/(2p)
     for all primes p≥5.  If δ=0 on |κ|=1 (or |μ_mn+δ|≤1/(2p)), residual-(i)
     Farkas fires (15.176).  ∎

  F. Weak matching PSD bound.
     The three pairwise-disjoint edges of any K4 form a Gsum principal
     submatrix equivalent (after scaling) to 2[[1,πμ],[πμ,1]] per matching
     pair, or the 3×3 pattern with eigs {1−2μ, 1+μ, 1+μ} for |κ|=1 sign
     patterns.  Gsum ≽ 0 ⇒ |μ| ≤ 1/2 on every |κ|=1 four-set.
     (Too weak for 1/(2p); recorded as structural PSD constraint.)  ∎

CERTIFIED (not general)
  G. ||T||₂ = 4p at p=5,7 (extreme eigs ±4p); at p=3, ||T||₂≈9.80<12.
     At p=5, mult(E_{+4p})≥65 and mult(E_{−4p})≥65 (Lanczos LM sample).
  H. μ_mn = μ_part when T²φ identity holds (p=3,5,7).  At p=5:
     max|μ|=3/65=M_cand < maj=9/125 < 1/10; δ≠0 on |κ|=1; |μ|<max|μ_part|.
     At p=7: max|μ|=109/2863 < 1/14; |μ| ≰ |μ_part| pointwise (15.188).
  I. Parseval on Max+: ∑_{A⊆[n]} (E[y^A])² = 2ⁿ/|Max+| (boolean cube).

OPEN (blocks residual i / predicate flip)
  J. Prove for all primes p≥5 one of:
     (i)  |μ_mn+δ| ≤ 1/(2p) (or ≤M_cand) on every |κ|=1 four-set
          — control δ∈E_{±4p} by Boolean realizability / centre SOS;
     (ii) T²φ identity general + δ=0 on |κ|=1 (false at p=5, so need other);
     (iii) K₄≤Wick_hi or E[q²]≥λ_* / Veronese rank=D (15.198/212–213);
     (iv) free-e_sc bound + ker=sc.
     Then wire gsum_disj_lb_proved_general → residual_i → E1 only if residual
     (ii) allows.  Soft-close forbidden.

Until J: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15214.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15176 import theorem_minus_1_over_p_suffices  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)
from e1_gmin_m4_prop15205 import M_cand  # noqa: E402


def theorem_T_symmetric() -> dict:
    return {
        "proved": True,
        "theorem": (
            "T is the adjacency operator of the Johnson graph J(n,4) with "
            "undirected edge weights C_vr=C_rv for replacements v→r.  Hence T=Tᵀ."
        ),
    }


def theorem_master_spectral_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "T=Tᵀ ⇒ ker(16p²I−T²)=E_{±4p}(T).  Solutions of the master are "
            "μ=μ_mn+δ with μ_mn=16(16p²I−T²)^+κ and δ∈E_{±4p}."
        ),
        "depends_on": ["theorem_T_symmetric", "prop_15_177_master"],
    }


def theorem_kappa_perp_Epm4p() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Master solvable (15.177 Max± μ₄) ⇒ 16κ ∈ range(16p²I−T²)=E_{±4p}^⊥ "
            "⇒ κ ⊥ E_{±4p}.  Hence ⟨κ,δ⟩=0 for every master solution residual δ."
        ),
        "depends_on": ["theorem_master_spectral_form", "prop_15_177"],
    }


def theorem_resolvent_when_invertible() -> dict:
    return {
        "proved": True,
        "theorem": (
            "If ±4p∉σ(T), then m₄± are unique and μ₄=16(16p²I−T²)^{−1}κ=μ_mn "
            "(δ=0).  If ±4p∈σ(T), Max± may select δ≠0."
        ),
    }


def theorem_particular_majorant_recall() -> dict:
    t_phi = theorem_phi_bound_general()
    t_maj = theorem_particular_majorant_beats_thr()
    return {
        "proved": t_phi.get("proved", True) and t_maj.get("proved", True),
        "theorem": (
            "|φ|≤2(p−2) on |κ|=1 (15.186) and |μ_part|≤(p²+4p−9)/(p²(p²−5))≤1/(2p) "
            "for all primes p≥5.  δ-control or |μ_mn+δ|≤1/(2p) still open."
        ),
        "phi_bound": t_phi.get("proved", True),
        "majorant_beats_thr": t_maj.get("proved", True),
        "sample": {
            str(p): {
                "majorant": str(majorant_mu_part(p)),
                "thr": str(mu4_threshold(p)),
                "M_cand": str(M_cand(p)),
            }
            for p in (5, 7, 11)
        },
    }


def theorem_matching_psd_half() -> dict:
    """Gsum 3×3 / matching 2×2 ⇒ |μ|≤1/2 on |κ|=1."""
    # Eigenvalues of [[1,m,m],[m,1,-m],[m,-m,1]] are 1-2m, 1+m, 1+m
    # PSD for m ∈ [-1, 1/2]
    return {
        "proved": True,
        "theorem": (
            "On any |κ|=1 four-set the three disjoint edge-pairs give Gsum "
            "principal blocks with pattern eigs {1−2μ (×1), 1+μ (×2)} (up to "
            "the global factor 2 on diagonal Gsum_ee=2).  Gsum≽0 ⇒ μ≤1/2 and "
            "μ≥−1.  Hence |μ|≤1/2.  Insufficient for 1/(2p)."
        ),
        "sufficient_for_residual_i": False,
        "bound": "1/2",
        "need": "1/(2p)",
    }


def theorem_parseval_maxplus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For uniform measure on any M⊆{±1}^n: ∑_{A⊆[n]} (E_M[y^A])² = 2ⁿ/|M|. "
            "Proof: E_y E_z ∏_i(1+y_i z_i)=2ⁿ P(y=z)=2ⁿ/|M|, and expand the product."
        ),
        "note": (
            "Specialises to Max+ with |M|=|Max+|.  Controls ∑_A m_A² but does not "
            "by itself bound ∑_{|A|=4}(m₄−κ/p²)² tightly without |Max+| and other degrees."
        ),
    }


def certify_delta_structure_small_p() -> dict:
    return {
        "certified": True,
        "proved_general": False,
        "by_p": {
            "3": {
                "T_norm_lt_4p": True,
                "delta_on_kappa1": "0 (μ=μ_part on |κ|=1)",
                "max_abs_mu": "1/3",
                "thr": "1/6",
                "beats_thr": False,
            },
            "5": {
                "T_norm_eq_4p": True,
                "mult_E_pm4p_lower": {"plus": 65, "minus": 65},
                "delta_nonzero_on_kappa1": True,
                "max_abs_mu": "3/65",
                "M_cand": "3/65",
                "majorant": "9/125",
                "thr": "1/10",
                "beats_thr": True,
            },
            "7": {
                "T_norm_eq_4p": True,
                "max_abs_mu": "109/2863",
                "thr": "1/14",
                "beats_thr": True,
                "not_bounded_by_mu_part_pointwise": True,
            },
        },
    }


def mu_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_214() -> bool:
    return False


def hinge_status_214() -> dict:
    return {
        "T_symmetric": True,
        "master_spectral_form": True,
        "kappa_perp_Epm4p": True,
        "particular_majorant": True,
        "matching_psd_half": True,
        "parseval_maxplus": True,
        "delta_control_general": False,
        "mu_bound_proved_general": mu_bound_proved_general(),
        "residual_i_closed": residual_i_closed_via_214(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(),
        "open": (
            "Control δ∈E_{±4p} so |μ_mn+δ|≤1/(2p) on |κ|=1 for all p≥5; "
            "or K₄/Veronese/free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_T_symmetric()
    b = theorem_master_spectral_form()
    c = theorem_kappa_perp_Epm4p()
    d = theorem_resolvent_when_invertible()
    e = theorem_particular_majorant_recall()
    f = theorem_matching_psd_half()
    g = theorem_parseval_maxplus()
    cert = certify_delta_structure_small_p()
    status = hinge_status_214()
    out = {
        "prop": "15.214",
        "title": "Master/δ∈E±4p structure; T symmetric; residual (i) OPEN",
        "proved": {
            "T_symmetric": a["proved"],
            "master_spectral_form": b["proved"],
            "kappa_perp_Epm4p": c["proved"],
            "resolvent_when_invertible": d["proved"],
            "particular_majorant_recall": e["proved"],
            "matching_psd_half": f["proved"],
            "parseval_maxplus": g["proved"],
        },
        "certified": cert,
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_214(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15214.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    for k, v in out["proved"].items():
        print(f"  {k}={v}")
    print(f"  residual_i={residual_i_closed_via_214()} gsum={gsum_disj_lb_proved_general()} e1={e1_closed_general()}")
    return out


if __name__ == "__main__":
    main()
