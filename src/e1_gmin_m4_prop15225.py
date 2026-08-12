#!/usr/bin/env python3
"""
Prop 15.225 — Constrained projectors P± on E_{±4p}; φ⊥E ⇒ μ_part=μ_mn;
δ-split and pointwise δ-room; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.214, 15.224)
  Master (16p² I − T²) μ = 16 κ on R^{binom(n,4)}, n=p²+1.
  ker(16p²I−T²)=E_{+4p}(T)⊕E_{−4p}(T)=:E_{±4p}.
  Solutions μ=μ_mn+δ with μ_mn=16(16p²I−T²)^+κ ⊥ E_{±4p}, δ∈E_{±4p}.
  μ_part=[(p²−1)κ−2φ]/(p²(p²−5)) solves the master on span{κ,φ} (15.224).
  Tφ=(n+2) star, star=−Tκ/6, T=Tᵀ, κ⊥E_{±4p} (15.214–215, 15.224).

PROVED Max+-free
  A. Constrained projectors on E_{±4p}.
        For v∈E_{±4p} write
          P_{+}v = (4p I + T)v / (8p),   P_{-}v = (4p I − T)v / (8p).
     Then on E_{+4p}: P_{+}=Id, P_{-}=0; on E_{−4p}: P_{+}=0, P_{-}=Id.
     Hence P_{+}+P_{-}=Id and P_{+}P_{-}=0 on E_{±4p}, and
        Tv = 4p (P_{+}−P_{-})v,   T²v = 16p² v.
     (These are the unique orthogonal projectors E_{±4p}→E_{±4p} because T
     acts as ±4p on the two summands.  They are *not* claimed as global
     polynomials that annihilate the rest of σ(T).)  ∎

  B. Unique δ-split.
     Every δ∈E_{±4p} decomposes uniquely as δ=δ_{+}+δ_{-} with
     δ_±=P_±δ∈E_{±4p}, ⟨δ_{+},δ_{-}⟩=0, ‖δ‖₂²=‖δ_{+}‖₂²+‖δ_{-}‖₂², and
     Tδ=4p(δ_{+}−δ_{-}).  ∎

  C. Max± realisation of the split.
     m₄^{+}=m₄^{+,mn}+k_{+} with k_{+}∈E_{+4p},
     m₄^{-}=m₄^{−,mn}+k_{-} with k_{-}∈E_{−4p},
     μ=½(m₄^{+}+m₄^{-})=μ_mn+δ with δ=½(k_{+}+k_{-}),
     so δ_{+}=k_{+}/2 and δ_{-}=k_{-}/2.  ∎

  D. φ ⊥ E_{±4p} (and star ⊥ E_{±4p}).
     Let v∈E_λ with λ∈{±4p}.  T=Tᵀ ⇒ ⟨Tφ,v⟩=λ⟨φ,v⟩.
     Tφ=(n+2)star=−(n+2)/6 · Tκ, so
        ⟨Tφ,v⟩=−(n+2)/6 · ⟨κ,Tv⟩=−(n+2)/6 · λ ⟨κ,v⟩=0
     because κ⊥E_{±4p}.  Thus λ⟨φ,v⟩=0, and λ≠0 ⇒ ⟨φ,v⟩=0.
     Same argument with star=−Tκ/6: ⟨star,v⟩=−λ/6 ⟨κ,v⟩=0.  ∎

  E. μ_part = μ_mn for every prime p≥5 with p²≠5.
     μ_part=aκ+bφ with a=(p²−1)/den, b=−2/den, den=p²(p²−5).
     By (D) and κ⊥E: μ_part ⊥ E_{±4p}.  By 15.224 D, μ_part solves the
     master.  The unique master solution orthogonal to E_{±4p} is μ_mn.
     Hence μ_part=μ_mn, and every Max± solution is
        μ = μ_part + δ,   δ∈E_{±4p}.  ∎
     (Upgrades the p=3,5,7 LSQR census of 15.186/214 to a general theorem.)

  F. Antipodal / φ-free particular coefficient.
        μ_avg_part := (p²−1)/(p²(p²−5)) · κ
     is the φ-average of μ_part (φ→−φ under C→−C with n fixed).  On |κ|=1:
        |μ_avg_part| = (p²−1)/(p²(p²−5)) ≤ 1/(2p)
     for all primes p≥5 (⇔ 2(p²−1)≤p(p²−5) ⇔ 0≤p³−2p²−5p+2; positive at
     p=5 and increasing).  ∎

  G. Pointwise δ-room for the Farkas threshold.
     On |κ|=1 with |φ|≤2(p−2): |μ_part|≤maj:=(p²+4p−9)/(p²(p²−5))≤1/(2p)
     (15.186 H).  Therefore if |δ_S| ≤ room_δ(p) for every |κ_S|=1 four-set,
     with
        room_δ(p) := 1/(2p) − maj
                   = p(p³ − 2p² − 13p + 18) / (2 p³ (p²−5))
                   = (p³ − 2p² − 13p + 18) / (2 p² (p²−5)),
     then |μ_S|=|μ_part_S+δ_S|≤1/(2p).  For all primes p≥5: room_δ(p)>0
     (same cubic as 15.186 H).  ∎

  H. Global L² orthogonality (structure for Aut-SOS / constrained CS).
     ⟨δ,κ⟩=⟨δ,φ⟩=⟨δ,star⟩=⟨δ,μ_part⟩=0.  Combined with (A)–(B): any
     Aut-averaged quadratic form on (δ_{+},δ_{-}) may use these linear
     constraints and Tδ=4p(δ_{+}−δ_{-}) without unconstrained CS.  ∎

CERTIFIED
  I. room_δ>0 and maj≤1/(2p) as Fraction for all primes 5≤p≤300.
  J. Census anchors unchanged: p=5 max|μ|=3/65; p=7 max|μ|=109/2863;
     |μ|≰|μ_part| pointwise at p=7 (so δ can share sign with μ_part).

OPEN (blocks residual i / predicate flip)
  K. Prove for all primes p≥5 one of:
     (i)   |δ_S|≤room_δ(p) on every |κ|=1 four-set (or |μ_part+δ|≤1/(2p)
           by constrained P± / Aut-SOS / Boolean moment constraints);
     (ii)  Path C / k₂/|Max+|≤gap_wick / K₄≤Wick_hi / free-e+ker=sc.
     Then wire residual_i / type_I / gsum via real imports.  Soft-close forbidden.

Until K: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15225.json
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
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)
from e1_gmin_m4_prop15224 import (  # noqa: E402
    theorem_T2_phi_closed,
    theorem_mu_part_solves_master_general,
)


def projector_plus_on_eigen(lambda_T: int, p: int) -> Fraction:
    """P_{+} eigenvalue: (4p + λ)/(8p) on a T-eigenvector with eigenvalue λ."""
    return Fraction(4 * p + lambda_T, 8 * p)


def projector_minus_on_eigen(lambda_T: int, p: int) -> Fraction:
    """P_{-} eigenvalue: (4p − λ)/(8p)."""
    return Fraction(4 * p - lambda_T, 8 * p)


def theorem_constrained_projectors() -> dict:
    """P± restrict to the orthogonal projectors E_{±4p}→E_{±λ}."""
    ok = True
    sample = {}
    for p in [pp for pp in range(5, 80) if is_prime(pp)]:
        # On E_{+4p}
        pp = projector_plus_on_eigen(4 * p, p)
        pm = projector_minus_on_eigen(4 * p, p)
        if pp != 1 or pm != 0:
            ok = False
        # On E_{-4p}
        pp2 = projector_plus_on_eigen(-4 * p, p)
        pm2 = projector_minus_on_eigen(-4 * p, p)
        if pp2 != 0 or pm2 != 1:
            ok = False
        # Sum / product on both
        if pp + pm != 1 or pp2 + pm2 != 1:
            ok = False
        if pp * pm != 0 or pp2 * pm2 != 0:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "P_plus_on_plus": str(pp),
                "P_minus_on_plus": str(pm),
                "P_plus_on_minus": str(pp2),
                "P_minus_on_minus": str(pm2),
            }
    return {
        "proved": ok,
        "theorem": (
            "On E_{±4p}: P_{+}=(4pI+T)/(8p), P_{-}=(4pI−T)/(8p) are the "
            "orthogonal projectors onto E_{+4p} and E_{−4p} respectively; "
            "P_{+}+P_{-}=Id, P_{+}P_{-}=0, T=4p(P_{+}−P_{-})."
        ),
        "by_p_sample": sample,
    }


def theorem_delta_split() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Every δ∈E_{±4p} splits uniquely as δ=δ_{+}+δ_{-} with δ_±=P_±δ, "
            "⟨δ_{+},δ_{-}⟩=0, ‖δ‖₂²=‖δ_{+}‖₂²+‖δ_{-}‖₂², Tδ=4p(δ_{+}−δ_{-})."
        ),
        "depends_on": ["theorem_constrained_projectors"],
    }


def theorem_maxpm_realisation() -> dict:
    return {
        "proved": True,
        "theorem": (
            "m₄^{+}=m₄^{+,mn}+k_{+} (k_{+}∈E_{+4p}), m₄^{-}=m₄^{−,mn}+k_{-} "
            "(k_{-}∈E_{−4p}), μ=μ_mn+δ with δ=½(k_{+}+k_{-}), δ_{+}=k_{+}/2, "
            "δ_{-}=k_{-}/2."
        ),
        "depends_on": ["theorem_delta_split", "prop_15_214_master_spectral"],
    }


def theorem_phi_perp_Epm4p() -> dict:
    """φ ⊥ E_{±4p} and star ⊥ E_{±4p} by T-symmetry + Tφ identity + κ⊥E."""
    # Pure structural theorem; algebraic dependencies recorded.
    t2 = theorem_T2_phi_closed()
    return {
        "proved": t2["proved"],
        "theorem": (
            "For every v∈E_λ, λ∈{±4p}: ⟨Tφ,v⟩=λ⟨φ,v⟩ and "
            "Tφ=(n+2)star=−(n+2)/6·Tκ ⇒ ⟨Tφ,v⟩=−(n+2)λ/6·⟨κ,v⟩=0 "
            "by κ⊥E_{±4p}.  Hence ⟨φ,v⟩=0.  Likewise ⟨star,v⟩=0."
        ),
        "depends_on": [
            "theorem_T_symmetric",
            "theorem_kappa_perp_Epm4p",
            "prop_15_215_Tphi_star",
            "prop_15_215_Tkappa_eq_minus_6_star",
        ],
    }


def theorem_mu_part_equals_mu_mn() -> dict:
    """μ_part = μ_mn for all primes p≥5 (p²≠5)."""
    master = theorem_mu_part_solves_master_general()
    phi_perp = theorem_phi_perp_Epm4p()
    ok = master["proved"] and phi_perp["proved"]
    sample = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        den = p * p * (p * p - 5)
        a = Fraction(p * p - 1, den)
        b = Fraction(-2, den)
        sample[str(p)] = {
            "a": str(a),
            "b": str(b),
            "mu_part": "a*kappa + b*phi",
            "perp_E": "kappa_perp and phi_perp => mu_part_perp",
        }
    return {
        "proved": ok,
        "theorem": (
            "μ_part=aκ+bφ solves the master (15.224) and is orthogonal to "
            "E_{±4p} (κ⊥E and φ⊥E).  Hence μ_part=μ_mn, and every master "
            "solution is μ=μ_part+δ with δ∈E_{±4p}."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_mu_part_solves_master_general",
            "theorem_phi_perp_Epm4p",
            "theorem_kappa_perp_Epm4p",
        ],
    }


def mu_avg_part_coeff(p: int) -> Fraction:
    """Coefficient of κ in the φ-averaged particular solution."""
    return Fraction(p * p - 1, p * p * (p * p - 5))


def theorem_avg_part_beats_thr() -> dict:
    ok = True
    sample = {}
    for p in range(5, 301):
        if not is_prime(p):
            continue
        c = mu_avg_part_coeff(p)
        thr = mu4_threshold(p)
        if c > thr:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101, 199, 293):
            sample[str(p)] = {
                "avg_coeff": str(c),
                "thr": str(thr),
                "beats": c <= thr,
            }
    return {
        "proved": ok,
        "theorem": (
            "On |κ|=1: |μ_avg_part|=(p²−1)/(p²(p²−5))≤1/(2p) for all primes "
            "p≥5.  Proof: 2(p²−1)≤p(p²−5)⇔p³−2p²−5p+2≥0; value 12 at p=5, "
            "derivative 3p²−4p−5>0 on [5,∞)."
        ),
        "by_p_sample": sample,
    }


def room_delta(p: int) -> Fraction:
    """1/(2p) − majorant_mu_part(p) = (p³−2p²−13p+18)/(2p²(p²−5))."""
    return Fraction(p**3 - 2 * p * p - 13 * p + 18, 2 * p * p * (p * p - 5))


def room_delta_via_thr_maj(p: int) -> Fraction:
    return mu4_threshold(p) - majorant_mu_part(p)


def theorem_room_delta_closed() -> dict:
    ok = True
    sample = {}
    maj_ok = theorem_particular_majorant_beats_thr()["proved"]
    phi_ok = theorem_phi_bound_general()["proved"]
    for p in range(5, 301):
        if not is_prime(p):
            continue
        r = room_delta(p)
        r2 = room_delta_via_thr_maj(p)
        if r != r2 or r <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101, 199, 293):
            sample[str(p)] = {
                "room_delta": str(r),
                "maj": str(majorant_mu_part(p)),
                "thr": str(mu4_threshold(p)),
                "positive": r > 0,
            }
    return {
        "proved": ok and maj_ok and phi_ok,
        "theorem": (
            "room_δ(p)=1/(2p)−(p²+4p−9)/(p²(p²−5))"
            "=(p³−2p²−13p+18)/(2p²(p²−5))>0 for all primes p≥5.  "
            "If |δ_S|≤room_δ(p) on every |κ|=1 four-set then "
            "|μ_part+δ|≤1/(2p) on |κ|=1 ⇒ residual-(i) Farkas (15.176)."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_186_majorant",
            "prop_15_186_phi_bound",
            "theorem_mu_part_equals_mu_mn",
        ],
        "sufficient_for_residual_i": False,
        "note": (
            "room_δ formula and positivity proved; the inequality "
            "|δ_S|≤room_δ(p) on |κ|=1 is OPEN."
        ),
    }


def theorem_L2_orthogonality() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For every master residual δ∈E_{±4p}: "
            "⟨δ,κ⟩=⟨δ,φ⟩=⟨δ,star⟩=⟨δ,μ_part⟩=0.  "
            "Proof: κ,φ,star ⊥ E (15.214 C + theorem_phi_perp_Epm4p); "
            "μ_part∈span{κ,φ}."
        ),
        "depends_on": [
            "theorem_phi_perp_Epm4p",
            "theorem_kappa_perp_Epm4p",
            "theorem_mu_part_equals_mu_mn",
        ],
        "note": (
            "Linear constraints for constrained CS / Aut-averaged SOS on "
            "(δ_{+},δ_{-}); does not by itself give the L^∞ room bound."
        ),
    }


def residual_i_closed_via_225() -> bool:
    """False until |δ|≤room_δ or another hinge is proved."""
    return False


def hinge_status_225() -> dict:
    return {
        "constrained_projectors": theorem_constrained_projectors()["proved"],
        "delta_split": theorem_delta_split()["proved"],
        "phi_perp_Epm4p": theorem_phi_perp_Epm4p()["proved"],
        "mu_part_equals_mu_mn": theorem_mu_part_equals_mu_mn()["proved"],
        "avg_part_beats_thr": theorem_avg_part_beats_thr()["proved"],
        "room_delta_closed": theorem_room_delta_closed()["proved"],
        "L2_orthogonality": theorem_L2_orthogonality()["proved"],
        "residual_i_closed_via_225": residual_i_closed_via_225(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|δ_S|≤room_δ(p) on |κ|=1 for all p≥5 (constrained P± / Aut-SOS); "
            "or Path C / k₂ gap / K₄ / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_constrained_projectors()
    b = theorem_delta_split()
    c = theorem_maxpm_realisation()
    d = theorem_phi_perp_Epm4p()
    e = theorem_mu_part_equals_mu_mn()
    f = theorem_avg_part_beats_thr()
    g = theorem_room_delta_closed()
    h = theorem_L2_orthogonality()
    status = hinge_status_225()
    out = {
        "prop": "15.225",
        "title": (
            "Constrained P± on E_{±4p}; φ⊥E ⇒ μ_part=μ_mn; "
            "δ-room formula; residual-(i) OPEN"
        ),
        "proved": {
            "constrained_projectors": a["proved"],
            "delta_split": b["proved"],
            "maxpm_realisation": c["proved"],
            "phi_perp_Epm4p": d["proved"],
            "mu_part_equals_mu_mn": e["proved"],
            "avg_part_beats_thr": f["proved"],
            "room_delta_closed": g["proved"],
            "L2_orthogonality": h["proved"],
            "residual_i_closed": residual_i_closed_via_225(),
        },
        "algebra": {
            "P_plus": "(4p I + T)/(8p) on E_pm4p",
            "P_minus": "(4p I - T)/(8p) on E_pm4p",
            "mu": "mu_part + delta, delta in E_pm4p",
            "room_delta": "(p^3 - 2p^2 - 13p + 18)/(2 p^2 (p^2 - 5))",
            "room_delta_p5": str(room_delta(5)),
            "room_delta_p7": str(room_delta(7)),
            "avg_coeff_p5": str(mu_avg_part_coeff(5)),
            "avg_coeff_p7": str(mu_avg_part_coeff(7)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_225(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "P±, φ⊥E, μ_part=μ_mn, and room_δ formula proved Max+-free.  "
            "Pointwise |δ|≤room_δ (or other residual-i hinge) still OPEN.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15225.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.225  residual_i={residual_i_closed_via_225()}")
    print(
        f"  P±={a['proved']} phi_perp={d['proved']} "
        f"mu_part=mu_mn={e['proved']} room_delta={g['proved']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
