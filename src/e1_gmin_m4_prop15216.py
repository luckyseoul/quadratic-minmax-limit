#!/usr/bin/env python3
"""
Prop 15.216 — Residual (i) K₄≤16n²−14n path: reductions + FATAL GAP on R≤2p.

Max+-free reductions (any conference of order n=p²+1 on Paley for Tκ spectrum):
  (15.215) Tκ on λ²=4(p²+15); η_part formula; η=η_part+δ, δ∈E_{4p}.
  (this)  κ∈V_λ⊕ker(T); ⟨κ,Tκ⟩=0; R_ke=128p/(3p²+17)≤2p (p≥5);
          crude ‖κ‖₂²/p⁴ ≤ thr-η-budget; Q≤10 ⇒ dual-eq empty (mass-min).

FATAL GAP (openai-referee cold AI-test 2026-08-07, 2× BLOCK):
  m₄=ke+δ with δ∈E_{4p}, ⟨ke,δ⟩=0 gives
    R = (R_ke‖ke‖² + 4p‖δ‖²)/(‖ke‖²+‖δ‖²)
  a convex combination of R_ke and **4p**.  R_ke≤2p implies only R≤4p,
  **not** R≤2p, unless ‖δ‖ is controlled so that
    (2p−R_ke)‖ke‖² ≥ 2p‖δ‖².
  The old claim "R_ke≤2p ⇒ R≤2p" was false algebra.  Soft-close forbidden.

OPEN hinge for this path: prove R(m₄)≤2p (control δ), or bound ‖η‖/K₄
without that Rayleigh inequality, or use another residual-(i) hinge.

Does **not** flip residual_i / type_I / gsum / E1 / L.
Writes evidence/e1_gmin_m4_prop15216.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15194 import dual_eq_need_kappa_e  # noqa: E402
from e1_gmin_m4_prop15197 import K4_threshold_for_Q10, Q_from_K4  # noqa: E402
from e1_gmin_m4_prop15198 import C_diag, wick_hi  # noqa: E402
from e1_gmin_m4_prop15215 import (  # noqa: E402
    eta_part_norm_sq,
    theorem_T3_kappa,
    theorem_Tphi_eq_nplus2_star,
    theorem_eta_part_formula,
    theorem_eta_part_lt_wick_budget,
    theorem_n3_count,
    theorem_star_values_by_exhaustion,
    theorem_sum_kappa_sq,
)


# ---------------------------------------------------------------------------
# A. sum_S star(S) κ(S) = 0 (⇒ ⟨κ, Tκ⟩=0)
# ---------------------------------------------------------------------------


def theorem_star_kappa_orthogonal() -> dict:
    """
    ∑_S star(S) κ(S) = 0 for any symmetric conference.
    Proof: expand star_v·π_m; each term is a product of four C-edges;
    for star_a·π_1 = C_ac C_ad C_cd the sum over distinct a,b,c,d factors as
    (n−3)·∑_{a,c,d distinct} C_ac C_ad C_cd and for fixed a≠c,
    ∑_{d≠a,c} C_ad C_cd = (C²)_{ac} − 0 = 0.  All vertex/matching pairs vanish.
    """
    # Exhaustion check on K4 labels is automatic; global sum uses C².
    # Structural smoke: 64 labels of star*κ average — not the full proof.
    edges = list(combinations(range(4), 2))
    total = 0
    for bits in product([-1, 1], repeat=6):
        lab = dict(zip(edges, bits))

        def E(i, j):
            return lab[tuple(sorted((i, j)))]

        matchings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
        kap = sum(E(a, b) * E(c, d) for (a, b), (c, d) in matchings)
        st = 0
        for v in range(4):
            pr = 1
            for u in range(4):
                if u != v:
                    pr *= E(v, u)
            st += pr
        total += st * kap
    # On abstract K4, sum over labelings of star*κ need not be 0 (finite);
    # the conference sum uses C² cancellation across the ambient graph.
    return {
        "proved": True,
        "theorem": (
            "∑_S star(S)κ(S)=0 for any symmetric conference C of order n.  "
            "Expand star_v·π_m; each monomial sums to 0 by C² (for a≠c, "
            "∑_{d≠a,c} C_ad C_cd=(C²)_{ac}=0).  Hence ⟨κ,Tκ⟩=−6∑star κ=0."
        ),
        "depends_on": ["C_squared_eq_n_minus_1_I", "Tkappa_eq_minus_6_star"],
        "k4_label_sum_star_kappa_not_used": total,  # diagnostic only
    }


def theorem_kappa_spectral_support() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Tκ∈V_λ:=E_λ⊕E_{−λ} with λ²=4(p²+15) (15.215 G).  T=Tᵀ is "
            "diagonalizable.  Write κ=P_V κ+P_other κ.  Then Tκ=T P_V κ+T P_other κ∈V_λ "
            "forces T P_other κ∈V_λ.  If P_other involves eigenvalues ∉{0,±λ}, the "
            "image under T would leave V_λ unless P_other κ=0.  Absorbing ker(T) into "
            "eigenvalue 0: κ∈V_λ⊕ker(T).  Consequently κ_±=(Tκ)_±/λ on ±λ."
        ),
        "depends_on": ["theorem_T3_kappa", "T_symmetric"],
    }


# ---------------------------------------------------------------------------
# B. R_ke ≤ 2p and R ≤ 2p
# ---------------------------------------------------------------------------


def R_ke(p: int) -> Fraction:
    """128 p / (3 p² + 17)."""
    return Fraction(128 * p, 3 * p * p + 17)


def theorem_R_ke_formula() -> dict:
    """
    On ke=κ/p²+η_part (spectral calculus on V_λ with ⟨κ,Tκ⟩=0 and κ_±=Tκ_±/λ):
      ⟨η,Tκ⟩ = ‖Tκ‖₂² · p / (3 p² (p²−5)),
      ⟨η,Tη⟩ = ‖Tκ‖₂² / p⁴ · 2p(p²+15)/(9(p²−5)²),
      ⟨κ,η⟩ = ‖Tκ‖₂² / (12 p² (p²−5)),
      ‖ke‖₂² and ⟨ke,T ke⟩ yield R_ke=128p/(3p²+17).
    """
    ok = True
    sample = {}
    # Cross-check against algebraic derivation for several p
    for p in [5, 7, 11, 13]:
        n = n_of(p)
        Tk2 = 6 * n * (n - 1) * (n - 2) * (n - 6)
        k2 = n * (n - 1) * (n - 2) * (n - 5) // 8
        ep2 = eta_part_norm_sq(p)
        cross = Fraction(Tk2 * p, 3 * p * p * (p * p - 5))
        Teta = Fraction(Tk2, p**4) * Fraction(2 * p * (p * p + 15), 9 * (p * p - 5) ** 2)
        Tke = 2 * cross / Fraction(p * p) + Teta
        kap_eta = Fraction(Tk2, p * p * 12 * (p * p - 5))
        ke2 = Fraction(k2, p**4) + ep2 + 2 * kap_eta / Fraction(p * p)
        R = Tke / ke2
        if R != R_ke(p):
            ok = False
        sample[str(p)] = {"R_ke": str(R), "closed": str(R_ke(p)), "match": R == R_ke(p)}
    return {
        "proved": ok,
        "theorem": (
            "R_ke := Rayleigh of T on κ/p²+η_part equals 128p/(3p²+17) for primes p≥5."
        ),
        "formula": "128p/(3p^2+17)",
        "by_p_sample": sample,
        "depends_on": [
            "theorem_star_kappa_orthogonal",
            "theorem_kappa_spectral_support",
            "theorem_eta_part_formula",
        ],
    }


def theorem_R_ke_le_2p() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        # 2p - R = 2p(3p²−47)/(3p²+17); 3p²−47≥28 at p=5
        holds = R_ke(p) <= 2 * p and (3 * p * p - 47) > 0
        if not holds:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "R_ke": str(R_ke(p)),
                "two_p": 2 * p,
                "le": R_ke(p) <= 2 * p,
            }
    return {
        "proved": ok,
        "theorem": (
            "R_ke=128p/(3p²+17)≤2p for all primes p≥5, since "
            "2p−R_ke=2p(3p²−47)/(3p²+17) and 3p²−47≥28>0."
        ),
        "by_p_sample": sample,
    }


def theorem_R_m4_le_2p() -> dict:
    """
    NOT proved.  Convex combination of R_ke≤2p and 4p is ≤4p, not ≤2p.

    Required for R≤2p: (2p−R_ke)‖ke‖₂² ≥ 2p‖δ‖₂².
    (Old writeup claimed (R_ke−2p)A ≤ −2p B is automatic from R_ke≤2p;
    that confuses A(R_ke−2p)≤0 with A(R_ke−2p)≤−2p B.)
    """
    return {
        "proved": False,
        "fatal_gap": (
            "R=(R_ke‖ke‖²+4p‖δ‖²)/(‖ke‖²+‖δ‖²) is a convex combination of "
            "R_ke and 4p.  R_ke≤2p ⇏ R≤2p when δ∈E_{4p} may be nonzero.  "
            "Need (2p−R_ke)‖ke‖²≥2p‖δ‖² or another η/K₄ bound."
        ),
        "theorem": (
            "OPEN: R(m₄)≤2p for Max+ m₄=ke+δ.  R is a convex combination of "
            "R_ke≤2p and 4p; this yields R≤4p only.  Sufficient condition "
            "(2p−R_ke)‖ke‖₂²≥2p‖δ‖₂² is unproved Max+-free."
        ),
        "depends_on": ["theorem_R_ke_le_2p"],
        "ai_test": "openai-referee deep_review + math_review BLOCK 2026-08-07",
    }


# ---------------------------------------------------------------------------
# C. ‖η‖₂² ≤ ‖κ‖₂²/p⁴ ≤ thr budget ⇒ K₄≤thr
# ---------------------------------------------------------------------------


def kappa_norm_sq(p: int) -> int:
    n = n_of(p)
    return n * (n - 1) * (n - 2) * (n - 5) // 8


def crude_eta_ub(p: int) -> Fraction:
    """‖κ‖₂²/p⁴."""
    return Fraction(kappa_norm_sq(p), p**4)


def thr_eta_budget(p: int) -> Fraction:
    """∑η² budget for K₄≤16n²−14n."""
    n = n_of(p)
    thr = K4_threshold_for_Q10(p)
    base = 12 * n * n - 48 * n + C_diag(p)
    return (thr - base) / 24


def theorem_eta_le_crude_when_R_le_2p() -> dict:
    return {
        "proved": True,  # conditional implication is correct
        "conditional": True,
        "theorem": (
            "From master (4pI−T)m₄=4κ/p: ‖η‖₂²=‖m₄‖₂²(R/(2p)−1)+‖κ‖₂²/p⁴ "
            "with η=m₄−κ/p² and R=⟨m₄,Tm₄⟩/‖m₄‖₂².  If R≤2p then "
            "R/(2p)−1≤0 so ‖η‖₂²≤‖κ‖₂²/p⁴.  (Antecedent R≤2p is OPEN.)"
        ),
        "depends_on": ["prop_15_177_master", "theorem_R_m4_le_2p"],
    }


def theorem_crude_le_thr_eta_budget() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        # Algebra: thr_bud − crude = (p²+1)(p²+9)/24 > 0
        diff = thr_eta_budget(p) - crude_eta_ub(p)
        expected = Fraction((p * p + 1) * (p * p + 9), 24)
        if diff != expected or diff <= 0:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "crude": str(crude_eta_ub(p)),
                "thr_eta_bud": str(thr_eta_budget(p)),
                "diff": str(diff),
                "formula_diff": str(expected),
            }
    return {
        "proved": ok,
        "theorem": (
            "‖κ‖₂²/p⁴ ≤ thr-η-budget for all primes p≥5, with difference "
            "(p²+1)(p²+9)/24>0.  (Explicit Fraction algebra from "
            "∑κ²=n(n−1)(n−2)(n−5)/8 and K₄ thr partition.)"
        ),
        "by_p_sample": sample,
    }


def theorem_K4_le_thr() -> dict:
    return {
        "proved": False,
        "theorem": (
            "OPEN for all primes p≥5: ∑_S η_S² ≤ thr-η-budget ⇒ K₄≤16n²−14n.  "
            "Crude ‖κ‖₂²/p⁴≤thr-η-budget is proved, but ‖η‖₂²≤crude needs R≤2p "
            "(OPEN: fatal Rayleigh gap on δ∈E_{4p})."
        ),
        "depends_on": [
            "theorem_eta_le_crude_when_R_le_2p",
            "theorem_crude_le_thr_eta_budget",
            "theorem_R_m4_le_2p",
            "prop_15_198_partition",
        ],
    }


def theorem_Q_le_10() -> dict:
    return {
        "proved": False,
        "theorem": (
            "K₄≤16n²−14n ⇔ Q_e≤10 (15.197 C).  K₄≤thr OPEN ⇒ Q_e≤10 OPEN."
        ),
        "depends_on": ["theorem_K4_le_thr", "prop_15_197"],
    }


# ---------------------------------------------------------------------------
# D. Q≤10 ⇒ dual-eq empty for all p≥5
# ---------------------------------------------------------------------------


def theorem_Q10_blocks_dual_eq_all_p() -> dict:
    """
    Under a≥−2, ∑_{disj}a=n−2, ∑a²≤10: only two-level k∈{0,1} is admissible
    for p≥5, and mass-min at k=1 is −2−α(n−2) > −2(2−α)=−4+2α.
    """
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        S = n - 2
        n_disj = (n - 2) * (n - 3) // 2
        # k=2: q-10 = 2(9n+2)/(n²-5n+2) > 0 for n≥26
        r2 = Fraction(S + 4, n_disj - 2)
        q2 = 4 * 2 + (n_disj - 2) * r2 * r2
        # k=1: q < 10
        r1 = Fraction(S + 2, n_disj - 1)
        q1 = 4 + (n_disj - 1) * r1 * r1
        alpha = alpha_particular(p)
        M = dual_eq_need_kappa_e(p)
        mass_k1 = -2 - alpha * S
        target = -2 * M
        margin = mass_k1 - target  # = 2 - alpha*n = 2 - 6/p
        holds = (
            q2 > 10
            and q1 <= 10
            and mass_k1 > target
            and margin == Fraction(2) - Fraction(6, p)
        )
        if not holds:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "q_k2": str(q2),
                "q_k1": str(q1),
                "mass_k1": str(mass_k1),
                "target": str(target),
                "margin": str(margin),
                "blocks": holds,
            }
    # Algebra of q2-10 and margin
    # q2-10 = 2(9n+2)/(n²-5n+2) proved by simplify
    n = 26
    q2_m10 = Fraction(2 * (9 * n + 2), n * n - 5 * n + 2)
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5 (n=p²+1≥26): under a≥−2, ∑a=n−2, ∑a²≤10, "
            "the only admissible two-level k with r≥−2 are k=0,1, because "
            "q(k=2)−10=2(9n+2)/(n²−5n+2)>0.  The worst mass-min is at k=1: "
            "−2−α(n−2), and −2−α(n−2)−(−2(2−α))=2−6/p>0.  Hence the "
            "box+mass minimum of ∑aκ is >−2(2−α), so dual-eq is empty (15.195)."
        ),
        "by_p_sample": sample,
        "q2_minus_10_at_n26": str(q2_m10),
        "depends_on": ["prop_15_195_mass_criterion", "prop_15_196_twolevel"],
    }


def theorem_dual_eq_empty_all_p() -> dict:
    return {
        "proved": False,
        "theorem": (
            "OPEN: dual-eq empty for all primes p≥5 via Q_e≤10.  The implication "
            "Q≤10 ⇒ dual-eq empty is proved (theorem_Q10_blocks_dual_eq_all_p); "
            "Q_e≤10 itself is OPEN (depends on K₄≤thr / R≤2p gap)."
        ),
        "depends_on": ["theorem_Q_le_10", "theorem_Q10_blocks_dual_eq_all_p"],
    }


# ---------------------------------------------------------------------------
# E. Residual (i) status + wiring (honest: OPEN)
# ---------------------------------------------------------------------------


def residual_i_dual_eq_empty_proved_general() -> bool:
    """False until R(m₄)≤2p or alternate K₄/dual-eq hinge is Max+-free proved."""
    return False


def residual_i_closed_via_216() -> bool:
    return residual_i_dual_eq_empty_proved_general()


def theorem_residual_i_closed() -> dict:
    # Re-check load-bearing algebra; full chain fails on R_m4_le_2p
    a = theorem_star_values_by_exhaustion()
    b = theorem_sum_kappa_sq()
    c = theorem_Tphi_eq_nplus2_star()
    d = theorem_T3_kappa()
    e = theorem_star_kappa_orthogonal()
    f = theorem_kappa_spectral_support()
    g = theorem_R_ke_formula()
    h = theorem_R_ke_le_2p()
    i = theorem_R_m4_le_2p()
    j = theorem_eta_le_crude_when_R_le_2p()
    k = theorem_crude_le_thr_eta_budget()
    ell = theorem_K4_le_thr()
    m = theorem_Q_le_10()
    n_ = theorem_Q10_blocks_dual_eq_all_p()
    o = theorem_dual_eq_empty_all_p()
    reductions_ok = all(
        t.get("proved", False)
        for t in (a, b, c, d, e, f, g, h, k, n_)
    )
    # j is conditional-true; i,ell,m,o are the OPEN hinge cluster
    all_ok = reductions_ok and i["proved"] and ell["proved"] and m["proved"] and o["proved"]
    return {
        "proved": False,
        "theorem": (
            "Residual (i) via K₄≤16n²−14n is OPEN.  Reductions OK: R_ke≤2p, "
            "crude≤thr-η-budget, Q≤10⇒dual-eq empty.  FATAL GAP: R(m₄)≤2p false "
            "as stated (convex combo with 4p).  Predicates residual_i/type_I stay "
            "False.  gsum_disj_lb False.  E1/L OPEN."
        ),
        "chain_ok": all_ok,
        "reductions_ok": reductions_ok,
        "fatal_gap": i.get("fatal_gap"),
        "depends_on": ["prop_15_215", "prop_15_195", "prop_15_196", "prop_15_197", "prop_15_198"],
    }


def hinge_status_216() -> dict:
    return {
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "residual_i_closed_via_216": residual_i_closed_via_216(),
        "K4_le_thr": False,
        "Q_le_10": False,
        "R_m4_le_2p": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "residual_ii_full": False,
        "L_open": True,
        "note": (
            "Residual (i) OPEN — fatal Rayleigh gap on δ∈E_{4p} (AI-test BLOCK).  "
            "E1/L OPEN.  gsum_disj_lb not flipped."
        ),
    }


def main() -> dict:
    t = theorem_residual_i_closed()
    status = hinge_status_216()
    out = {
        "prop": "15.216",
        "title": "Residual (i) K₄ thr path OPEN — fatal Rayleigh gap on δ∈E_{4p}",
        "proved": {
            "star_kappa_orthogonal": theorem_star_kappa_orthogonal()["proved"],
            "kappa_spectral_support": theorem_kappa_spectral_support()["proved"],
            "R_ke_formula": theorem_R_ke_formula()["proved"],
            "R_ke_le_2p": theorem_R_ke_le_2p()["proved"],
            "R_m4_le_2p": theorem_R_m4_le_2p()["proved"],
            "eta_le_crude_conditional": theorem_eta_le_crude_when_R_le_2p()["proved"],
            "crude_le_thr_budget": theorem_crude_le_thr_eta_budget()["proved"],
            "K4_le_thr": theorem_K4_le_thr()["proved"],
            "Q_le_10": theorem_Q_le_10()["proved"],
            "Q10_blocks_dual_eq": theorem_Q10_blocks_dual_eq_all_p()["proved"],
            "dual_eq_empty": theorem_dual_eq_empty_all_p()["proved"],
            "residual_i_closed": t["proved"],
        },
        "fatal_gap": theorem_R_m4_le_2p().get("fatal_gap"),
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_216(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close_E1_L": False,
        "soft_close_residual_i": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15216.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    for k, v in out["proved"].items():
        print(f"  {k}={v}")
    print(f"  residual_i={residual_i_closed_via_216()}")
    print(f"  gsum={gsum_disj_lb_proved_general()} e1={e1_closed_general()}")
    return out


if __name__ == "__main__":
    main()
