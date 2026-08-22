#!/usr/bin/env python3
"""
Leftover 1 hinge — Max+-free reductions for
  (i)  GLOBAL QVAR on the full Max+ ensemble (all k mixed; λ=0 unsplit)
  (ii) principal ||δ||² room after QVAR
Do **not** number this unless both estimates are proved for all primes p≥5.
Per-stratum QVAR on every k≥7 is false at (41,7) and is **not** equivalent
to (i).

Does **not** treat G_{u,disj} as a Gram.  Aut-Schur / Gsum / pairing stay
their own units.  L follows the four-leftover AND and may close if leftover
1, leftover 2/3, and lemma D are all actually proved.  p=13 orbits are not
a close.

IDENTITIES (proved Max+-free Fraction)
  A. λ_exc = 32 E|Z_ψ|² / [q(q-1)]  (15.589 E).  QVAR
        E|Z_ψ|² ≥ 3q(q-1)/16  ⇔  λ_exc ≥ 6.
     Fail: drop the 16; replace 32 by 16.
  B. ||A_ψ||_HS² = q(q-1)/32, tr A_ψ=0 (15.589 H).  Spherical average
        V_sph = q(q-1)(q+1)/(4(q+5))
     strictly exceeds the QVAR threshold for every prime p≥5:
        V_sph − 3q(q-1)/16 = q(q-1)(q-11)/(16(q+5)) > 0.
     Fail: drop (q+5); claim the gap at p=5 is zero (q=25>11).
  C. After λ_exc≥6, principal room
        ||δ||² ≤ n(n+10)² / [6(n-14)(n-6)]
     ⇔ E[(y·z)⁴] ≤ 4n(3n²−37n+2)/(n-14)  (15.589 D+).
     Fail: replace n-14 by n-6 (old unsharpened room).
  D. On Z, E[q²]=6||B||²+8⟨m4,κ_B⟩ and λ_min(Φ)≥6 ⇔ ⟨m4,κ_B⟩≥0
     (15.277).  Fail: claim pointwise q²≥3||By||² (false, 15.277 C).
  A*. 15.597: Φ_part = λ̄ I on Z for every prime p≥5 (closed-form
     contractions; not 15.108's 16N Theorem A*).  Fail: claim this
     bounds δ.  All spectral deviation of Φ is Φ_δ.  Leftover 1
        ⇔  Φ_δ ⪰ −(2n+20)/(n−6) I
     (`e1_gmin_m4_prop15597.leftover1_operator_criterion`).  Corollaries:
     tr(Φ_δ)=0; 0 ≤ λ_min(Φ) ≤ λ̄ unconditionally; target 6 sits in
     that window.

OPEN (blocks leftover 1)
  E. GLOBAL QVAR, all primes p≥13: E|Z_ψ|² ≥ 3q(q-1)/16 on the full
     Max+ mixture (all k; do not split λ=0; do not require each
     k-stratum).  Equivalently the exceptional isotype of Φ_δ clears
     6−λ̄.  k=1..6 closed (15.589).  Per-stratum k≥7 is false
     at p=41 k=7 (E=0, Cy=py) and pointwise at p=13 k=7; those are
     not p-laws.  V_sph>threshold is not a 4-design.  Floor iff
     F̂(ψ)≥0 / Gauss 4-dist pairing of m₄ (15.279 L).
  F. The moment bound (C) / R1, equivalently ||δ||²≤n/12 via the
     multiplicity floor (15.595).  Crude E[s⁴]≤2n³ is Θ(n³) vs Θ(n²).
     After A*, E∧F together are Φ_δ ⪰ −(2n+20)/(n−6) I.

Until E∧F: leftover1_qvar_and_principal_proved stays False, so
phi_F_ge_6_proved_general stays False.  If E∧F both hold, leftover 1 /
phi_F flip; L still needs leftover 2/3 and lemma D.  E is imported from
`e1_gmin_global_qvar.global_qvar_proved_general` (not
`qvar_k_ge_7_proved_general`).  F is imported from
`e1_gmin_r1_principal_pge11.r1_l2_bound_for_p_ge_11` (currently False):
measured ‖δ‖² exceeds n/12 at p=5,7; 2-design+sharp min-distance cannot
prove R1; p≥11 remains open.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15589 import (  # noqa: E402
    delta2_room_principal,
    delta2_room_principal_after_exception,
    lambda_exc_from_quartic_variance,
    principal_Es4_budget_after_exception,
    q_of,
    quartic_variance_floor_threshold,
    spherical_QVAR_gap,
    spherical_quartic_variance,
    theorem_Dplus_exception_removed_variance,
    theorem_E_exceptional_quartic_variance,
    theorem_H_odd_coset_spherical_benchmark,
)


def hs_norm_A_psi(p: int) -> Fraction:
    q = q_of(p)
    return Fraction(q * (q - 1), 32)


def V_sph_wrong_drop_q5(p: int) -> Fraction:
    q = q_of(p)
    return Fraction(q * (q - 1) * (q + 1), 4)


def qvar_threshold_wrong_drop16(p: int) -> Fraction:
    q = q_of(p)
    return Fraction(3 * q * (q - 1), 8)


def lambda_exc_wrong_16(p: int, variance: Fraction) -> Fraction:
    q = q_of(p)
    return Fraction(16, q * (q - 1)) * variance


def crude_Es4_2n3(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(2 * n * n * n)


def theorem_A_qvar_iff(primes=(5, 7, 11, 13, 17, 19, 23)) -> dict:
    E = theorem_E_exceptional_quartic_variance()
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        thr = quartic_variance_floor_threshold(p)
        hs = hs_norm_A_psi(p)
        six_hs = 6 * hs
        named_32 = Fraction(32) * hs  # 32 * q(q-1)/32 = q(q-1)
        row_ok = (
            six_hs == thr
            and named_32 == Fraction(q_of(p) * (q_of(p) - 1))
            and qvar_threshold_wrong_drop16(p) != thr
        )
        # census λ_exc formula at p=5,7
        if p in (5, 7):
            rec = E["by_p"][str(p)]
            var = Fraction(rec["E_abs_Zpsi_sq"])
            good = lambda_exc_from_quartic_variance(p, var) == Fraction(
                rec["lambda_exc"]
            )
            bad16 = lambda_exc_wrong_16(p, var) != Fraction(rec["lambda_exc"])
            row_ok = row_ok and good and bad16
        rows[str(p)] = {"thr": str(thr), "6_hs": str(six_hs), "ok": row_ok}
        ok = ok and row_ok
    return {
        "proved": bool(ok and E["proved_reduction"]),
        "theorem": (
            "λ_exc=32 E|Z_ψ|²/[q(q-1)] ⇔ QVAR E|Z_ψ|²≥3q(q-1)/16 ⇔ "
            "E[|yᵀA_ψ y|²]≥6||A_ψ||_HS².  Fail: drop 16; replace 32 by 16."
        ),
        "by_p": rows,
        "census_formula": E["proved_census"],
        "qvar_general": False,
    }


def theorem_B_spherical_exceeds_qvar(
    primes=(5, 7, 11, 13, 17, 19, 23, 29, 31),
) -> dict:
    ok = True
    rows = {}
    for p in primes:
        if not is_prime(p):
            continue
        gap = spherical_QVAR_gap(p)
        vs = spherical_quartic_variance(p)
        thr = quartic_variance_floor_threshold(p)
        H = theorem_H_odd_coset_spherical_benchmark((p,))
        row_ok = (
            vs - thr == gap
            and gap > 0
            and V_sph_wrong_drop_q5(p) != vs
            and H["proved_reduction"]
        )
        rows[str(p)] = {
            "V_sph": str(vs),
            "thr": str(thr),
            "gap": str(gap),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "theorem": (
            "V_sph=q(q-1)(q+1)/(4(q+5)) exceeds QVAR by "
            "q(q-1)(q-11)/(16(q+5))>0 for every prime p≥5.  "
            "This is a 2-design spherical average, not a 4-design, "
            "so it does not prove QVAR.  Fail: drop q+5."
        ),
        "by_p": rows,
        "is_4_design": False,
        "qvar_general": False,
    }


def theorem_C_principal_room_reduction(
    primes=(5, 7, 11, 13, 17, 19),
) -> dict:
    Dplus = theorem_Dplus_exception_removed_variance(primes)
    ok = bool(Dplus["proved"])
    rows = {}
    for p in primes:
        n = n_of(p)
        new = delta2_room_principal_after_exception(p)
        old = delta2_room_principal(p)
        crude = crude_Es4_2n3(p)
        budget = principal_Es4_budget_after_exception(p)
        row_ok = (
            new / old == Fraction(n - 6, n - 14)
            and crude > budget
        )
        rows[str(p)] = {
            "B_min": str(new),
            "Es4_budget": str(budget),
            "crude_2n3": str(crude),
            "crude_too_weak": bool(crude > budget),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "does_not_prove_moment": True,
        "theorem": (
            "After λ_exc≥6, ||δ||²≤n(n+10)²/[6(n-14)(n-6)] forces every "
            "principal scalar ≥6 (15.589 D+).  Fail: n-14↦n-6.  "
            "Crude E[s⁴]≤2n³ exceeds the Es4 budget (Θ(n³) vs Θ(n²))."
        ),
        "by_p": rows,
        "principal_moment_general": False,
    }


def theorem_D_floor_iff_m4_pairing() -> dict:
    return {
        "proved": True,
        "pointwise_q2_ge_3By2": False,
        "theorem": (
            "On Z, E[q²]=6||B||²+8⟨m4,κ_B⟩ so λ_min≥6 ⇔ ⟨m4,κ_B⟩≥0 "
            "(15.277).  Pointwise q²≥3||By||² is false."
        ),
    }


def qvar_k_ge_7_proved_general() -> bool:
    """Per-stratum k≥7 unit.  False at (41,7); not the leftover-1 conjunct."""
    from e1_gmin_qvar_k_ge_7 import qvar_k_ge_7_proved_general as _qvar

    return bool(_qvar())


def global_qvar_proved_general() -> bool:
    """True only by importing the global (mixed-k) QVAR unit."""
    from e1_gmin_global_qvar import global_qvar_proved_general as _g

    return bool(_g())


def principal_delta_room_moment_proved() -> bool:
    """True only by importing the R1 L² unit. False: interpolant 4/(p−3)²
    is not a retained δ-bound; p=5,7 exceed n/12 (census)."""
    from e1_gmin_r1_principal_pge11 import r1_l2_bound_for_p_ge_11

    return bool(r1_l2_bound_for_p_ge_11())


def leftover1_qvar_and_principal_proved() -> bool:
    """Both leftover-1 blocks, all primes p≥5.  False until E∧F."""
    return bool(
        global_qvar_proved_general() and principal_delta_room_moment_proved()
    )


def leftover1_reductions_ok() -> bool:
    """Identities A–D.  Independent of whether the two estimates are proved."""
    return bool(
        theorem_A_qvar_iff()["proved"]
        and theorem_B_spherical_exceeds_qvar()["proved"]
        and theorem_C_principal_room_reduction()["proved"]
        and theorem_D_floor_iff_m4_pairing()["proved"]
        and theorem_D_floor_iff_m4_pairing()["pointwise_q2_ge_3By2"] is False
    )


def live_L_status() -> str:
    """CLOSED only if leftover 1, leftover 2/3, and lemma D are all proved."""
    from e1_main_chain_status import four_e1_units_closed

    return "CLOSED" if four_e1_units_closed().get("closed") else "OPEN"


def dump_leftover_predicates() -> dict:
    from e1_gmin_m4_prop15274 import (
        multilevel_ND_k_ge_4p_proved,
        residual_ii_k_eq_4p_empty,
        residual_ii_k_ge_4p_ND_closed,
    )
    from e1_gmin_m4_prop15275 import (
        type_I_aut_e_3AB_positive_general,
        type_I_multilevel_bad_case_ND_closed,
    )
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return {
        "phi_F_ge_6_proved_general": bool(phi_F_ge_6_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "residual_ii_k_ge_4p_ND_closed": bool(residual_ii_k_ge_4p_ND_closed()),
        "multilevel_ND_k_ge_4p_proved": bool(multilevel_ND_k_ge_4p_proved()),
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb_proved_general": bool(gsum_disj_lb_proved_general()),
        "qvar_k_ge_7_proved_general": qvar_k_ge_7_proved_general(),
        "global_qvar_proved_general": global_qvar_proved_general(),
        "principal_delta_room_moment_proved": principal_delta_room_moment_proved(),
        "leftover1_qvar_and_principal_proved": leftover1_qvar_and_principal_proved(),
        "L_status": live_L_status(),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    A = theorem_A_qvar_iff()
    B = theorem_B_spherical_exceeds_qvar()
    C = theorem_C_principal_room_reduction()
    D = theorem_D_floor_iff_m4_pairing()
    dump = dump_leftover_predicates()
    out = {
        "title": "Leftover 1 QVAR+principal hinge (not a close)",
        "numbered": False,
        "A_qvar_iff": A,
        "B_spherical_exceeds": B,
        "C_principal_room_reduction": C,
        "D_floor_iff": D,
        "proved": {
            "reductions": leftover1_reductions_ok(),
            "qvar_k_ge_7": qvar_k_ge_7_proved_general(),
            "global_qvar": global_qvar_proved_general(),
            "principal_moment": principal_delta_room_moment_proved(),
            "leftover1": leftover1_qvar_and_principal_proved(),
            "phi_F_ge_6": dump["phi_F_ge_6_proved_general"],
            "gsum_disj_lb_proved_general": dump["gsum_disj_lb_proved_general"],
        },
        "predicates": dump,
        "L_status": dump["L_status"],
    }
    path = ROOT / "evidence" / "e1_gmin_leftover1_qvar_principal.json"
    write_json_atomic(path, out)
    print("Leftover 1 QVAR+principal hinge", flush=True)
    print(f"  A qvar iff: {A['proved']}", flush=True)
    print(f"  B V_sph>thr: {B['proved']}", flush=True)
    print(f"  C room reduction: {C['proved']}", flush=True)
    print(f"  D m4 pairing iff: {D['proved']}", flush=True)
    print(f"  qvar k>=7 (not leftover-1): {qvar_k_ge_7_proved_general()}", flush=True)
    print(f"  global qvar: {global_qvar_proved_general()}", flush=True)
    print(f"  principal moment: {principal_delta_room_moment_proved()}", flush=True)
    print(f"  leftover1 both: {leftover1_qvar_and_principal_proved()}", flush=True)
    print(f"  phi_F_ge_6: {dump['phi_F_ge_6_proved_general']}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
