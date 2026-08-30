#!/usr/bin/env python3
"""
Prop 15.193 — Residual (ii) exhaustiveness audit (honest).

Auditor request: under residual-(ii) hypotheses, does failure of weak ND
force both (1) S∈{2,4} on every Max+ vector and (2) f_e=3−S everywhere?

VERDICT (this prop): **exhaustiveness is NOT proved.**
  Affine dual two-level branch is CLOSED by Prop 15.179 freeze-to-tight.
  The bounded even range through 4p−2 is CLOSED by separate ND: (ii-a)
  15.237 + (ii-b) 15.236 + affine 15.179.  The even-k≥4p multi-level
  remainder is still OPEN (15.274).

Does **not** soft-close residual (i) / E1 / L. Soft-close forbidden.

SETUP (residual (ii): deep freeness-fail s₊=2, k≥3p)
  s₊=2: min Max+ score is 2; scores even when k even.
  Freeness-fail for edge e: f_e ≡ +1 on U₂ := {y∈Max+: S(y)=2}.
  Weak ND Max+-witness: exists y with S=2, f_e=−1 ⇒ S_H=1 ⇒ Φ(H)≥Φ−2.
  So freeness-fail ⇒ that Max+ witness of weak ND is unavailable.

WHAT IS PROVED (score algebra, Max+-free Fraction):
  A. Freeness-fail ⇔ f_e ≡ +1 on U₂ (definition).
  B. Under freeness-fail and s₊=2: on every Max+ vector,
        S_H = S + f_e ≥ 3
     (on U₂: S_H=3; on S≥4: S_H≥S−1≥3). Hence no Max+ vector has S_H≤1,
     so the 15.171.C Max+ weak-ND witness fails.
  C. Two-level mass lower bound: if scores ≥2 step 2, then
        N₂/N ≥ 2 − k/(2p)
     (from E[S]≥4−2a). Auto-freeness when lb > thr=(p+1)/(2p), i.e. k≤3p−2
     (15.168.F). For k≥3p, lb ≤ 1/2 < thr, so freeness *can* fail.
  D. Affine dual two-level freeness-fail (S∈{2,4} and f_e=3−S):
     S_H≡3 ⇒ k=3p−1, impossible for k≥3p (15.179 freeze). Branch CLOSED.

WHAT IS NOT PROVED (exhaustiveness still False; bounded ND closed elsewhere):
  E. Exhaustiveness. Freeness-fail does **not** by itself force
     support(S)⊆{2,4} on Max+. Multi-level in the bounded even range is
     closed by ND in 15.237, not by forcing two-level; k≥4p remains open.
  F. Affine completion. Even under two-level S∈{2,4} and freeness-fail
     (f_e≡+1 on U₂), f_e need not be ≡−1 on U₄. Non-affine two-level is
     closed by ND in 15.236 within its stated bounded range, not by freeze.
  G. Prop 15.171 dual two-level Gsum Farkas assumes affine f_e=3−S on
     both Max±; it is not needed after 15.236/15.237.

PREDICATES
  residual_ii_affine_branch_closed()  = True   (15.179)
  residual_ii_exhaustiveness_proved() = False  (this prop)
  residual_ii_bounded_even_k_le_4p_minus_2_closed()
                                      = affine ∧ (ii-a) ∧ (ii-b)
                                        (15.179 + 15.237 + 15.236)
  residual_ii_full_closed()           = bounded result ∧ live even-k≥4p gate
  deep_s2_freeness_fail_k_ge_3p_ND_closed() wires to residual_ii_full_closed

Residual (ii) subcases for k≥3p freeness-fail:
  (ii-a) multi-level: CLOSED by 15.237
  (ii-b) two-level non-affine: CLOSED by 15.236
  Affine two-level (ii-c): CLOSED by 15.179

Writes evidence/e1_gmin_m4_prop15193.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15720 import (  # noqa: E402
    required_bitight_levels_empty_all_primes,
)
from e1_gmin_m4_prop15168 import freeness_threshold, deep_s2_freeness_lb  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15179 import (  # noqa: E402
    residual_ii_dual_twolevel_affine_closed,
    freeze_S_H_equals_3,
    theorem_freeze_all_primes,
)
from e1_gmin_m4_prop15236 import (  # noqa: E402
    residual_ii_a_ND_closed,
    residual_ii_b_ND_closed,
)


# ---------------------------------------------------------------------------
# A–B: freeness-fail forces Max+ weak-ND witness failure (proved)
# ---------------------------------------------------------------------------


def freeness_fail_forces_f_e_plus_on_U2() -> dict:
    """Definition restated: freeness-fail ⇔ f_e≡+1 on U₂."""
    return {
        "proved": True,
        "theorem": (
            "Deep freeness-fail for edge e means f_e≡+1 on U₂={S=2}⊆Max+. "
            "Equivalently: no Max+ vector with S=2 and f_e=−1."
        ),
    }


def freeness_fail_kills_Max_plus_weak_ND_witness() -> dict:
    """
    Under freeness-fail + s₊=2: S_H ≥ 3 on every Max+ vector,
    so no Max+ y with Q_H(y)=Φ−2 via S_H=1.
    """
    # Check discrete score algebra for S≥2 even, f_e=±1 under freeness-fail
    samples = []
    ok = True
    for S in (2, 4, 6, 8, 10):
        if S == 2:
            f_e = 1  # freeness-fail
            S_H = S + f_e
            assert S_H == 3
            samples.append({"S": S, "f_e": f_e, "S_H": S_H, "ge_3": S_H >= 3})
            if S_H < 3:
                ok = False
        else:
            # freeness-fail constrains only U₂; on S≥4 both signs possible
            for f_e in (-1, 1):
                S_H = S + f_e
                samples.append({"S": S, "f_e": f_e, "S_H": S_H, "ge_3": S_H >= 3})
                if S_H < 3:
                    ok = False
    return {
        "proved": ok,
        "samples": samples,
        "theorem": (
            "Freeness-fail + s₊=2 ⇒ on Max+: S_H=3 on U₂ and S_H≥3 on S≥4. "
            "Hence no Max+ vector has S_H≤1; the 15.171.C Max+ weak-ND "
            "witness (S=2,f_e=−1⇒S_H=1) is unavailable."
        ),
    }


# ---------------------------------------------------------------------------
# C: mass lower bound / freeness-fail possible for k≥3p (proved algebra)
# ---------------------------------------------------------------------------


def freeness_fail_possible_when_k_ge_3p(p: int) -> dict:
    """For two-level mass lb: k≥3p ⇒ N₂/N lb ≤ 1/2 < thr ⇒ freeness can fail."""
    thr = freeness_threshold(p)
    rows = []
    ok = True
    for k in (3 * p - 2, 3 * p - 1, 3 * p, 3 * p + 1, 3 * p + 2):
        lb = deep_s2_freeness_lb(p, k)
        can_fail = lb <= thr
        auto = lb > thr
        rows.append(
            {
                "k": k,
                "N2_N_lb": str(lb),
                "thr": str(thr),
                "auto_freeness": auto,
                "freeness_can_fail": can_fail,
            }
        )
        if k <= 3 * p - 2 and not auto:
            ok = False
        if k >= 3 * p and not can_fail:
            ok = False
    return {
        "p": p,
        "proved": ok,
        "rows": rows,
        "theorem": (
            "N₂/N ≥ 2−k/(2p) under s₊=2 even scores. Auto-freeness iff "
            "lb>thr, i.e. k≤3p−2. For k≥3p, lb≤1/2<(p+1)/(2p)=thr so "
            "freeness-fail is not mass-obstructed."
        ),
    }


# ---------------------------------------------------------------------------
# E–F: exhaustiveness gaps (explicitly NOT proved)
# ---------------------------------------------------------------------------


def multi_level_not_ruled_out() -> dict:
    """
    Freeness-fail does not force support ⊆ {2,4}.
    Multi-level configurations remain residual-(ii) candidates.
    """
    # Illustrative: scores {2,4,6} with masses a,b,c; E[S]=2a+4b+6c = k/p
    # and a ≤ thr can hold with c>0 (not two-level).
    examples = []
    for p, k in ((5, 16), (5, 18), (7, 22), (7, 24)):
        # k even ≥ 3p
        if k < 3 * p or k % 2 == 1:
            continue
        thr = freeness_threshold(p)
        target = Fraction(k, p)
        # pick c=1/10 mass at 6, b residual two-level-like
        c = Fraction(1, 10)
        # 2a+4b+6c = target, a+b+c=1 ⇒ 2a+4(1-a-c)+6c = target
        # 2a + 4 − 4a − 4c + 6c = target ⇒ 4 − 2a + 2c = target
        # 2a = 4 + 2c − target ⇒ a = 2 + c − target/2
        a = 2 + c - target / 2
        b = 1 - a - c
        feasible = a >= 0 and b >= 0 and c > 0 and a <= thr
        examples.append(
            {
                "p": p,
                "k": k,
                "mass_2": str(a),
                "mass_4": str(b),
                "mass_6": str(c),
                "E_S": str(target),
                "a_le_thr": a <= thr,
                "multi_level_mass_feasible": feasible,
            }
        )
    any_feasible = any(e["multi_level_mass_feasible"] for e in examples)
    return {
        "exhaustiveness_proved": False,
        "multi_level_mass_feasible_examples_exist": any_feasible,
        "examples": examples,
        "theorem": (
            "Score first-moment + freeness thr do not force two-level: "
            "there exist (p,k≥3p) with positive mass at S=6, N₂/N≤thr, "
            "and E[S]=k/p. Multi-level freeness-fail is not mass-impossible. "
            "No shipped score inequality kills (ii-a)."
        ),
    }


def non_affine_two_level_not_ruled_out() -> dict:
    """
    Two-level S∈{2,4} + f_e≡+1 on U₂ does not force f_e≡−1 on U₄.
    Freeze requires full affine f_e=3−S.
    """
    # On two-level: a = 2 − k/(2p). Freeness-fail: f_e=1 on S=2.
    # Non-affine: let mass of (S=4, f_e=+1) be positive.
    # Then S_H is 3 on U₂, and {3,5} on U₄ — not constantly 3.
    # First moment E[S_H]=(k+1)/p still holds; freeze contradiction needs
    # S_H≡3, which needs f_e=−1 on all U₄.
    freeze = freeze_S_H_equals_3()
    return {
        "exhaustiveness_proved": False,
        "affine_forces_S_H_const_3": freeze["proved"],
        "non_affine_S_H_not_const": True,
        "theorem": (
            "Under S∈{2,4} and freeness-fail alone: f_e≡+1 on U₂ only. "
            "If some y with S=4 has f_e=+1 then S_H(y)=5 ≠ 3, so S_H is "
            "not constantly 3 and 15.179 freeze does not fire. Non-affine "
            "two-level freeness-fail (ii-b) is not closed by freeze."
        ),
    }


def desired_exhaustiveness_lemma_statement() -> dict:
    """
    The lemma that would finish residual (ii) with 15.179 — still unproved.
    """
    return {
        "statement": (
            "Exhaustiveness lemma (desired, NOT proved): Under residual-(ii) "
            "hypotheses (deep s₊=2 undercutter, freeness-fail for e, k≥3p), "
            "failure of Max+ weak ND implies (1) S∈{2,4} on every Max+ "
            "vector and (2) f_e=3−S on every Max+ vector "
            "(and the dual Max− affine form if needed for dual two-level)."
        ),
        "proved": False,
        "if_proved_then": (
            "15.179 freeze closes residual (ii) fully for k≥3p; fail-eq "
            "k=3p−1 already empty under bi-tight."
        ),
        "current_status": (
            "Only (partial) freeness-fail ⇒ f_e≡+1 on U₂ is proved. "
            "(1) and (2) are open; see multi_level_not_ruled_out and "
            "non_affine_two_level_not_ruled_out."
        ),
    }


# ---------------------------------------------------------------------------
# Predicates
# ---------------------------------------------------------------------------


def residual_ii_exhaustiveness_proved() -> bool:
    """Full exhaustiveness weak-ND-fail ⇒ affine two-level — not proved."""
    return False


def residual_ii_affine_branch_closed() -> bool:
    """15.179 dual two-level affine freeness-fail closed for k≥3p."""
    return residual_ii_dual_twolevel_affine_closed()


def residual_ii_bounded_even_k_le_4p_minus_2_closed() -> bool:
    """
    Historical bounded split: affine branch plus ND for (ii-a) and (ii-b)
    in the even range 3p+1 <= k <= 4p-2.
    """
    return bool(
        residual_ii_affine_branch_closed()
        and residual_ii_a_ND_closed()
        and residual_ii_b_ND_closed()
    )


def residual_ii_full_closed() -> bool:
    """Full residual-(ii) gate, including the live even-k>=4p remainder."""
    from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed

    return bool(
        residual_ii_bounded_even_k_le_4p_minus_2_closed()
        and residual_ii_k_ge_4p_ND_closed()
    )


def residual_ii_open_subcases() -> list[str]:
    out = []
    if not residual_ii_a_ND_closed():
        out.append(
            "(ii-a) multi-level freeness-fail s₊=2, k≥3p "
            "(S support not forced ⊆{2,4}; no two-level slope)"
        )
    if not residual_ii_b_ND_closed():
        out.append(
            "(ii-b) two-level S∈{2,4} non-affine freeness-fail "
            "(f_e≡+1 on U₂, not f_e≡−1 on U₄)"
        )
    if not residual_ii_affine_branch_closed():
        out.append("(ii-c) dual two-level affine freeness-fail (should be closed by 15.179)")
    from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed

    if not residual_ii_k_ge_4p_ND_closed():
        out.append("non-Walsh residual (ii), even k≥4p: multi-level leftover")
    return out


def hinge_status_193() -> dict:
    from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed

    return {
        "residual_ii_affine_branch_k_ge_3p": (
            "CLOSED by 15.179 freeze-to-tight first moment"
            if residual_ii_affine_branch_closed()
            else "OPEN"
        ),
        "residual_ii_b_ND": (
            "CLOSED by 15.236 (dichotomy + dual-bad-case empty)"
            if residual_ii_b_ND_closed()
            else "OPEN"
        ),
        "residual_ii_a_ND": (
            "CLOSED"
            if residual_ii_a_ND_closed()
            else "OPEN — multi-level"
        ),
        "residual_ii_exhaustiveness": (
            "OPEN — freeness-fail does not force S∈{2,4} and f_e=3−S"
        ),
        "residual_ii_full": (
            "CLOSED" if residual_ii_full_closed() else
            "OPEN: bounded k≤4p−2 is closed; multi-level even k≥4p remains"
        ),
        "type_I_multilevel": (
            "CLOSED"
            if type_I_multilevel_bad_case_ND_closed()
            else "OPEN — multi-level Max− bad case"
        ),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "bound_proved_general": False,
        "open_subcases": residual_ii_open_subcases(),
        "open": (
            "Close the multi-level even-k≥4p remainder (15.274). The bounded "
            "range through 4p−2 is already closed by 15.179/236/237."
        ),
    }


def main() -> dict:
    freeze_th = theorem_freeze_all_primes([p for p in range(5, 40) if is_prime(p)])
    A = freeness_fail_forces_f_e_plus_on_U2()
    B = freeness_fail_kills_Max_plus_weak_ND_witness()
    C_rows = {
        str(p): freeness_fail_possible_when_k_ge_3p(p) for p in (5, 7, 11, 13, 17)
    }
    C_ok = all(r["proved"] for r in C_rows.values())
    E = multi_level_not_ruled_out()
    F = non_affine_two_level_not_ruled_out()
    lemma = desired_exhaustiveness_lemma_statement()
    affine = residual_ii_affine_branch_closed()
    exhaust = residual_ii_exhaustiveness_proved()
    full = residual_ii_full_closed()
    bt = required_bitight_levels_empty_all_primes()
    out = {
        "title": (
            "Prop 15.193 residual (ii) exhaustiveness audit: affine branch "
            "CLOSED (15.179); full residual (ii) OPEN at even k≥4p"
        ),
        "proved": {
            "freeness_fail_eq_f_e_plus_on_U2": A["proved"],
            "freeness_fail_kills_Max_plus_weak_ND_witness": B["proved"],
            "freeness_fail_mass_possible_k_ge_3p": C_ok,
            "affine_dual_twolevel_freeze_15_179": freeze_th["proved"] and affine,
            "residual_ii_affine_branch_closed": affine,
            "residual_ii_exhaustiveness_proved": exhaust,
            "residual_ii_b_ND_closed": residual_ii_b_ND_closed(),
            "residual_ii_a_ND_closed": residual_ii_a_ND_closed(),
            "residual_ii_full_closed": full,
            "multi_level_ruled_out": False,
            "non_affine_two_level_ND_closed": residual_ii_b_ND_closed(),
            "non_affine_two_level_exhaustiveness": False,
            "residual_i_closed": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
            "bi_tight_empty": bt,
        },
        "algebra": {
            "A_freeness_def": A,
            "B_weak_ND_witness": B,
            "C_mass": C_rows,
            "E_multi_level": E,
            "F_non_affine": F,
            "desired_lemma": lemma,
            "freeze_15_179": freeze_th,
        },
        "hinge": hinge_status_193(),
        "F3": (
            "Do not claim residual (ii) full CLOSED from the bounded 15.179/236/237 "
            "split; even k≥4p and E1/L remain open."
        ),
        "audit_verdict": (
            "Residual (ii): affine and bounded even range closed; "
            "multi-level even k≥4p remains open."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15193.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.193 residual (ii) exhaustiveness audit")
    print(f"  freeness-fail ⇒ f_e≡+1 on U₂: {A['proved']}")
    print(f"  freeness-fail kills Max+ weak ND witness: {B['proved']}")
    print(f"  mass allows freeness-fail k≥3p: {C_ok}")
    print(f"  affine branch closed (15.179): {affine}")
    print(f"  exhaustiveness proved: {exhaust}")
    print(f"  residual (ii) FULL closed: {full}")
    print(f"  open subcases: {residual_ii_open_subcases()}")
    print(f"  multi-level mass feasible examples: {E['multi_level_mass_feasible_examples_exist']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
