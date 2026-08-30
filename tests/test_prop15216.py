"""Tests for Prop 15.216 — residual (i) K₄ thr path; R≤2p OPEN (fatal gap)."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15216 import (  # noqa: E402
    R_ke,
    crude_eta_ub,
    hinge_status_216,
    residual_i_closed_via_216,
    theorem_K4_le_thr,
    theorem_Q10_blocks_dual_eq_all_p,
    theorem_Q_le_10,
    theorem_R_ke_formula,
    theorem_R_ke_le_2p,
    theorem_R_m4_le_2p,
    theorem_crude_le_thr_eta_budget,
    theorem_dual_eq_empty_all_p,
    theorem_eta_le_crude_when_R_le_2p,
    theorem_kappa_spectral_support,
    theorem_residual_i_closed,
    theorem_star_kappa_orthogonal,
    thr_eta_budget,
)


def test_star_kappa_orthogonal():
    assert theorem_star_kappa_orthogonal()["proved"]


def test_kappa_spectral_support():
    assert theorem_kappa_spectral_support()["proved"]


def test_R_ke_formula():
    t = theorem_R_ke_formula()
    assert t["proved"]
    assert R_ke(5) == Fraction(160, 23)
    assert R_ke(7) == Fraction(224, 41)


def test_R_ke_le_2p():
    assert theorem_R_ke_le_2p()["proved"]
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert R_ke(p) <= 2 * p


def test_R_m4_le_2p_open_fatal_gap():
    """Convex combo of R_ke and 4p does not yield R≤2p — AI-test BLOCK."""
    t = theorem_R_m4_le_2p()
    assert t["proved"] is False
    assert "fatal_gap" in t
    assert "4p" in t["fatal_gap"]


def test_eta_crude_conditional_and_budget():
    # Implication R≤2p ⇒ η≤crude is conditional-true; crude≤thr proved.
    assert theorem_eta_le_crude_when_R_le_2p()["proved"]
    assert theorem_eta_le_crude_when_R_le_2p().get("conditional") is True
    assert theorem_crude_le_thr_eta_budget()["proved"]
    for p in (5, 7, 11, 13, 17, 19):
        assert crude_eta_ub(p) <= thr_eta_budget(p)
        assert thr_eta_budget(p) - crude_eta_ub(p) == Fraction(
            (p * p + 1) * (p * p + 9), 24
        )


def test_K4_and_Q_open():
    assert theorem_K4_le_thr()["proved"] is False
    assert theorem_Q_le_10()["proved"] is False


def test_Q10_blocks_all_p():
    """Q≤10 ⇒ dual-eq empty is proved independently of K₄ thr."""
    t = theorem_Q10_blocks_dual_eq_all_p()
    assert t["proved"]
    for p in (5, 7, 11, 13):
        assert t["by_p_sample"][str(p)]["blocks"]


def test_dual_eq_empty_open():
    assert theorem_dual_eq_empty_all_p()["proved"] is False


def test_k4_path_still_open():
    t = theorem_residual_i_closed()
    assert t["proved"] is False
    assert theorem_K4_le_thr()["proved"] is False


def test_dual_eq_via_gplus():
    assert residual_i_dual_eq_empty_proved_general() is True
    assert residual_i_closed_via_216() is True


def test_type_I_closed_not_via_gsum():
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True


def test_e1_closed():
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_216()
    assert h["residual_i_closed_via_216"] is True
    assert h["Q_le_10"] is False
    assert h["R_m4_le_2p"] is False
    assert h["e1_closed_general"] is True
    assert h["L_open"] is False
