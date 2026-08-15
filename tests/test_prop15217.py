"""Tests for Prop 15.217 — Φ identity; m4₂ residual-(i) reduction."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15217 import (  # noqa: E402
    Phi_closed,
    delta_room_for_R,
    hinge_status_217,
    ke2_of,
    m4_target_R,
    residual_i_closed_via_217,
    theorem_Phi_identity,
    theorem_R_implies_K4_thr,
    theorem_R_le_2p_iff_m4_bound,
    theorem_eta_formula_with_Phi,
    theorem_m4_bound_open,
)


def test_Phi_identity():
    t = theorem_Phi_identity()
    assert t["proved"]
    assert Phi_closed(5) == 1950
    assert Phi_closed(7) == 14700
    # n(n-1)(n-2)/8
    assert Phi_closed(5) == Fraction(26 * 25 * 24, 8)


def test_R_iff_m4_bound():
    t = theorem_R_le_2p_iff_m4_bound()
    assert t["proved"]
    for p in (5, 7, 11, 13, 17, 19):
        assert 2 * Phi_closed(p) / Fraction(p * p) == m4_target_R(p)
        assert delta_room_for_R(p) == m4_target_R(p) - ke2_of(p)
        assert delta_room_for_R(p) > 0


def test_eta_formula():
    assert theorem_eta_formula_with_Phi()["proved"]


def test_conditional_R_implies_dual_eq():
    t = theorem_R_implies_K4_thr()
    assert t["proved"]  # the implication algebra
    assert "OPEN" in t["conditional_on"] or "m₄" in t["conditional_on"]


def test_m4_bound_still_open():
    t = theorem_m4_bound_open()
    assert t["proved"] is False
    assert t["K4_equivalence_algebra_ok"]
    assert t["pathc_residual_sufficient"]


def test_residual_i_217_path_open_live_dual_eq_closed():
    assert residual_i_closed_via_217() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_217()
    assert h["Phi_identity"] is True
    assert h["R_iff_m4_bound"] is True
    assert h["m4_bound"] is False
    assert h["residual_i_closed_via_217"] is False
