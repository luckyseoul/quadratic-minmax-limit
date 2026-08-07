"""Tests for Prop 15.215 — K4/Wick Tκ–η_part algebra; residual (i) OPEN."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
)
from e1_gmin_m4_prop15198 import eta_sumsq_budget_for_wick_hi  # noqa: E402
from e1_gmin_m4_prop15215 import (  # noqa: E402
    K4_bound_proved_general,
    eta_part_norm_sq,
    hinge_status_215,
    residual_i_closed_via_215,
    theorem_T3_kappa,
    theorem_Tkappa_norm,
    theorem_Tphi_eq_nplus2_star,
    theorem_eta_part_formula,
    theorem_eta_part_lt_wick_budget,
    theorem_n3_count,
    theorem_reduction_K4_via_eta,
    theorem_star_values_by_exhaustion,
    theorem_sum_kappa_sq,
)


def test_star_exhaustion():
    t = theorem_star_values_by_exhaustion()
    assert t["proved"]
    assert t["n_kappa1_labels"] == 48
    assert t["n_kappa3_labels"] == 16


def test_sum_kappa_sq():
    assert theorem_sum_kappa_sq()["proved"]


def test_n3_count():
    t = theorem_n3_count()
    assert t["proved"]
    # census anchors
    assert t["by_p_sample"]["3"]["n3"] == 30
    assert t["by_p_sample"]["5"]["n3"] == 3250
    assert t["by_p_sample"]["7"]["n3"] == 53900


def test_Tkappa_norm():
    assert theorem_Tkappa_norm()["proved"]
    for p, expected in [(3, 17280), (5, 1872000), (7, 31046400)]:
        n = n_of(p)
        assert 6 * n * (n - 1) * (n - 2) * (n - 6) == expected


def test_Tphi():
    assert theorem_Tphi_eq_nplus2_star()["proved"]


def test_T3_kappa():
    assert theorem_T3_kappa()["proved"]


def test_eta_part_formula():
    t = theorem_eta_part_formula()
    assert t["proved"]
    assert eta_part_norm_sq(5) == Fraction(2912, 100)  # 29.12 = 728/25?
    # 5*24*26*28/(6*25*20) = 5*24*26*28/3000 = 3360/3000*26/26 wait
    assert eta_part_norm_sq(5) == Fraction(5 * 24 * 26 * 28, 6 * 25 * 20)


def test_eta_part_lt_budget():
    t = theorem_eta_part_lt_wick_budget()
    assert t["proved"]
    for p in (5, 7, 11, 13, 17, 19):
        assert eta_part_norm_sq(p) < eta_sumsq_budget_for_wick_hi(p)


def test_reduction_open():
    r = theorem_reduction_K4_via_eta()
    assert r["proved_implications"]
    assert r["K4_bound_proved_general"] is False


def test_predicates_false():
    assert K4_bound_proved_general() is False
    assert residual_i_closed_via_215() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_215()
    assert h["eta_part_lt_wick_budget"]
    assert h["delta_E4p_bound"] is False
    assert h["residual_i_closed"] is False
