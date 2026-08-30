"""Tests for Prop 15.213 — Veronese/λ_* attack; residual (i) still OPEN."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15213 import (  # noqa: E402
    certify_cost_D_sum_ne_ratio,
    certify_lambda_star_floor_small_p,
    free_e_bound_proved_general,
    hinge_status_213,
    lambda_star_floor_proved_general,
    residual_i_closed_via_213,
    theorem_dual_four_thirds_budget,
    theorem_four_design_shortcut_dead,
    theorem_lambda_star_sufficiency,
    theorem_spherical_2_design,
    veronese_full_rank_proved_general,
)


def test_spherical_2_design_proved():
    assert theorem_spherical_2_design()["proved"]


def test_lambda_star_sufficiency_proved():
    t = theorem_lambda_star_sufficiency()
    assert t["proved"]
    for p in ("5", "7", "11"):
        assert t["by_p_sample"][p]["positive"]


def test_four_design_shortcut_dead():
    t = theorem_four_design_shortcut_dead()
    assert t["proved"]
    assert t["n_distinct_ge_2"]
    assert t["mult_sum_eq_D"]
    assert t["n_distinct"] == 3
    assert t["mult_sum"] == 65


def test_dual_four_thirds_budget():
    assert theorem_dual_four_thirds_budget()["proved"]


def test_lambda_star_floor_census():
    c = certify_lambda_star_floor_small_p()
    assert c["certified"]
    assert not c["proved_general"]
    assert c["by_p"]["5"]["sharp"]
    assert c["by_p"]["5"]["ker_dim_eq_n"]
    assert Fraction(c["by_p"]["5"]["min_rayleigh"]) == Fraction(80, 13)


def test_cost_D_ratio_census():
    c = certify_cost_D_sum_ne_ratio()
    assert c["max_ratio_observed"] < c["four_thirds"]
    for row in c["by_p_session"].values():
        assert row["lt_4_3"]


def test_predicates_remain_false():
    assert free_e_bound_proved_general() is False
    assert veronese_full_rank_proved_general() is False
    assert lambda_star_floor_proved_general() is False
    assert residual_i_closed_via_213() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status_open():
    h = hinge_status_213()
    assert h["spherical_2_design_proved"]
    assert h["four_design_shortcut_dead"]
    assert h["lambda_star_floor_general"] is False
    assert h["residual_i_closed"] is False
    assert h["gsum_disj_lb_proved_general"] is False
