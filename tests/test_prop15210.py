"""Tests for Prop 15.210 — D(C) cost budget; hinge open."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15210 import (  # noqa: E402
    certify_cost_D_budget,
    cost_D_bound_proved_general,
    free_e_bound_proved_general,
    hinge_status_210,
    residual_i_closed_via_210,
    theorem_eight_over_p_beats_need,
    theorem_implication_costD_ker_sc,
)


def test_eight_over_p_proved():
    t = theorem_eight_over_p_beats_need()
    assert t["proved"]
    assert t["by_p_sample"]["5"]["beats"]


def test_implication_recorded():
    t = theorem_implication_costD_ker_sc()
    assert t["proved"]
    assert "ker_eq_sc_general" in t["hypotheses_open"]
    assert "cost_D_le_8_over_p_general" in t["hypotheses_open"]


def test_cost_D_census_small():
    """D(C) cost ≤ 8/p at p=5,7 (fast)."""
    c = certify_cost_D_budget([5, 7])
    assert c["certified"]
    assert c["proved_general"] is False
    assert c["max_cost_times_p"] < 8.0
    for p in ("5", "7"):
        row = c["by_p"][p]
        assert row["cost_le_8_over_p"]
        assert row["cost_lt_need"]
        assert row["blocks_if_ker_sc"]


def test_predicates_honest_false():
    assert free_e_bound_proved_general() is False
    assert cost_D_bound_proved_general() is False
    assert residual_i_closed_via_210() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_210()
    assert h["residual_i_closed"] is False
    assert h["cost_D_le_8_over_p_general"] is False
    assert h["ker_eq_sc_general"] is False
