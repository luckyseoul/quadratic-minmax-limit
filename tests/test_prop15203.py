"""Prop 15.203: dual construction Max+-free; residual_i stays open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15203 import (
    dual_construction_cost,
    free_e_bound_proved_general,
    residual_i_closed_via_203,
    theorem_alpha_n_minus_1_beats_need,
)


def test_alpha_n_minus_1_beats_need():
    r = theorem_alpha_n_minus_1_beats_need()
    assert r["proved"] is True


def test_dual_construction_p5_blocks_if_ker_sc():
    d = dual_construction_cost(5)
    assert d["ok"] is True
    assert d["rnorm"] < 1e-10
    assert d["minW"] >= -1e-10
    assert d["degstd"] < 1e-8
    assert d["blocks_if_ker_eq_sc"] is True
    assert d["cost"] < d["need"]


def test_dual_construction_p3_does_not_block():
    d = dual_construction_cost(3)
    assert d["ok"] is True
    # dual-eq feasible at p=3; construction must not falsely claim block
    assert d["blocks_if_ker_eq_sc"] is False


def test_predicates_remain_open():
    assert free_e_bound_proved_general() is False
    assert residual_i_closed_via_203() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
