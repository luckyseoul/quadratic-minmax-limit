"""Tests for Prop 15.359 — D quadratic dead; no small 3-feature name."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15359 import (
    D_quad,
    count_p11_maxplus,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_D_quad_hits_live_and_not_p5_swap():
    A = prove_A()
    assert A["proved"] is True
    assert D_quad(3) == 1
    assert D_quad(5) == 13
    assert D_quad(7) == 409
    assert D_quad(11) == 2353
    assert D_quad(5) != 1


def test_p11_census_exceeds_quad_prediction():
    B = prove_B()
    assert B["proved"] is True
    rec = count_p11_maxplus()
    assert rec["n_maxplus"] == 85052
    assert rec["n_maxplus"] > 22 * D_quad(11)
    assert rec["D_lb"] > D_quad(11)
    assert rec["n_not"] == 58080


def test_no_small_cubic_or_3feature():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_cubic"] == 0
    assert C["n_prod3"] == 0
    assert C["n_lin3"] == 0
    assert C["fake_cubic_p3"] != 1


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["D_quad_formula"] is True
    assert out["proved"]["p11_census_kills_D_quad"] is True
    assert out["proved"]["no_small_cubic_or_3feature"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
