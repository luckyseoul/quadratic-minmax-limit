"""Tests for Prop 15.400 — floor λ≥6 is not pointwise on Max+."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15400 import (
    lambda_rows,
    main,
    one_d_weight_lb,
    pointwise_floor,
    prove_A,
    prove_B,
    prove_open,
    s1d_named,
)


def test_not_pointwise_p5_two_level():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["n_lambda_0"] == 100
    assert A["rows"]["5"]["n_lambda_s1d"] == n_1d(5) == 30
    assert A["rows"]["5"]["n_lt6"] == 100
    lam = lambda_rows(5, 2)
    assert (lam < 6 - 1e-9).sum() == 100
    # fail-when-wrong: not every Max+ is ≥6
    assert not (lam >= 6 - 1e-12).all()


def test_p7_NL_does_not_vanish():
    A = prove_A()
    assert A["rows"]["7"]["n_lt6"] > 0
    assert A["rows"]["7"]["n_lambda_0"] != HPLUS[7] - n_1d(7)
    assert A["rows"]["7"]["n0_eq_NL"] is False


def test_p5_mix_is_1d_indicator():
    B = prove_B()
    assert B["proved"] is True
    assert one_d_weight_lb(5) == Fraction(80, 13)
    assert s1d_named(5) == Fraction(80, 3)
    assert abs(B["NL_share"]) < 1e-8


def test_one_d_lb_misses_at_p7():
    C = prove_open()
    assert C["proved"] is False
    assert one_d_weight_lb(7) < 6
    assert pointwise_floor() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["not_pointwise"] is True
    assert out["proved"]["p5_mix_is_1d_indicator"] is True
    assert out["proved"]["pointwise_floor"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
