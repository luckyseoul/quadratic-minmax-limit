"""Tests for Prop 15.314 — named T-types, U-split, F_1 floor rewrite."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15314 import (
    F1_named,
    main,
    named_sizes,
    prove_F1_and_floor,
    prove_open,
    prove_type_sizes,
    prove_U_split,
)


def test_type_sizes_named_and_not_constant_6():
    A = prove_type_sizes()
    assert A["proved"] is True
    assert named_sizes(5)["++"] == 6
    assert named_sizes(7)["++"] == 6
    assert named_sizes(11)["++"] == 12
    assert named_sizes(13)["++"] == 22
    assert named_sizes(5)["N*"] == 4
    assert named_sizes(7)["N*"] == 12
    assert named_sizes(7)["--N1"] == 4
    assert named_sizes(5)["--N1"] == 0
    assert A["six_constant_false"] is True


def test_U_split_mixed_vs_squares():
    B = prove_U_split()
    assert B["proved"] is True
    assert B["bothplus_on_p1mod4"] is False
    assert B["swap_odd_is_wrong"] is True
    assert B["by_p"]["5"]["bothplus"] == 0
    assert B["by_p"]["7"]["split_ok"] is True


def test_F1_recovers_live_floor():
    C = prove_F1_and_floor()
    assert C["proved"] is True
    assert C["drop_special_j_fails"] is True
    assert C["by_p"]["5"]["Qpp"] == "48/13"
    assert C["by_p"]["5"]["Qn"] == "32/13"
    assert abs(C["by_p"]["5"]["named_min"] - float(Fraction(80, 13))) < 1e-9
    assert C["by_p"]["7"]["Qpp"] == "1544/409"
    assert abs(C["by_p"]["7"]["named_min"] - float(Fraction(3072, 409))) < 1e-9


def test_F1_special_j_not_generic():
    # p=7 j=4 = (p+1)/2 is not 16-2 Q++
    p = 7
    Qpp = float(Fraction(1544, 409))
    Qmm = float(Fraction(1376, 409))
    special = F1_named(p, 4, Qpp, Qmm)
    generic = 16.0 - 2.0 * Qpp
    assert abs(special - generic) > 0.1


def test_floor_still_open_general():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["type_sizes"] is True
    assert out["proved"]["U_split"] is True
    assert out["proved"]["F1_rewrite"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
