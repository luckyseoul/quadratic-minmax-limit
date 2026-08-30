"""Tests for Prop 15.304 — λ_j as torus-Gauss square with den |U|."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15301 import dim_Vplus
from e1_gmin_m4_prop15304 import (
    lambda_from_hatN,
    lambda_from_square,
    live_lambda,
    main,
    prove_hatN_identity,
    prove_open,
    prove_square_identity,
)


def test_hatN_matches_live_lambda():
    A = prove_hatN_identity()
    assert A["proved"] is True
    assert A["drop_delta_hit"] >= 1
    for p, j in ((5, 2), (7, 0), (7, 2)):
        live = live_lambda(p, j)
        pred = lambda_from_hatN(p, j, use_delta=True)
        assert abs(pred - live) < 1e-8


def test_drop_delta_breaks():
    live = live_lambda(5, 2)
    bad = lambda_from_hatN(5, 2, use_delta=False)
    assert abs(bad - live) > 1.0


def test_square_named_den_U():
    B = prove_square_identity()
    assert B["proved"] is True
    assert B["den1_hit"] >= 1
    assert B["dimVplus_hit"] >= 1
    assert B["named_den"].startswith("p+1")
    for p, j in ((5, 2), (5, 0), (7, 0), (7, 4)):
        live = float(live_lambda(p, j).real)
        pred = lambda_from_square(p, j, den=p + 1)
        assert abs(pred - live) < 1e-8


def test_wrong_den_breaks():
    live5 = float(live_lambda(5, 2).real)
    assert abs(lambda_from_square(5, 2, den=1) - live5) > 0.5
    live7 = float(live_lambda(7, 0).real)
    assert abs(lambda_from_square(7, 0, den=dim_Vplus(7)) - live7) > 0.5
    assert dim_Vplus(7) == 25
    assert dim_Vplus(7) != 7 + 1


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["live_square_beats_threshold"] is True
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert C["S_k_named_in_p"] is False
    assert C["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["hatN_contraction"] is True
    assert out["proved"]["square_named_den_U"] is True
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
