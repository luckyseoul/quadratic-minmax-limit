"""Prop 15.199 — frame K4 identity; Wick Rayleigh; residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15198 import wick_hi
from e1_gmin_m4_prop15199 import (
    A2_budget_for_wick,
    K4_bound_proved_general,
    certify_frame_constants,
    hinge_status_199,
    lambda_max_need_for_wick,
    main,
    residual_i_closed_via_frame,
    sigma2_Xc,
    theorem_frame_identity,
    theorem_lambda_max_path_dead,
    theorem_wick_target_rayleigh,
)


def test_frame_and_wick_algebra():
    assert theorem_frame_identity()["proved"] is True
    th = theorem_wick_target_rayleigh()
    assert th["proved"] is True
    for p in (5, 7, 11, 13):
        n = n_of(p)
        assert sigma2_Xc(p) == n * (n - 2)
        assert A2_budget_for_wick(p) == 8 * n * n + 48 * n
        assert lambda_max_need_for_wick(p) == Fraction(8 * (n + 6), n - 2)
        # consistency with Wick: 4n² + need_A2 = Wick_hi
        assert 4 * n * n + A2_budget_for_wick(p) == wick_hi(p)


def test_lambda_max_path_documented_dead():
    th = theorem_lambda_max_path_dead()
    assert th["proved_dead"] is True
    assert lambda_max_need_for_wick(3) == 16
    assert lambda_max_need_for_wick(5) == Fraction(32, 3)


def test_census_frame_p3_p5_drive_shipped():
    """Drive certify_frame_constants on real Max+ path."""
    c3 = certify_frame_constants(3)
    assert c3["norm_match"] is True
    assert c3["ySy_match"] is True
    assert c3["K4_le_Wick"] is True
    assert abs(c3["K4"] - 1680) < 1e-6
    assert abs(c3["E_R"] - 16.0) < 1e-6  # equality case

    c5 = certify_frame_constants(5)
    assert c5["norm_match"] is True
    assert c5["ySy_match"] is True
    assert c5["K4_le_Wick"] is True
    # E[R] from K4=120400/13
    n = 26
    ER = (Fraction(120400, 13) - 4 * n * n) / (n * (n - 2))
    assert abs(c5["E_R"] - float(ER)) < 1e-6
    assert c5["E_R"] < float(lambda_max_need_for_wick(5))


def test_predicates_open():
    assert K4_bound_proved_general() is False
    assert residual_i_closed_via_frame() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_199()
    assert h["frame_identity_K4"] is True
    assert h["lambda_max_path_dead_p5"] is True
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["frame_identity"] is True
    assert out["proved"]["K4_bound_general"] is False
    assert out["proved"]["L_closed"] is False
