"""Tests for Prop 15.325 — A_U spectrum on Paley θ+; two-point too weak."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15325 import (
    AU_eigs_named,
    AU_eigs_squares,
    main,
    prove_common_neighbors,
    prove_open,
    prove_spectrum,
    prove_twopoint_weak,
    theta_plus_norms,
    two_point_var_ub_float,
)


def test_spectrum_is_Kl_on_Omega_sign_not_squares_at_p5():
    A = prove_spectrum()
    assert A["proved"] is True
    assert A["squares_all_p_false"] is True
    e5 = sorted(AU_eigs_named(5))
    assert abs(e5[0] - (1 - 5**0.5)) < 1e-8
    assert abs(e5[1] - (1 + 5**0.5)) < 1e-8
    assert sorted(round(x, 8) for x in AU_eigs_named(5)) != sorted(
        round(x, 8) for x in AU_eigs_squares(5)
    )
    # p≡3: θ+ norms are squares, lists agree
    assert sorted(round(x, 8) for x in AU_eigs_named(7)) == sorted(
        round(x, 8) for x in AU_eigs_squares(7)
    )
    assert len(theta_plus_norms(5)) == 2
    assert len(theta_plus_norms(13)) == 6


def test_common_neighbors_N_constant():
    B = prove_common_neighbors()
    assert B["proved"] is True
    assert B["fiber_split_false"] is True
    for lo, hi in B["by_p"]["5"].values():
        assert lo == hi


def test_twopoint_is_180_and_misses_floor():
    C = prove_twopoint_weak()
    assert C["proved"] is True
    assert C["twopoint_closes_p5_false"] is True
    assert abs(two_point_var_ub_float(5) - 180.0) < 0.05
    assert two_point_var_ub_float(5) > 360 / 7


def test_majorant_and_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["A_square_named"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["AU_spectrum_named"] is True
    assert out["proved"]["nn_N_constant"] is True
    assert out["proved"]["twopoint_too_weak"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
