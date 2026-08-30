"""Tests for Prop 15.296 — X_A is ψ-mode of u; naive partners dead."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import bool_coll_off
from e1_gmin_m4_prop15296 import (
    main,
    prove_Q_lt_Coll_p5,
    prove_XA_is_psi_mode,
    prove_not_two_level,
    prove_open,
    prove_sum_u_constant,
)


def test_XA_is_psi_mode():
    A = prove_XA_is_psi_mode([5])
    assert A["proved"] is True
    assert A["sample"]["5"]["max_diff"] < 1e-6
    assert A["invert_hit"] < A["invert_tot"]


def test_sum_u_constant():
    B = prove_sum_u_constant([5])
    assert B["proved"] is True
    assert B["sample"]["5"]["want"] == 25 * 24


def test_not_two_level():
    C = prove_not_two_level()
    assert C["proved"] is True
    assert C["sample"]["5"]["n_unique"] >= 5
    assert C["sample"]["5"]["two_level"] is False
    assert C["sample"]["7"]["two_level"] is False


def test_Q_lt_Coll_p5():
    D = prove_Q_lt_Coll_p5()
    assert D["proved"] is True
    assert D["min_Q_off_over_q2"] < D["coll_off_over_q2"]
    # 32/13 vs 3-10/25
    assert abs(D["coll_off_over_q2"] - (3 - 10 / 25)) < 1e-9
    assert bool_coll_off(5) == 3 * 25 * 25 - 10 * 25


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["XA_is_psi_mode"] is True
    assert out["proved"]["sum_u_constant"] is True
    assert out["proved"]["u_not_two_level"] is True
    assert out["proved"]["Q_lt_Coll_p5"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
