"""Tests for Prop 15.531 — axis n_R chart-dependent; lin-form = n_free."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15527 import split_free_orbits
from e1_gmin_m4_prop15531 import (
    RANDOM_LIN_AFF,
    main,
    n_X_wrong_98,
    prove_A,
    prove_B,
    prove_open,
)


def test_lin_form_is_all_free_orbits():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["n_R_lin"] == 4
    assert A["by_p"]["5"]["n_free"] == 4
    assert A["by_p"]["7"]["n_R_lin"] == 114
    assert A["by_p"]["7"]["n_X_lin"] == 0
    assert A["by_p"]["7"]["n_X_lin"] != n_X_wrong_98(7)


def test_axis_nR_not_translation_invariant_at_p7():
    A = prove_A()
    assert A["by_p"]["7"]["n_R_axis_minkey"] == 16
    assert A["by_p"]["7"]["n_R_axis_any"] == 20
    assert split_free_orbits(7)["n_R"] == 16


def test_random_ksets_not_all_linaff():
    B = prove_B()
    assert B["proved"] is True
    assert RANDOM_LIN_AFF[7]["hits"] != RANDOM_LIN_AFF[7]["n"]
    assert RANDOM_LIN_AFF[5]["hits"] != RANDOM_LIN_AFF[5]["n"]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_free_named_in_p"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
