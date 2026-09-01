"""Tests for Prop 15.534 — I is Gauss × cubic S, not a χ-type law."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15523 import G_chi
from e1_gmin_m4_prop15534 import (
    I_pred,
    I_pred_G_eq_p,
    I_pred_drop_chi,
    S_table,
    field_tables,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_G_S_formula_hits_live_I():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["formula_ok"] is True
    assert A["rows"]["7"]["formula_ok"] is True
    assert A["rows"]["5"]["S0"] == -1
    assert A["rows"]["5"]["S1"] == -1
    assert A["rows"]["7"]["S0"] == -1
    # Orin samples: I is not a single p-law
    assert abs(A["p5_I11"] + 30) < 1e-6
    assert abs(A["p7_I11"] - 98) < 1e-6
    assert 98 == 2 * 7 * 7
    assert -30 != 2 * 5 * 5


def test_drop_chi_and_G_eq_p_fail():
    A = prove_A()
    rec = A["p5_I_1p"]
    assert rec["pred"] == rec["live"]
    assert rec["drop"] != rec["pred"]
    assert rec["drop"] == -rec["pred"]
    F = field_tables(5)
    S = S_table(F)
    assert I_pred(F, S, 1, 1) == -30
    assert I_pred_G_eq_p(F, S, 1, 1) == 30
    assert I_pred_G_eq_p(F, S, 1, 1) != I_pred(F, S, 1, 1)
    assert G_chi(5) == -5
    assert I_pred_drop_chi(F, S, 1, 5) != I_pred(F, S, 1, 5)


def test_chi_types_do_not_pin_S():
    B = prove_B()
    assert B["proved"] is True
    assert set(B["witness"]["S"]) == {-6, 2}
    assert B["p7_n_split"] >= 1
    # fail-when-wrong: type constancy
    assert len(B["p5_split"]) >= 1


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["Delta_conn_p_rational"] is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["I_equals_G_chi_S"] is True
    assert out["proved"]["chi_types_do_not_pin_S"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
