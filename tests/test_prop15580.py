"""Tests for Prop 15.580 — Max− Ω_− Galois + F; DEAD as Type I kill."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15564 import F_named
from e1_gmin_m4_prop15580 import main, prove_A, prove_B, prove_open


def test_maxminus_galois_support():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["z0_err"] < 1e-6
    assert A["rows"]["5"]["nsupp_Om_plus_max"] == 0
    assert A["rows"]["5"]["kcnt"] == {"1": 30, "3": 100}
    assert A["rows"]["7"]["kcnt"] == {"1": 140, "3": 1176, "4": 4410}
    assert A["rows"]["5"]["n_partial_line"] == 0
    # fail-when-wrong: nsupp=18 legal at p=5
    assert 18 % 4 != 0


def test_F_same_on_Omega_minus():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["rows"]["5"]["F"] - F_named(5)) < 1e-4
    assert abs(B["rows"]["7"]["F"] - F_named(7)) < 1e-4
    assert abs(B["rows"]["11"]["F"] - F_named(11)) < 1e-4
    assert F_named(5) == -2 * (3 * 25 + 2)
    # fail: drop +2
    assert F_named(5) != -2 * (3 * 25)


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["path_dead"] is True
    assert C["type_I_multilevel_bad_case_ND_closed"] is False
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["maxminus_omega_minus_galois"] is True
    assert out["proved"]["F_on_Omega_minus"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
