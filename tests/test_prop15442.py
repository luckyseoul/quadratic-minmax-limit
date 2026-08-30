"""Tests for Prop 15.442 — p=7 NL multiplicity is not a p≡3 law."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15442 import (
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    same_orbit,
    stab_of,
    ystar_odot_second_stab4,
)
from e1_gmin_m4_prop15138 import construct_ystar


def test_second_stab4_not_ystar_orbit():
    rec = ystar_odot_second_stab4(11)
    assert rec["found"] is True
    assert rec["stab"] == 4
    y0 = construct_ystar(11)["y"]
    assert stab_of(11, y0) == 4
    assert same_orbit(11, y0, rec["y"]) is False
    A = prove_A()
    assert A["proved"] is True
    assert A["same_orbit_as_ystar"] is False


def test_c_lo_and_p7_not_19():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["c_lo"] == 1
    assert B["by_p"]["7"]["c_lo"] == 14
    assert B["by_p"]["7"]["c_lo"] != 19
    assert B["by_p"]["7"]["n_stab2"] == 3
    assert B["by_p"]["11"]["c_lo"] == 18
    assert B["by_p"]["11"]["stabs"].count(4) >= 2


def test_p7_pattern_dead():
    C = prove_C()
    assert C["proved"] is True
    assert C["m4_p11"] >= 2
    assert C["p7_pattern"] == [2, 2, 2, 2, 4, 8]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.442"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["second_stab4_p11"] is True
    assert out["proved"]["c_lo_table"] is True
    assert out["proved"]["p7_pattern_not_p3_law"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
