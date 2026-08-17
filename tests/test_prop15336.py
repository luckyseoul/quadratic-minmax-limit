"""Tests for Prop 15.336 — isolated p=7 NL orbit; still 15.x."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15336 import (
    ISOL_Q,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_nc_p5_connected_p7_isolated():
    A = prove_A()
    assert A["proved"] is True
    assert A["by"]["5"]["n_isolated"] == 0
    assert A["by"]["5"]["deg_min"] >= 1
    assert A["by"]["7_isol"]["n"] == 1176
    assert A["by"]["7_isol"]["all_zero"] is True


def test_isolated_Q_is_8_3_not_ensemble():
    B = prove_B()
    assert B["proved"] is True
    assert B["Q"]["++"] == "8/3"
    assert B["Q"]["N*"] == "4"
    assert B["Q"]["--N1"] == "4"
    assert ISOL_Q["++"] != LIVE_QPP[7]
    assert ISOL_Q["++"] != Fraction(40, 9)


def test_one_full_line_involution():
    C = prove_C()
    assert C["proved"] is True
    assert C["full_min"] == 1
    assert C["full_max"] == 1
    assert C["stab"] == 2
    assert C["invol_minus"] is True


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["isolated_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["nc_isolates_1176_at_p7"] is True
    assert out["proved"]["isolated_Q_wick_off_pp"] is True
    assert out["proved"]["one_full_line_and_involution"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
