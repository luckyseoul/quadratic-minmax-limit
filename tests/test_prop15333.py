"""Tests for Prop 15.333 — split-c dead at p=11; no linlin for c(p)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15333 import (
    LIVE_C,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    split_c_guess,
    three_p_minus_2,
    ystar_stab,
)


def test_three_p_minus_2_is_p7_only():
    A = prove_A()
    assert A["proved"] is True
    assert three_p_minus_2(7) == 19
    assert three_p_minus_2(5) == 13
    assert three_p_minus_2(5) != LIVE_C[5]
    assert split_c_guess(3) == 0
    assert split_c_guess(5) == 1
    assert split_c_guess(7) == 19


def test_ystar_stab_kills_split_at_p11():
    B = prove_B()
    assert B["proved"] is True
    assert ystar_stab(5)["stab"] == 6
    assert ystar_stab(7)["stab"] == 4
    assert ystar_stab(11)["stab"] == 4
    assert ystar_stab(11)["is_pplus1"] is False
    assert ystar_stab(11)["stab"] != 12
    assert B["split_dead_as_multiplicity"] is True


def test_named_orbits_are_lower_bound():
    C = prove_C()
    assert C["proved"] is True
    assert C["ystar_orbit_p5"] == 100
    assert C["ystar_orbit_p7"] == 588
    assert C["apnc_empty_p13"] is True


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Hplus_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["three_p_minus_2_p5_only"] is True
    assert out["proved"]["split_dead_p11"] is True
    assert out["proved"]["named_orbits_are_LB"] is True
    assert out["proved"]["Hplus_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
