"""Tests for Prop 15.450 — Hoffman ∩ S0 is a 19-free c-window."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15409 import live_u, u_window
from e1_gmin_m4_prop15450 import (
    c_2p_plus_5,
    c_3p_minus_2,
    c_window_ends,
    in_c_window,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    u_from_c,
    u_from_c_wrong,
)


def test_hoffman_lattice():
    assert u_from_c(5, 1) == 2 == live_u(5)
    assert u_from_c(7, 19) == 57 == live_u(7)
    assert u_from_c_wrong(5, 1) == 1
    assert u_from_c_wrong(5, 1) != live_u(5)
    A = prove_A()
    assert A["proved"] is True


def test_p5_unique_and_interpolants_die():
    assert c_window_ends(5) == (1, 1, 1)
    assert in_c_window(5, 1) is True
    assert in_c_window(5, c_3p_minus_2(5)) is False
    assert in_c_window(5, c_2p_plus_5(5)) is False
    assert u_from_c(5, 13) == 26
    assert u_window(5) == (2, 2)
    B = prove_B()
    assert B["proved"] is True
    assert B["bad_3p2_u5"] == 26


def test_p13_window_misses_interpolants():
    c0, c1, n = c_window_ends(13)
    assert n == 92
    assert c0 == 551
    assert c1 == 642
    assert in_c_window(13, 37) is False
    assert in_c_window(13, 31) is False
    ulo, uhi = u_window(13)
    assert ulo <= 3306 <= uhi
    assert not (ulo <= 222 <= uhi)
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"]["13"]["n"] == 92


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["n_hits_13"] == 92
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.450"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["u_from_c"] is True
    assert out["proved"]["c_window_s0"] is True
    assert out["proved"]["interpolants_miss"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["c_window_5"] == [1, 1, 1]
    assert out["algebra"]["C"]["by_p"]["13"]["c_min"] == 551
