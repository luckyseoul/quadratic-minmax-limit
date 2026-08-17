"""Tests for Prop 15.331 — NL Aut_∞ stabs | p+1; 409|Φ_24(7) p=7-only."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15308 import G_order, G_order_no_frob
from e1_gmin_m4_prop15331 import (
    HPLUS,
    Phi4,
    Phi12,
    Phi24,
    c_from_stabs,
    main,
    nl_stabs,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_orbit_stab_recovers_Hplus():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["n_1d"] == 30
    assert A["by_p"]["7"]["n_1d"] == n_1d(7) == 140
    assert A["by_p"]["5"]["Nplus"] == HPLUS[5]
    assert A["by_p"]["7"]["recon"] == HPLUS[7] - n_1d(7)
    assert G_order_no_frob(5) != G_order(5)


def test_nl_stabs_divide_pplus1():
    B = prove_B()
    assert B["proved"] is True
    assert nl_stabs(3) == []
    assert nl_stabs(5) == [6]
    assert sorted(nl_stabs(7)) == [2, 2, 2, 2, 4, 8]
    assert c_from_stabs(5, nl_stabs(5)) == 1
    assert c_from_stabs(7, nl_stabs(7)) == 19
    assert 3 not in nl_stabs(7)
    assert 28 not in nl_stabs(7)
    assert 20 not in nl_stabs(5)
    assert 8 % 3 != 0


def test_409_divides_Phi24_7_not_equal():
    C = prove_C()
    assert C["proved"] is True
    assert Phi24(7) == 409 * 73 * 193
    assert Phi24(7) != 409
    assert Phi12(5) == 601
    assert Phi12(5) % 13 != 0
    assert Phi4(5) == 26
    assert Phi4(5) % 13 == 0
    assert (3 * 7 - 2) == 19
    assert (3 * 5 - 2) != 1


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Hplus_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["orbit_stab_sum"] is True
    assert out["proved"]["nl_stabs_div_pplus1"] is True
    assert out["proved"]["cyclo_409_is_p7_only"] is True
    assert out["proved"]["Hplus_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
