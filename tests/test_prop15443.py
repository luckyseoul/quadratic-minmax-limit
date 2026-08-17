"""Tests for Prop 15.443 — c=c_NC+4·1_{p≡3} at p=5,7."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15411 import LIVE_QNL
from e1_gmin_m4_prop15439 import Q_NL_from_c
from e1_gmin_m4_prop15443 import (
    YSTARS8,
    c_NC,
    c_from_split,
    extra_pp,
    main,
    nc_closure_nl_stabs,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    stab_of,
    ystar_odot_tc,
)


def test_stab8_and_c_NC():
    y8 = ystar_odot_tc(7, *YSTARS8)
    assert y8 is not None
    assert stab_of(7, y8) == 8
    assert tuple(sorted(nc_closure_nl_stabs(5))) == (6,)
    assert sorted(nc_closure_nl_stabs(7)) == [2, 2, 2, 4, 8]
    assert c_NC(5) == 1
    assert c_NC(7) == 15
    assert c_NC(7) != 19
    A = prove_A()
    assert A["proved"] is True
    assert A["drop8"] == 14


def test_split_and_fails():
    assert extra_pp(5) == 0
    assert extra_pp(7) == 4
    assert c_from_split(5) == 1
    assert c_from_split(7) == 19
    assert c_NC(7) + 0 != 19
    assert c_NC(5) + 4 == 5 != 1
    B = prove_B()
    assert B["proved"] is True
    assert B["drop4"] == 15
    assert B["add4_p5"] == 5


def test_Q_NL_from_this_c():
    assert Q_NL_from_c(5, c_from_split(5)) == LIVE_QNL[5]
    assert Q_NL_from_c(7, c_from_split(7)) == LIVE_QNL[7]
    assert Q_NL_from_c(7, c_NC(7)) != LIVE_QNL[7]
    C = prove_C()
    assert C["proved"] is True


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.443"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["c_NC_closure"] is True
    assert out["proved"]["c_split_live"] is True
    assert out["proved"]["Q_NL_from_split"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["c_NC7"] == 15
    assert out["algebra"]["B"]["c7"] == 19
