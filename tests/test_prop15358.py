"""Tests for Prop 15.358 — catalog miss; quadratic dead; 2-type above ub."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15350 import Qpp_from_L
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15358 import (
    H_quadratic,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_catalog_and_wrong_guesses():
    A = prove_A()
    assert A["proved"] is True
    assert A["n_D_hits"] == 0
    assert A["n_QNL_hits"] == 0
    assert A["wrong"]["3p-2_at_p5"] == 13
    assert A["wrong"]["dimVp_at_p7"] == 25
    assert A["wrong"]["3p-2_at_p5"] != 1


def test_quadratic_not_2p_at_p11():
    B = prove_B()
    assert B["proved"] is True
    assert H_quadratic(3) == 6
    assert H_quadratic(5) == 130
    assert H_quadratic(7) == 5726
    assert H_quadratic(11) % 22 != 0


def test_twotype_above_ub():
    C = prove_C()
    assert C["proved"] is True
    assert Qpp_from_L(13) > Qpp_floor_ub(13)
    assert Q_1d_pp_named(13) > Qpp_floor_ub(13)
    assert Qpp_from_L(13) != Qpp_floor_ub(13)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert D["Q_NL_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["catalog_miss"] is True
    assert out["proved"]["quadratic_not_2p_multiple"] is True
    assert out["proved"]["twotype_above_ub_p1"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
