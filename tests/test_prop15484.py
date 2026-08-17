"""Tests for Prop 15.484 — TR L++ is not a formula in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15484 import (
    L_tr_drop12,
    L_tr_interp,
    agl_tr_maxplus_L,
    live_cache_L_tr,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_interp_hits_live():
    A = prove_A()
    assert A["proved"] is True
    assert L_tr_interp(5) == 24
    assert L_tr_interp(7) == Fraction(1337, 12)
    assert live_cache_L_tr(5) == 24
    assert live_cache_L_tr(7) == Fraction(1337, 12)
    assert L_tr_drop12(7) != Fraction(1337, 12)


def test_interp_dies_p11():
    B = prove_B()
    assert B["proved"] is True
    pred = L_tr_interp(11)
    assert pred == Fraction(9089, 12)
    Ls = [Fraction(s) for s in B["values"]]
    assert Ls
    assert pred not in Ls


def test_not_function_of_p():
    C = prove_C()
    assert C["proved"] is True
    assert C["nuniq_p11"] >= 2
    Ls = agl_tr_maxplus_L(11, max_tries=200, seed=1, want=4)
    assert len(set(Ls)) >= 2


def test_raw_T():
    D = prove_D()
    assert D["proved"] is True
    assert D["raw5"] == "18"


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.484"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["interp_dies_p11"] is True
    assert out["proved"]["TR_L_not_fn_of_p"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
