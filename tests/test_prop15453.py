"""Tests for Prop 15.453 — ns-net N_r = |H+| at p=7."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15453 import (
    N4_P7,
    NS_NET7,
    main,
    mobius_pair,
    n4_pred,
    ns_net_slopes,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_n3_p5_is_Hplus():
    A = prove_A()
    assert A["proved"] is True
    assert A["n2"] == 2040
    assert set(A["n3"].values()) == {130}
    assert 30 not in A["n3"].values()
    assert A["n_1d"] == n_1d(5) == 30


def test_n4_type_split():
    B = prove_B()
    assert B["proved"] is True
    assert B["n5726"] == 9
    assert B["n4844"] == 6
    assert B["diff"] == 882
    assert n4_pred(1, 2) == 5726
    assert n4_pred(1, 3) == 4844
    assert N4_P7[(1, 3)] == 4844


def test_nr7_is_Hplus():
    C = prove_C()
    assert C["proved"] is True
    sl = ns_net_slopes(7)
    assert set(int(s) for s in sl) == set(NS_NET7)
    assert C["all_5726"] is True
    assert C["Hplus"] == HPLUS[7] == 5726
    # fail-when-wrong: {0,∞,1,3} is the other type
    assert n4_pred(1, 3) == 4844
    s = mobius_pair(1, 2, 5, 7)
    t = mobius_pair(1, 2, 6, 7)
    assert N4_P7[tuple(sorted((s, t)))] == 5726


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.453"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["n3_p5_is_Hplus"] is True
    assert out["proved"]["n4_p7_type_split"] is True
    assert out["proved"]["nr7_is_Hplus"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
