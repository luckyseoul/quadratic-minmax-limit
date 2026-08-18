"""Tests for Prop 15.513 — even pairing S1=-16; Q_τ still unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15513 import (
    main,
    n_J2_named,
    n_J2_wrong,
    n_Jm1_named,
    n_Jm3_named,
    prove_A,
    prove_B,
    prove_open,
    s1_coefs,
)


def test_jcnt_named():
    A = prove_A()
    assert A["proved"] is True
    assert n_J2_named(5) == 4
    assert n_Jm1_named(5) == 2
    assert n_Jm3_named(5) == 4
    assert n_J2_named(5) != n_J2_wrong(5)
    assert A["rows"]["13"]["n2"] == 60


def test_s1_is_minus_16():
    B = prove_B()
    assert B["proved"] is True
    c_M, c_0 = s1_coefs(5)
    assert c_M == 0
    assert c_0 == -16
    assert B["rows"]["17"]["c_0"] == "-16"


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["D_named_in_p"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
