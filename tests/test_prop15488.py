"""Tests for Prop 15.488 — ridge CSP names |H+| at p=3,5."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15488 import (
    F_named,
    csp_count,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_p3_named_csp():
    A = prove_A()
    assert A["proved"] is True
    assert A["N"] == 6
    assert len(F_named(3)) == 4
    assert csp_count(3, flip_T1=True) != 6
    assert csp_count(3, F_named(3, drop_level=True)) != 6


def test_p5_named_csp():
    B = prove_B()
    assert B["proved"] is True
    assert B["N"] == 130
    assert len(F_named(5)) == 31
    assert csp_count(5, flip_T1=True) != 130
    assert csp_count(5, F_named(5, drop_level=True)) != 130
    assert csp_count(5, F_named(5, drop_affine=True)) != 130


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["N"] == 6
    assert out["B"]["N"] == 130
