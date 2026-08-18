"""Tests for Prop 15.514 — S2 quadratic in M; Q_τ still unnamed."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15514 import (
    main,
    prove_A,
    prove_B,
    prove_open,
    s2_from_jcnt,
    s2_named,
    s2_vertex_named,
    s2_vertex_wrong,
)


def test_s2_coefs_match_jcnt():
    A = prove_A()
    assert A["proved"] is True
    assert s2_named(5) == (Fraction(10), Fraction(-96), Fraction(256))
    assert s2_named(5) == s2_from_jcnt(5)
    assert s2_named(5)[0] != 0
    assert s2_vertex_named(5) == Fraction(24, 5)
    assert s2_vertex_named(5) != s2_vertex_wrong(5)


def test_pos4_majorant_too_loose():
    B = prove_B()
    assert B["proved"] is True
    assert Fraction(B["s2_pos4_ub"]) > Fraction(B["s2_need"])


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
