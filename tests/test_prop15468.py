"""Tests for Prop 15.468 — Max+ χ-pair law; 3-point not constant."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15468 import (
    S_chi_drop_infty,
    S_chi_named,
    chi_plus_named,
    k_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_algebra_pair_law():
    A = prove_A()
    assert A["proved"] is True
    assert S_chi_named(5) == 30
    assert S_chi_named(7) == 84
    assert chi_plus_named(5) == 30
    assert chi_plus_named(7) == 126
    assert S_chi_drop_infty(5) == Fraction(65, 2)
    assert S_chi_drop_infty(5) != 30


def test_cache_pairs():
    B = prove_B()
    assert B["proved"] is True
    assert chi_plus_named(5) != Fraction(k_named(5) * (k_named(5) - 1), 4)


def test_triples_split():
    C = prove_C()
    assert C["proved"] is True
    assert C["nunique"] == 16
    assert C["nunique"] != 1


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.468"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["pair_law"] is True
    assert out["proved"]["triple_not_constant"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
