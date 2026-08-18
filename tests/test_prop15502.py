"""Tests for Prop 15.502 — pairing corollary; S0 vs forms; Aut orbits."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15497 import pair_named
from e1_gmin_m4_prop15502 import (
    D_S0_from_forms,
    main,
    pair_drop_shift,
    pair_from_Q_forms,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_pairing_is_formal_corollary():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"][5]["n_match"] == A["by_p"][5]["n_even"]
    assert A["by_p"][13]["n_match"] == A["by_p"][13]["n_even"]
    assert pair_from_Q_forms(5, 2, 13) == pair_named(5, 2, 13)
    assert pair_drop_shift(5, 2, 13) != pair_named(5, 2, 13)


def test_S0_forms_integer_only_p5():
    B = prove_B()
    assert B["proved"] is True
    assert D_S0_from_forms(5) == 13
    assert D_S0_from_forms(13).denominator != 1
    assert D_S0_from_forms(13) != Fraction(13)


def test_classes_not_single_orbits_past_p5():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"][5]["n_bad"] == 0
    assert C["by_p"][13]["n_bad"] > 0
    assert C["by_p"][13]["sub_n"] != C["by_p"][13]["sub_orb"]
    assert C["by_p"][13]["sub_n"] == 10
    assert 13 - 3 > 2


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
