"""Tests for Prop 15.448 — every Type+ L_1D slot named in p."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15298 import mu_minus, mu_plus
from e1_gmin_m4_prop15414 import Lpar_1d_named
from e1_gmin_m4_prop15415 import Lspl_1d_named
from e1_gmin_m4_prop15448 import (
    L1d_from_weights,
    L1d_named,
    Lns_mix_named,
    Lpm1_1d_inverted,
    Lpm1_1d_named,
    Lsq_11_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    slot_of,
    weights_named,
)


def test_chi_Fp_trivial():
    A = prove_A()
    assert A["proved"] is True
    F = _kernel_field(5)
    assert int(F["chi"][3]) == 1
    assert pow(3, 2, 5) == 4  # χ_p(3)=−1


def test_weights_and_fail_n10():
    B = prove_B()
    assert B["proved"] is True
    wt = weights_named(5, "ns_mm")
    assert wt["n10"] == 0
    assert wt["n01"] == Fraction(1, 3)
    bad = dict(wt)
    bad["n10"] = bad["n01"]
    assert L1d_from_weights(5, **bad) != L1d_from_weights(5, **wt)


def test_named_forms_and_lifts():
    assert Lpm1_1d_named(5) == Fraction(125, 3)
    assert Lpm1_1d_inverted(5) == Fraction(125, 4)
    assert Lpm1_1d_inverted(5) != Lpm1_1d_named(5)
    assert L1d_named(5, "Sub") == Lpar_1d_named(5) == Fraction(100, 3)
    assert Lsq_11_named(5) == 25
    assert Lsq_11_named(7) == Lspl_1d_named(7) == Fraction(539, 5)
    assert Lsq_11_named(13) == Lspl_1d_named(13)
    assert Lns_mix_named(5) == mu_plus(5) * mu_minus(5) == Fraction(25, 2)
    C = prove_C()
    assert C["proved"] is True
    assert "sq_11" not in C["live"]["5"]
    assert C["live"]["5"]["pm1"] == ["125/3"]
    assert C["live"]["7"]["sq_11"] == ["539/5"]


def test_slot_split_is_chi_not_tor():
    F = _kernel_field(7)
    # extra-++ 21,28 are sq_11 (Lspl), not a single Tor value
    assert slot_of(F, 7, 21) == "sq_11"
    assert slot_of(F, 7, 2) == "Sub"


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.448"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["chi_Fp_trivial"] is True
    assert out["proved"]["weights_named"] is True
    assert out["proved"]["L1d_slots_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["wrong_pm1_p5"] == "125/4"
