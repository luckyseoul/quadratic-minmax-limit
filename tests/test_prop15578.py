"""Tests for Prop 15.578 — 2χ Fejér 4-point; p=7 buckets; mix open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15578 import (
    E4_2chi_named,
    main,
    mu_2chi_named,
    prove_A,
    prove_B,
    prove_open,
)


def test_twochi_fejer_and_fails():
    A = prove_A()
    assert A["proved"] is True
    assert E4_2chi_named(5) == Fraction(5, 1)
    assert E4_2chi_named(7) == Fraction(14, 3)
    assert mu_2chi_named(5) == Fraction(2000, 1)
    assert mu_2chi_named(7) != Fraction(126224, 15)
    assert abs(float(E4_2chi_named(7)) - 4.0) > 1e-9
    assert abs(float(E4_2chi_named(7)) - 5.0) > 1e-9
    assert A["by_p"]["7"]["r1_not_p4"] is True


def test_p7_buckets_and_mix_open():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["E4_A3"] - 140 / 3) < 1e-6
    assert abs(B["E4_A6"] - 35 / 3) < 1e-6
    assert abs(B["E4_B4w"] - 133 / 9) < 1e-6
    assert abs(B["E4_B4u"] - 133 / 9) > 0.1
    assert abs(B["E4_B4_dropAP"] - 133 / 9) > 0.1
    C = prove_open()
    assert C["proved"] is False
    assert C["mix_is_plaw"] is False
    assert C["mu_full_le_muk3_Nstar_p7"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.578"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15578.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "2800" in text
