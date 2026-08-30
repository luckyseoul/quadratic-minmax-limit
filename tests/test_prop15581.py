"""Tests for Prop 15.581 — p=5 named Q++/q² meets 15.507 floor."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15507 import Qpp_ub_J2
from e1_gmin_m4_prop15581 import (
    Qpp_p5_drop_k3,
    Qpp_p5_named,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_named_qpp_is_48_13_and_drop_k3_fails():
    A = prove_A()
    assert A["proved"] is True
    assert Qpp_p5_named() == Fraction(48, 13)
    assert Qpp_p5_drop_k3() == Fraction(16, 3)
    assert Qpp_p5_drop_k3() > Qpp_ub_J2(5)
    assert Qpp_p5_named() != Qpp_ub_J2(5)
    assert Qpp_p5_named() != Fraction(26, 7)


def test_floor_and_flags():
    B = prove_B()
    assert B["proved"] is True
    assert Qpp_p5_named() < Qpp_ub_J2(5)
    assert Qpp_ub_J2(5) == Fraction(26, 7)
    assert Fraction(26, 7) - Fraction(48, 13) == Fraction(2, 91)
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.581"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15581.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "15.568" in text
