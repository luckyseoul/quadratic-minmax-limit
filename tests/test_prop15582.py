"""Tests for Prop 15.582 — 1D vanishes off F_p; p=5 Q_N*=32/13."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15582 import (
    Qnstar_p5_named,
    Qnstar_p5_with_1d,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_1d_vanishes_off_Fp():
    A = prove_A()
    assert A["proved"] is True
    for p in ("5", "7"):
        assert abs(A["by_p"][p]["mu_1d_N*"]) < 1e-6
        assert abs(A["by_p"][p]["mu_1d_N1"]) < 1e-6
        assert A["by_p"][p]["ne_sub"] is True


def test_qnstar_p5_and_open():
    B = prove_B()
    assert B["proved"] is True
    assert Qnstar_p5_named() == Fraction(32, 13)
    assert Qnstar_p5_with_1d() == Fraction(48, 13)
    assert Qnstar_p5_named() != Qnstar_p5_with_1d()
    C = prove_open()
    assert C["proved"] is False
    assert C["k3_2line_eq_1line"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["sixteen_p_minus_4_not_imported"] is True
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.582"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15582.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "16(p−4)" in text or "16(p-4)" in text
