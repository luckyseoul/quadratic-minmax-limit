"""Tests for Prop 15.575 — μ_1d is the Type+ Johnson 4-point."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15575 import (
    main,
    mu_1d_16p4_over_3,
    mu_1d_const_eps,
    mu_1d_drop4,
    mu_1d_named,
    prove_A,
    prove_open,
)


def test_mu_1d_johnson_hits_and_const_fails():
    A = prove_A()
    assert A["proved"] is True
    assert mu_1d_named(5) == Fraction(10000, 3)
    assert mu_1d_named(7) == Fraction(124852, 5)
    assert abs(float(mu_1d_named(5)) - A["by_p"]["5"]["live"]) < 1e-6
    assert abs(float(mu_1d_named(7)) - A["by_p"]["7"]["live"]) < 1e-6
    assert abs(float(mu_1d_named(5)) - float(mu_1d_const_eps(5))) > 1.0
    assert abs(float(mu_1d_named(7)) - float(mu_1d_16p4_over_3(7))) > 1.0
    assert abs(float(mu_1d_named(7)) - float(mu_1d_drop4(7))) > 1.0
    assert A["by_p"]["5"]["const_fails"] is True
    assert A["by_p"]["7"]["sixteen_fails_or_p5"] is True
    assert A["by_p"]["11"]["match"] is True
    assert A["by_p"]["13"]["match"] is True


def test_flags_and_open():
    B = prove_open()
    assert B["proved"] is False
    assert B["Q_sub_closed_in_p"] is False
    assert B["mu_full_named_in_p"] is False
    assert B["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.575"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15575.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "Aut_e DEAD" in text
    assert "16p⁴/3" in text or "16p4/3" in text
