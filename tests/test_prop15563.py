"""Tests for Prop 15.563 — full-Ω not Fejer; A_full open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15563 import main, prove_A, prove_B, prove_open


def test_p7_two_energy_types_not_equi():
    A = prove_A()
    assert A["proved"] is True
    r7 = A["by_p"]["7"]
    assert r7["n_types"] == 2
    assert r7["patterns"]["(1, 1, 1, 3)"] == 1764
    assert r7["patterns"]["(1, 1, 2, 2)"] == 2646
    assert r7["n_full"] == 4410
    assert A["by_p"]["5"]["n_types"] == 1


def test_full_not_fejer_not_two_level():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["7"]["min_fej_shape"] > 0.05
    assert B["by_p"]["7"]["n_two_level"] == 0
    assert B["by_p"]["5"]["min_fej_shape"] < 1e-8


def test_flags_untouched_A_full_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["A_full_named_in_p"] is False
    assert C["Q3_double_named_in_p"] is False
    assert C["Q3_n3_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.563"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15563.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
