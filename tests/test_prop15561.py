"""Tests for Prop 15.561 — n_k3 and Q3_generic named."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15561 import (
    A1d_named,
    main,
    n_k3_C_m2,
    n_k3_drop_pm1,
    n_k3_named,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    q3_generic_drop_nk3,
    q3_generic_named,
)


def test_n_k3_named():
    A = prove_A()
    assert A["proved"] is True
    assert n_k3_named(5) == 100
    assert n_k3_named(7) == 1176
    assert n_k3_C_m2(5) != 100
    assert n_k3_drop_pm1(7) != 1176
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["drop_pm1_fails"] is True
        assert row["C_m2_fails"] is True


def test_A_k3_vanishes_on_1line():
    B = prove_B()
    assert B["proved"] is True
    assert abs(A1d_named(7) + 274.4) < 1e-9
    for p in (5, 7):
        row = B["by_p"][str(p)]
        assert row["match"] is True
        assert row["A1d_fails"] is True


def test_Q3_generic_named_and_flags():
    C = prove_C()
    assert C["proved"] is True
    assert abs(q3_generic_named(7, 5726) + 1248520.0) < 1e-6
    assert abs(q3_generic_drop_nk3(7, 5726) + 1248520.0) > 1.0
    assert C["by_p"]["5"]["live"] is None
    assert C["by_p"]["7"]["match"] is True
    assert C["by_p"]["7"]["drop_nk3_fails"] is True
    D = prove_open()
    assert D["proved"] is False
    assert D["Q3_double_named_in_p"] is False
    assert D["Q3_n3_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.561"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15561.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
