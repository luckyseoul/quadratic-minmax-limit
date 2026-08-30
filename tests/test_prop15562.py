"""Tests for Prop 15.562 — k=3 Fejer names A_k3 and Q3_n3 piece."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15562 import (
    A_k3_n3_named,
    A_k3_n3_wrong_pminus2,
    S_drop24,
    S_named,
    S_wrong_pm1,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    q3_n3_k3_drop_third,
    q3_n3_k3_named,
)


def test_S_is_minus_p2_minus_1_over_24():
    A = prove_A()
    assert A["proved"] is True
    assert abs(S_named(5) + 1.0) < 1e-12
    assert abs(S_named(7) + 2.0) < 1e-12
    assert abs(S_named(5) - S_wrong_pm1(5)) > 0.5
    assert abs(S_named(5) - S_drop24(5)) > 0.5
    for p in (5, 7, 11, 13, 17, 19, 23):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["pm1_fails"] is True
        assert row["drop24_fails"] is True


def test_A_k3_dbl_zero_and_A_k3_n3_named():
    B = prove_B()
    C = prove_C()
    assert B["proved"] is True
    assert C["proved"] is True
    assert abs(A_k3_n3_named(5) + 250.0) < 1e-9
    assert abs(A_k3_n3_named(7) + 686.0 / 3.0) < 1e-9
    assert abs(A_k3_n3_named(7) - A_k3_n3_wrong_pminus2(7)) > 1.0
    for p in (5, 7):
        assert B["by_p"][str(p)]["match"] is True
        assert C["by_p"][str(p)]["match"] is True
        assert C["by_p"][str(p)]["pminus2_fails"] is True


def test_k3_Q3_n3_piece_and_flags():
    D = prove_D()
    assert D["proved"] is True
    assert abs(q3_n3_k3_named(5) + 25000.0) < 1e-6
    assert abs(q3_n3_k3_named(7) + 268912.0) < 1e-4
    assert abs(q3_n3_k3_named(5) - q3_n3_k3_drop_third(5)) > 1.0
    # p=5: piece is all of Q3_n3; p=7: piece is not
    assert D["by_p"]["5"]["p5_Q3n3_or_p7_not_full"] is True
    assert D["by_p"]["7"]["p5_Q3n3_or_p7_not_full"] is True
    E = prove_open()
    assert E["proved"] is False
    assert E["Q3_double_named_in_p"] is False
    assert E["Q3_n3_named_in_p"] is False
    assert E["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.562"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is True
    assert out["E"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15562.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
