"""Tests for Prop 15.558 — J_all / J_T named; Q3_T still live."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15558 import (
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    q302_drop_plus1,
    q302_named,
)


def test_J_all_eight_term():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        row = A["by_p"][str(p)]
        assert row["ok"] is True
        assert row["drop111"] > 1.0
        assert row["drop_eighth"] > 1.0


def test_J_T_via_S_Omega():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7, 11, 13):
        row = B["by_p"][str(p)]
        assert row["ok"] is True
        assert row["J1_eq_Jall_err"] > 1.0


def test_G3_identity_and_flags():
    C = prove_C()
    assert C["proved"] is True
    for p in (5, 7):
        row = C["by_p"][str(p)]
        assert row["ok"] is True
        assert row["err_G02_is_free"] > 1.0
        assert row["err_one_Q3"] > 1.0
    D = prove_open()
    assert D["proved"] is False
    assert D["Q3_named_in_p"] is False
    assert D["NUM_SUM_16pA_general"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.558"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15558.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text


def test_Q3_02_certified_formula():
    """Field F=-2(3p²+2) at p≤23 ⇒ Q3_02=-4N(2p²+1)/p. Not a 2-point fit."""
    assert abs(q302_named(5, 130) + 5304.0) < 1e-9
    assert abs(q302_named(7, 5726) + 323928.0) < 1e-9
    assert abs(q302_named(5, 130) - q302_drop_plus1(5, 130)) > 1.0
    assert abs(q302_named(7, 5726) - q302_drop_plus1(7, 5726)) > 1.0
