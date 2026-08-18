"""Tests for Prop 15.564 — F=−2(3p²+2) names Q3_02 for every odd p."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15558 import LIVE_N, q302_drop_plus1, q302_named
from e1_gmin_m4_prop15564 import (
    F_named,
    F_wrong_I0,
    LIVE_Q302,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    q302_from_F,
)


def test_I_equals_two():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["I0_fails"] is True


def test_chi_minus_one_and_G2():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7, 11):
        row = B["by_p"][str(p)]
        assert row["match"] is True
        assert row["G2_neg_fails"] is True
        assert row["chi_m1"] == 1


def test_F_names_Q302_and_flags():
    C = prove_C()
    D = prove_D()
    assert C["proved"] is True
    assert D["proved"] is True
    for p in (5, 7, 11):
        assert C["by_p"][str(p)]["match"] is True
        assert C["by_p"][str(p)]["drop6_fails"] is True
        assert D["by_p"][str(p)]["match_F"] is True
        assert D["by_p"][str(p)]["I0_fails"] is True
        assert abs(F_named(p) - F_wrong_I0(p)) > 1.0
    for p in (5, 7):
        N = LIVE_N[p]
        assert abs(q302_named(p, N) - LIVE_Q302[p]) < 1e-6
        assert abs(q302_from_F(p, N) - LIVE_Q302[p]) < 1e-6
        assert abs(q302_named(p, N) - q302_drop_plus1(p, N)) > 1.0
        assert D["by_p"][str(p)]["match_Q302"] is True
        assert D["by_p"][str(p)]["drop_plus1_fails"] is True
    E = prove_open()
    assert E["proved"] is False
    assert E["Q3_double_named_in_p"] is False
    assert E["A_full_dbl_named_in_p"] is False
    assert E["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.564"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is True
    assert out["E"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15564.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
