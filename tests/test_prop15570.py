"""Tests for Prop 15.570 — A3 family mean is a Jacobi sum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import twochi_amp_named
from e1_gmin_m4_prop15570 import (
    a3_drop_chi2,
    a3_named_closed,
    a3_only_6amp,
    a3_wrong_16p2,
    a6_amp_times_p_minus_3_over_4,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_a3_closed_hits_and_16p2_dies():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13, 17, 19):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert abs(row["got"] - a3_named_closed(p)) < 1e-6
    assert abs(a3_named_closed(7) - 16.0 * 49.0) < 1e-9
    assert abs(a3_named_closed(7) - 6.0 * twochi_amp_named(7)) < 1e-9
    assert abs(a3_named_closed(5) - a3_wrong_16p2(5)) > 1.0
    assert abs(a3_named_closed(11) - a3_wrong_16p2(11)) > 1.0
    assert abs(a3_named_closed(5) - a3_only_6amp(5)) > 1.0
    assert abs(a3_named_closed(11) - a3_drop_chi2(11)) > 1.0
    assert A["by_p"]["5"]["sixteen_fails"] is True
    assert A["by_p"]["11"]["sixteen_fails"] is True
    assert A["by_p"]["5"]["six_amp_fails_or_p3"] is True
    assert A["by_p"]["11"]["drop_chi2_fails_or_ok"] is True


def test_a6_guess_dies_and_flags():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["by_p"]["7"]["err"]) < 1e-6
    assert abs(B["by_p"]["11"]["err"]) < 1e-6
    assert B["by_p"]["19"]["guess_dies_p19"] is True
    assert abs(B["by_p"]["19"]["A6"] - a6_amp_times_p_minus_3_over_4(19)) > 1.0
    C = prove_open()
    assert C["proved"] is False
    assert C["A6_named_in_p"] is False
    assert C["A_full_dbl_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.570"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15570.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "Aut_e DEAD" in text
    assert "16p²" in text or "16p2" in text
