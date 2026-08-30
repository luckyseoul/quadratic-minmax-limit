"""Tests for Prop 15.572 — A6 family mean is a Jacobi sum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import twochi_amp_named
from e1_gmin_m4_prop15571 import a6_guess_3mchi2
from e1_gmin_m4_prop15572 import (
    a6_drop_chi3_deg,
    a6_named_closed,
    a6_t_minus_12_everywhere,
    deg_p2_named,
    main,
    prove_A,
    prove_open,
    t_gen_named,
)


def test_a6_closed_hits_and_guesses_die():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13, 17, 23, 29, 73):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert abs(row["got"] - float(a6_named_closed(p))) < 1e-6
    from fractions import Fraction

    assert a6_named_closed(5) == 1000
    assert a6_named_closed(7) == Fraction(392, 3)
    assert abs(float(a6_named_closed(7)) - twochi_amp_named(7)) < 1e-9
    assert abs(float(a6_named_closed(23)) - a6_guess_3mchi2(23)) > 1.0
    assert abs(float(a6_named_closed(23)) - float(a6_t_minus_12_everywhere(23))) > 1.0
    assert abs(float(a6_named_closed(13)) - float(a6_drop_chi3_deg(13))) > 1.0
    assert A["by_p"]["23"]["amp_guess_fails_or_ok"] is True
    assert A["by_p"]["13"]["chi3_fails_or_ok"] is True
    assert t_gen_named(11) == -12
    assert t_gen_named(7) == -2
    assert deg_p2_named(5) == 40


def test_flags_and_open():
    B = prove_open()
    assert B["proved"] is False
    assert B["A_full_dbl_named_in_p"] is False
    assert B["Q_tau_named"] is False
    assert B["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.572"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15572.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "Aut_e DEAD" in text
    assert "amp(3−χ(2))/2" in text or "amp(3-χ(2))/2" in text
