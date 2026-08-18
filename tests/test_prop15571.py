"""Tests for Prop 15.571 — A6 expansions; amp guesses die."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import gauss_table, prod_dbl_k
from e1_gmin_m4_prop15571 import (
    a6_deg_named,
    a6_guess_3mchi2,
    a6_guess_pm3_over4,
    a6_nondeg_named,
    main,
    pair_k,
    prove_A,
    prove_B,
    prove_open,
)


def test_a6_expansions_match_and_p5_3pt_fails():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["fail_3pt_on_deg"] is True
    G, omega = gauss_table(7)
    tau = G[1]
    got = prod_dbl_k(pair_k(7, {1, 2, 3}), 7, G, omega)
    assert abs(got - a6_nondeg_named(7, 1, 2, tau)) < 1e-8
    G5, omega5 = gauss_table(5)
    tau5 = G5[1]
    got5 = prod_dbl_k(pair_k(5, {1, 4}), 5, G5, omega5)
    assert abs(got5 - a6_deg_named(5, 1, tau5)) < 1e-8
    assert abs(got5 - a6_nondeg_named(5, 1, 4, tau5)) > 1.0


def test_amp_guesses_die_and_flags():
    B = prove_B()
    assert B["proved"] is True
    assert abs(B["by_p"]["7"]["live"] - a6_guess_3mchi2(7)) < 1e-6
    assert abs(B["by_p"]["23"]["live"] - a6_guess_3mchi2(23)) > 1.0
    assert abs(B["by_p"]["19"]["live"] - a6_guess_pm3_over4(19)) > 1.0
    C = prove_open()
    assert C["proved"] is False
    assert C["A6_named_in_p"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.571"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15571.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
