"""Tests for Prop 15.567 — 2χ double product is a complete Gauss sum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import (
    chip,
    gauss_table,
    main,
    prod_dbl_k,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    twochi_amp_drop_m1,
    twochi_amp_named,
    twochi_amp_wrong_8p2_over_3,
    twochi_named,
)


def test_zhat_is_two_G_khat():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["drop2_fails"] is True


def test_twochi_amp_is_named_complete_sum():
    B = prove_B()
    assert B["proved"] is True
    assert abs(twochi_amp_named(5)) < 1e-12
    assert abs(twochi_amp_named(13)) < 1e-12
    assert abs(twochi_amp_named(7) - 8.0 * 49.0 / 3.0) < 1e-12
    assert abs(twochi_amp_named(11) + 24.0 * 121.0 / 5.0) < 1e-12
    # fail 8p²/3 at p=11
    assert abs(twochi_amp_named(11) - twochi_amp_wrong_8p2_over_3(11)) > 1.0
    # fail drop (1−χ(−1)) at p=5
    assert abs(twochi_amp_drop_m1(5)) > 1.0
    G, omega = gauss_table(7)
    k = [0] * 7
    k[0] = 1
    k[3] = -1  # χ(3)=−1
    got = prod_dbl_k(tuple(k), 7, G, omega)
    assert abs(got.real - twochi_named(7, 3)) < 1e-8
    assert abs(got.real - twochi_amp_named(7)) > 1.0  # drop χ
    assert chip(7, 3) == -1


def test_live_k_and_flags():
    C = prove_C()
    assert C["proved"] is True
    assert C["match"] is True
    assert C["uniform_fails"] is True
    D = prove_open()
    assert D["proved"] is False
    assert D["A_full_dbl_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.567"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15567.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
