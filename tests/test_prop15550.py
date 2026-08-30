"""Tests for Prop 15.550 — S is Kloosterman; K_all rewrite not a p-law."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15550 import (
    S_named,
    S_of,
    S_wrong_drop4,
    S_wrong_lam,
    S_wrong_sign,
    kloosterman,
    main,
    pref_S,
    prove_A,
    prove_B,
    prove_open,
)


def test_S_is_kl_pref_one():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13, 17, 19):
        row = A["by_p"][str(p)]
        assert row["ok"] is True
        assert row["chi_m1"] == 1
        assert abs(row["pref_re"] - 1.0) < 1e-8
        assert abs(row["S0_re"] + 1.0) < 1e-8
        assert row["maxerr_S_eq_K"] < 1e-8
        assert abs(row["G1_re"] - row["G2_re"]) < 1e-8
    F = _kernel_field(5)
    assert abs(S_of(F, 0) - kloosterman(F, 1, 0)) < 1e-10
    assert abs(S_of(F, 0) + 1) < 1e-10
    assert abs(pref_S(F) - 1) < 1e-10
    assert abs(S_of(F, 3) - S_named(F, 3)) < 1e-10


def test_fail_when_wrong_dies():
    B = prove_B()
    assert B["proved"] is True
    F = _kernel_field(5)
    got = S_of(F, 2)
    assert abs(got - S_named(F, 2)) < 1e-10
    assert abs(got - S_wrong_sign(F, 2)) > 1.0
    assert abs(got - S_wrong_drop4(F, 2)) > 1.0
    assert abs(got - S_wrong_lam(F, 2)) > 1.0
    for p in (5, 7, 11, 13, 17, 19):
        row = B["by_p"][str(p)]
        assert row["ok"] is True
        for key in ("G", "mG", "chiG", "drop4", "lam", "sign"):
            assert row[key] > 1.0


def test_flags_untouched_Kall_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["K_all_named_in_p"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    assert C["J_named_ok"] is True
    assert C["rewrite_ok"] is True
    rec = C["rewrite"]
    assert rec["match_249050"] is True
    assert rec["term0_not_bulk"] is True
    assert rec["plus_fails"] is True
    assert abs(rec["term0"] - 650.0) < 1e-4
    assert abs(rec["K_all"] - 249050.0) < 1e-4
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.550"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    assert out["proved"]["K_all_named_in_p"] is False
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15550.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
