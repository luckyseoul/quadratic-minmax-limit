"""Tests for Prop 15.555 — G3 is Φ of F_p-norms; Q3 Fourier of Φ."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15555 import main, prove_A, prove_B, prove_open


def test_G3_depends_only_on_norms():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["ok"] is True
        assert row["n_split"] == 0
        assert row["n_constant"] == row["n_keys"]
        assert row["chi_split"] > 0
        assert row["chi_max_nuniq"] >= 2
    assert A["by_p"]["5"]["nuniq_Phi"] == 7
    assert A["by_p"]["7"]["nuniq_Phi"] == 14


def test_Q3_is_collision_plus_Phi_fourier():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7):
        row = B["by_p"][str(p)]
        assert row["Phi_rebuild_maxerr"] < 1e-6
        assert row["ok"] is True
        for kind, rec in row["samples"].items():
            assert rec["match"] is True
            assert rec["drop_fails"] is True
            if rec["live"] is not None:
                assert abs(rec["raw"] - rec["live"]) < 1e-2


def test_flags_untouched_Q3_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Phi_named_in_p"] is False
    assert C["Q3_named_in_p"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.555"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15555.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
