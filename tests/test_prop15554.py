"""Tests for Prop 15.554 — ŷ support, Q3 line types, K_all split."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15554 import (
    LIVE_Q3,
    LIVE_SPLIT,
    main,
    prove_A,
    prove_B,
    prove_open,
    triple_kind,
)


def test_ypoint_supported_on_sum_zero():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["ok"] is True
        assert row["m1_max"] < 1e-6
        assert row["m2_off_max"] < 1e-6
        assert row["q3_nont_max"] < 1e-4
        assert row["q3_trip_max"] > 1.0
        assert abs(row["q3_nont_max"] - row["q3_trip_max"]) > 1.0


def test_Q3_constant_on_line_types():
    B = prove_B()
    assert B["proved"] is True
    assert triple_kind(5, 1) == "double"
    assert triple_kind(5, 6) == "n3"
    assert triple_kind(7, 1) == "double"
    assert triple_kind(7, 2) == "generic"
    assert triple_kind(7, 17) == "n3"
    r5 = B["by_p"]["5"]
    assert r5["double"]["n"] > 0
    assert abs(r5["double"]["mean"] - LIVE_Q3[5]["double"]) < 1e-4
    assert abs(r5["n3"]["mean"] - LIVE_Q3[5]["n3"]) < 1e-4
    assert r5["generic"]["n"] == 0
    assert abs(r5["n3"]["mean"] - r5["double"]["mean"]) > 1.0
    r7 = B["by_p"]["7"]
    assert abs(r7["double"]["mean"] - LIVE_Q3[7]["double"]) < 1e-4
    assert abs(r7["generic"]["mean"] - LIVE_Q3[7]["generic"]) < 1e-4
    assert abs(r7["n3"]["mean"] - LIVE_Q3[7]["n3"]) < 1e-4
    assert abs(r7["n3"]["mean"] - r7["double"]["mean"]) > 1.0
    assert r7["Q2"]["ok"] is True


def test_flags_untouched_bulk_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["split_certified"] is True
    assert C["K_all_named_in_p"] is False
    assert C["Q3_named_in_p"] is False
    assert C["Q4_is_Wick"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    for p, want in LIVE_SPLIT.items():
        rec = C["splits"][str(p)]
        assert rec["match"] is True
        assert rec["wick_fails"] is True
        assert abs(rec["Kall"] - want["Kall"]) < 1e-2
        assert abs(rec["WickOO"] - rec["bulk"]) > 1.0
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.554"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15554.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
