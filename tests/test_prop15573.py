"""Tests for Prop 15.573 — Q_sub exclusive 1D/k=3/full mix."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15573 import main, prove_A, prove_open, q_sub_mix


def test_exclusive_mix_hits_live_Q_sub():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["dup_fails"] is True
        hit = q_sub_mix(p)
        assert abs(hit["mix"] - hit["live"]) < 1e-9
        assert "1d" in hit["parts"]
        assert "k3" in hit["parts"]
    dup = q_sub_mix(5, double_count=True)
    hit5 = q_sub_mix(5)
    assert abs(dup["mix"] - hit5["live"]) > 1.0
    assert "full_dup" in dup["parts"]


def test_flags_and_Q_tau_still_open():
    B = prove_open()
    assert B["proved"] is False
    assert B["Q_tau_named_in_p"] is False
    assert B["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.573"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15573.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p-5)/15" in text or "(p−5)/15" in text
