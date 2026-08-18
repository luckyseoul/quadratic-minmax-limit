"""Tests for Prop 15.557 — pair-sum of Φ_free is named."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15557 import (
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    sum_G3_1d_named,
    sum_G3_named,
    sum_G3_named_drop_2Np,
    sum_Phi_free_named,
)


def test_sum_G3_is_2N_over_p_minus_2pN():
    A = prove_A()
    assert A["proved"] is True
    assert abs(sum_G3_named(5, 130) + 1248.0) < 1e-9
    assert abs(sum_G3_named(7, 5726) + 78528.0) < 1e-9
    assert abs(sum_G3_named(5, 130) - sum_G3_named_drop_2Np(5, 130)) > 1.0
    for p in (5, 7):
        row = A["by_p"][str(p)]
        assert row["match"] is True
        assert row["drop_fails"] is True


def test_sum_G3_1d_is_named_binomial():
    B = prove_B()
    assert B["proved"] is True
    assert abs(sum_G3_1d_named(5) + 288.0) < 1e-9
    assert abs(sum_G3_1d_named(7) + 1920.0) < 1e-9
    m = 3
    want5 = -(5 - 1) * (5 + 1) ** 2 * comb(5, m) / 5
    assert abs(sum_G3_1d_named(5) - want5) < 1e-12
    assert abs(sum_G3_1d_named(5) - 0.0) > 1.0
    for p in (5, 7):
        row = B["by_p"][str(p)]
        assert row["match"] is True
        assert row["drop_C_fails"] is True


def test_sum_Phi_free_named_and_flags():
    C = prove_C()
    assert C["proved"] is True
    assert abs(sum_Phi_free_named(5, 130) + 960.0) < 1e-9
    assert abs(sum_Phi_free_named(7, 5726) + 76608.0) < 1e-9
    for p in (5, 7):
        row = C["by_p"][str(p)]
        assert row["match"] is True
        assert row["zero_fails"] is True
    D = prove_open()
    assert D["proved"] is False
    assert D["Phi_free_pointwise_named"] is False
    assert D["NUM_SUM_16pA_general"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.557"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15557.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
