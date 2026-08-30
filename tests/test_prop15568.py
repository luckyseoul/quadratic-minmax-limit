"""Tests for Prop 15.568 — named p=7 k-measure hits live A_full,T."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15567 import twochi_amp_named
from e1_gmin_m4_prop15568 import (
    fam_A3_named,
    fam_A6_named,
    interpolate_p_minus_5_over_15,
    main,
    mix_from_types,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    twochi_tilt_mean,
    want_A,
    want_B,
    want_mix,
)


def test_typeA_named_hits_and_uniform_fails():
    A = prove_A()
    assert A["proved"] is True
    assert A["match"] is True
    assert A["unif_fails"] is True
    assert A["both_signs_fail"] is True
    assert A["p11_fails"] is True
    assert abs(A["A"] - want_A(7)) < 1e-6
    assert abs(A["A3"] - 16.0 * 7**2) < 1e-6
    assert abs(A["A6"] - twochi_amp_named(7)) < 1e-6
    assert abs(A["t2"] - twochi_amp_named(7) / 9.0) < 1e-12
    assert A["nA3"] == 42
    assert A["nA6"] == 21
    assert len(fam_A3_named(7)) == 42
    assert len(fam_A6_named(7)) == 21


def test_typeB_named_hits_and_cuts_fail():
    B = prove_B()
    assert B["proved"] is True
    assert B["match"] is True
    assert B["unif_fails"] is True
    assert B["drop_ap_fails"] is True
    assert B["flip_mid_fails"] is True
    assert B["p11_fails"] is True
    assert abs(B["B"] - want_B(7)) < 1e-6
    assert abs(B["B4"] + 14896.0 / 27.0) < 1e-6
    assert abs(B["t2"] - twochi_tilt_mean(7, 1.0 / 3.0)) < 1e-12
    assert B["n_supp"] == 105
    assert B["sum_w"] == 189


def test_mix_and_flags():
    C = prove_C()
    assert C["proved"] is True
    assert abs(C["mix"] - want_mix(7)) < 1e-6
    assert abs(C["mix"] - mix_from_types(want_A(7), want_B(7))) < 1e-6
    assert C["p11_fails"] is True
    assert C["interp_fails"] is True
    assert abs(interpolate_p_minus_5_over_15(7) - want_mix(7)) > 1.0
    D = prove_open()
    assert D["proved"] is False
    assert D["A_full_dbl_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["prop"] == "15.568"
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "src" / "e1_gmin_m4_prop15568.py"
    text = src.read_text(encoding="utf-8")
    assert "phi_F_ge_6_proved_general" not in text
    assert "from e1_gmin_m4_prop15278" not in text
    assert "(p−5)/15" in text or "(p-5)/15" in text
    assert "Aut_e DEAD" in text
