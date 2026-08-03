"""Tests for Prop 15.159 — Φ spectrum structure + 16N chain predicate."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15100 import d_of, n_of
from e1_gmin_m4_prop15159 import (
    SPECTRUM_P5,
    SPECTRUM_P7,
    certify_live_spectrum,
    dim_Z,
    dual_gap_eigenvalues,
    lambda_max_design_threshold,
    mu_bar,
    prove_dual_gap_G,
    prove_spectrum_p5,
    prove_spectrum_p7,
    prove_theorem_C,
    psl_irrep_dim_d,
    sixteen_N_from_dual_gap,
    sixteen_N_from_mult_d_and_kappa96,
)

ROOT = Path(__file__).resolve().parents[1]


def test_threshold_algebra_below_16_and_above_mu():
    """Shipped Fraction algebra: 16(d−2)/d < 16 and ≥ μ̄ for p≥5."""
    r = prove_theorem_C([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    for p in (5, 7, 11, 13):
        thr = lambda_max_design_threshold(p)
        assert thr < 16
        assert thr >= mu_bar(p)
        d = d_of(p)
        assert thr == Fraction(16 * (d - 2), d)


def test_spectrum_p5_exact_forms():
    """Exact p=5 spectrum table: λ_max=16(d−2)/d, mult=d, mean=μ̄."""
    r = prove_spectrum_p5()
    assert r["mult_sum_ok"] is True
    assert r["equals_16_dminus2_over_d"] is True
    assert r["mult_top_eq_d"] is True
    assert r["mean_eq_mu_bar"] is True
    assert max(SPECTRUM_P5) == Fraction(176, 13)
    assert SPECTRUM_P5[Fraction(176, 13)] == 13
    assert sum(SPECTRUM_P5.values()) == dim_Z(5)


def test_spectrum_p7_exact_forms():
    """Exact p=7 spectrum: den=N/(4p), mult(λ_max)=d, λ_max<16."""
    r = prove_spectrum_p7()
    assert r["mult_sum_ok"] is True
    assert r["mult_top_eq_d"] is True
    assert r["lambda_max_lt_16"] is True
    assert r["lambda_max_lt_threshold"] is True
    assert r["den_N_over_4p"] == 409
    assert max(SPECTRUM_P7) == Fraction(4320, 409)
    assert sum(SPECTRUM_P7.values()) == dim_Z(7)


def test_psl_irrep_dim_matches_d():
    """PSL(2,p²) has irrep dimension d=(p²+1)/2."""
    for p in (3, 5, 7, 11, 13):
        assert psl_irrep_dim_d(p) == d_of(p)
        assert psl_irrep_dim_d(p) == n_of(p) // 2


def test_sixteen_N_chain_predicate_real_entry():
    """
    Drive shipped sixteen_N_from_mult_d_and_kappa96.
    With mult≥d and κ²≤96n (certified at p=5,7 by residual census), 16N holds.
    Without κ certificate, returns false — not hard-coded success.
    """
    ok5 = sixteen_N_from_mult_d_and_kappa96(13, 13, True, 5)
    assert ok5["sixteen_N"] is True
    assert ok5["bi_tight_empty_via_15_61"] is True

    ok7 = sixteen_N_from_mult_d_and_kappa96(25, 25, True, 7)
    assert ok7["sixteen_N"] is True

    # Missing κ bound ⇒ chain does not fire
    fail = sixteen_N_from_mult_d_and_kappa96(13, 13, False, 5)
    assert fail["sixteen_N"] is False

    # mult too small
    fail2 = sixteen_N_from_mult_d_and_kappa96(12, 13, True, 5)
    assert fail2["sixteen_N"] is False


def test_dual_gap_G_structure():
    """
    Theorem G: dual gap G=(d/32)(16I−Φ) has eigs {1,2,4} at p=5;
    G ≽ I at p=5,7; G≽I ⇒ 16N predicate.
    """
    g = prove_dual_gap_G()
    assert g["certified_p5_p7"] is True
    assert g["proved_for_all_p"] is False
    assert g["p5"]["equals_1_2_4"] is True
    assert g["p5"]["G_ge_I"] is True
    assert g["p7"]["G_ge_I"] is True
    assert g["p7"]["strictly_above_4"] is True

    g5 = dual_gap_eigenvalues(SPECTRUM_P5, 13)
    assert g5 == {Fraction(1): 13, Fraction(2): 26, Fraction(4): 26}

    assert sixteen_N_from_dual_gap(Fraction(1), 5)["sixteen_N"] is True
    assert sixteen_N_from_dual_gap(Fraction(3475, 818), 7)["sixteen_N"] is True
    assert sixteen_N_from_dual_gap(Fraction(1, 2), 5)["sixteen_N"] is False


@pytest.mark.parametrize("p", [5, 7])
def test_live_phi_spectrum_matches_structure(p: int):
    """Live Max+ eigensystem: mult(λ_max)=d and λ_max matches table."""
    r = certify_live_spectrum(p)
    if r.get("skipped"):
        pytest.skip(r.get("reason", "no Max+ cache"))
    assert r["mult_top_eq_d"] is True
    assert r["lambda_max_lt_16"] is True
    if p == 5:
        assert abs(r["lambda_max"] - float(Fraction(176, 13))) < 1e-8
        assert r["lambda_max_le_threshold"] is True
    if p == 7:
        assert abs(r["lambda_max"] - float(Fraction(4320, 409))) < 1e-8


def test_evidence_json_honest_open():
    """Evidence file must not claim general residual closed."""
    path = ROOT / "evidence" / "e1_gmin_m4_prop15159.json"
    if not path.is_file():
        from e1_gmin_m4_prop15159 import main

        main()
    e = json.loads(path.read_text())
    assert e["L_status"] == "OPEN"
    assert e["proved"]["residual_for_all_p_ge_5"] is False
    assert e["status"]["general_residual"] is False
    assert e["proved"]["spectrum_p5"] is True
    assert e["proved"]["spectrum_p7"] is True
