"""Tests for Prop 15.114 — γ-calculus; Tf_y; ∑γ,∑γ²; ED4 types."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15114 import (  # noqa: E402
    Tf_energy_formula,
    certify_ed4_types,
    certify_gamma_and_Tf,
    is_prime,
    main,
    prove_gamma_algebra,
    sum_gamma2_formula,
    sum_gamma_formula,
)


def test_gamma_algebra_closed_forms():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_gamma_algebra(primes)
    assert r["proved"] is True
    # spot checks
    assert sum_gamma_formula(5) == Fraction(6, 5) * 14950
    assert sum_gamma2_formula(3) == Fraction(1440)
    assert sum_gamma2_formula(5) == Fraction(93600)
    assert sum_gamma2_formula(7) == Fraction(1411200)
    # energy consistency at p=5
    e = Tf_energy_formula(5)
    assert e > 0


def test_cert_gamma_Tf_p3_p5():
    for p in (3, 5):
        c = certify_gamma_and_Tf(p, n_y=2)
        if c.get("skipped"):
            continue
        assert c["all_ok"] is True
        for s in c["samples"]:
            assert s["Tf_multiplicative_ok"] is True
            assert s["sum_gamma_ok"] is True
            assert s["sum_gamma2_ok"] is True
            assert s["Tf_energy_ok"] is True


def test_ed4_types_and_residual_cert():
    for p in (3, 5, 7):
        e = certify_ed4_types(p)
        if e.get("skipped"):
            continue
        assert e["global_le_suf"] is True
        assert e["delta2_le_rho_min2"] is True
        assert e["max_y_le_suf"] is True
        if p == 7:
            # multi-type is the structural finding
            assert e["not_2point_homogeneous"] is True
            assert e["n_types"] == 3
        if p in (3, 5):
            # homogeneous at small p
            assert e["n_types"] == 1


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Tf_multiplicative"] is True
    assert out["proved"]["sum_gamma2"] is True
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15114.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.114"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
