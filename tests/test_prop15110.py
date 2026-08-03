"""Tests for Prop 15.110 — closed Max+ identities; ρ_min < budget p≥7."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15110 import (  # noqa: E402
    certify_delta_le_rho_min,
    e4_formula,
    e4_newton,
    is_prime,
    main,
    prove_e4_formula,
    prove_rho_min_lt_budget,
    prove_sum_kappa_prod_algebra,
    rho_min_minus_budget,
    sum_kappa_prod_formula,
)


def test_sum_kappa_prod_algebra():
    primes = [p for p in range(3, 50) if is_prime(p)]
    r = prove_sum_kappa_prod_algebra(primes)
    assert r["proved"] is True


def test_e4_formula_matches_newton():
    primes = [p for p in range(3, 40) if is_prime(p)]
    r = prove_e4_formula(primes)
    assert r["proved"] is True
    for p in primes:
        n = p * p + 1
        assert e4_newton(p + 1, n) == e4_formula(p)
        assert e4_newton(-(p + 1), n) == e4_formula(p)


def test_rho_min_lt_budget_all_primes_ge_7():
    primes = [p for p in range(7, 100) if is_prime(p)]
    r = prove_rho_min_lt_budget(primes)
    assert r["proved"] is True
    for p in primes:
        assert rho_min_minus_budget(p) < 0
        assert rho_min2(p) < delta2_budget_hyp(p)


def test_rho_min_gt_budget_at_p5():
    # Equality case: ρ_min² > budget at p=5 (δ takes the rest to equality)
    assert rho_min2(5) > delta2_budget_hyp(5)
    assert rho_min_minus_budget(5) > 0


def test_closed_forms_spot():
    assert sum_kappa_prod_formula(5) == Fraction(1950)
    assert e4_formula(5) == Fraction(-90)
    assert e4_formula(7) == Fraction(-308)


def test_cert_delta_le_rm_p3_p5_p7():
    for p in (3, 5, 7):
        c = certify_delta_le_rho_min(p)
        if c.get("skipped"):
            continue
        assert c["delta2_le_rho_min2"] is True
        assert c["delta2_le_budget"] is True


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["sum_kappa_prod"] is True
    assert out["proved"]["rho_min_lt_budget_p_ge_7"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["sum_rho2_le_T_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15110.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.110"
    assert data["L_status"] == "OPEN"
