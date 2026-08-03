"""Tests for Prop 15.108 — residual-Gram/Schur; algebraic Thm A; Parseval T_ρ."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15105 import room_hyp_of  # noqa: E402
from e1_gmin_m4_prop15108 import (  # noqa: E402
    certify_parseval,
    certify_popp_phi_conversion,
    is_prime,
    kappa2_stratum,
    main,
    prove_m4_parseval_form,
    prove_parseval_target,
    prove_residual_gram_identities,
    prove_theorem_A_polynomial,
    rho_target,
    thmA_gap_rational,
)


def test_thmA_gap_closed_form_matches_defs():
    primes = [p for p in range(3, 50) if is_prime(p)]
    r = prove_theorem_A_polynomial(primes)
    assert r["proved"] is True
    for p in primes:
        assert r["by_p"][str(p)]["match"] is True
        assert r["by_p"][str(p)]["nonneg"] is True


def test_thmA_gap_factorization_at_samples():
    # 128(p-3)(p+3)(p^4-12p^2-5)/((p^2-5)^2 (p^2+1)^2)
    for p in [5, 7, 11, 13, 17, 19, 23]:
        num = 128 * (p - 3) * (p + 3) * (p**4 - 12 * p * p - 5)
        den = (p * p - 5) ** 2 * (p * p + 1) ** 2
        assert thmA_gap_rational(p) == Fraction(num, den)
        assert thmA_gap_rational(p) > 0
    assert thmA_gap_rational(3) == 0


def test_parseval_target_closed_forms():
    primes = [p for p in range(3, 40) if is_prime(p)]
    r = prove_parseval_target(primes)
    assert r["proved"] is True
    for p in primes:
        T = rho_target(p)
        assert T == rho_min2(p) + delta2_budget_hyp(p)
        if p > 3:
            assert delta2_budget_hyp(p) == room_hyp_of(p) / 24
            assert T > rho_min2(p)


def test_residual_gram_flat_identity():
    primes = [p for p in range(5, 40) if is_prime(p)]
    r = prove_residual_gram_identities(primes)
    assert r["proved"] is True


def test_m4_kappa_stratum():
    primes = [p for p in range(3, 30) if is_prime(p)]
    r = prove_m4_parseval_form(primes)
    assert r["proved"] is True
    for p in primes:
        n = p * p + 1
        assert kappa2_stratum(p) == Fraction(
            n * (n - 1) * (n - 2) * (n - 5), 8
        )


def test_cert_parseval_p3_p5():
    for p in (3, 5):
        c = certify_parseval(p)
        if c.get("skipped"):
            continue
        assert c["rho2_le_T"] is True
        assert c["orth_le_room_hyp"] is True


def test_cert_parseval_p7_if_available():
    c = certify_parseval(7)
    if c.get("skipped"):
        return
    assert c["rho2_le_T"] is True
    assert c["orth_le_room_hyp"] is True
    # strict at p=7
    ratio = Fraction(c["ratio_orth_over_room_hyp"])
    assert ratio < 1
    assert ratio > 0


def test_popp_phi_conversion_p5():
    c = certify_popp_phi_conversion(5)
    if c.get("skipped"):
        return
    assert c["lambda2_le_4_over_N"] is True
    assert c["L_phi_le_16"] is True
    # exact: λ2=11/845, 4N λ2=176/13
    assert abs(c["four_N_lambda2"] - 176 / 13) < 1e-6


def test_general_status_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["theorem_A_polynomial"] is True
    assert out["proved"]["parseval_target_formula"] is True
    assert out["proved"]["sum_rho2_le_T_for_all_p"] is False
    assert out["proved"]["sixteen_N_for_all_p"] is False
    assert out["general_status"]["orth_le_room_hyp_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15108.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.108"
    assert data["L_status"] == "OPEN"
