"""Tests for Prop 15.117 — Path C hyp residual; ρ_min pairings."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import (  # noqa: E402
    b_dot_f_formula,
    certify_bgf_and_pairings,
    certify_hyp_residual,
    is_prime,
    main,
    prove_hyp_budget_form,
    prove_pairings,
    prove_rho_vs_hyp,
    prove_slack_identity,
    rho_min_dot_f_if_bgf_zero,
    rho_min_dot_m4,
    room_hyp_over_24,
)


def test_hyp_budget_closed_form():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_hyp_budget_form(primes)
    assert r["proved"] is True
    assert room_hyp_over_24(5) == Fraction(1536, 65)
    assert room_hyp_over_24(5) == delta2_budget_hyp(5)
    assert room_hyp_over_24(7) == delta2_budget_hyp(7)


def test_rho_vs_hyp_and_slack():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_rho_vs_hyp(primes)
    assert r["proved"] is True
    assert rho_min2(5) > room_hyp_over_24(5)
    assert rho_min2(7) < room_hyp_over_24(7)
    s = prove_slack_identity(primes)
    assert s["proved"] is True


def test_pairings_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_pairings(primes)
    assert r["proved"] is True
    assert b_dot_f_formula(5) == Fraction(2 * (5**4 - 1), 5)
    assert rho_min_dot_f_if_bgf_zero(5) == Fraction(208, 5)
    assert rho_min_dot_f_if_bgf_zero(7) == Fraction(800, 11)
    assert rho_min_dot_m4(5) > 0


def test_cert_and_main_open():
    c3 = certify_bgf_and_pairings(3)
    if not c3.get("skipped"):
        assert c3["bgf_pointwise_zero"] is True
        assert c3["bf_constant_match"] is True
        assert c3["rmf_match_closed"] is True

    for p in (3, 5, 7):
        h = certify_hyp_residual(p)
        if h.get("skipped"):
            continue
        assert h["delta2_le_hyp"] is True
        assert h["delta2_le_rho_min2"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["hyp_budget_form"] is True
    assert out["proved"]["pairings"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    assert out["proved"]["bgf_pointwise_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15117.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.117"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
