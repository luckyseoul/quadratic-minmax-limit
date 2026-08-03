"""Tests for Prop 15.102: resolvent δ-calculus; orth=24‖δ‖²."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, kappa_hyp, n_of, proj_kappa2, room_orth  # noqa: E402
from e1_gmin_m4_prop15102 import (  # noqa: E402
    C_diag,
    b2_of,
    delta2_budget_96n,
    delta2_budget_hyp,
    mu2_of,
    n3_of,
    prove_invertible_case,
    prove_kappa_min_eq_proj,
    prove_orth_delta_identity,
    prove_reduction,
    prove_rho_min_closed_form,
    prove_source_and_mu2_form,
    rho_min2,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_source_mu2_form():
    out = prove_source_and_mu2_form(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        assert mu2_of(p) == 4 * (p * p + 15)
        assert 16 * p * p - mu2_of(p) == 12 * (p * p - 5)
        assert 16 * p * p + mu2_of(p) == 20 * (p * p + 3)
        n = n_of(p)
        assert n3_of(p) == Fraction(n * (n - 1) * (n - 2) * (n - 6), 96)
        assert b2_of(p) == Fraction(576) * n3_of(p) / Fraction(p**4)


def test_rho_min_closed_form():
    out = prove_rho_min_closed_form(_primes(40))
    assert out["proved"] is True
    assert rho_min2(3) == Fraction(200, 9)
    assert rho_min2(5) == Fraction(728, 25)
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        claimed = Fraction(5 * n * (p * p - 1) * (p * p + 3), 6 * p * p * (p * p - 5))
        assert rho_min2(p) == claimed
        # balanced formula
        b2 = b2_of(p)
        mu2 = mu2_of(p)
        bal = b2 * Fraction(16 * p * p + mu2, (16 * p * p - mu2) ** 2)
        assert rho_min2(p) == bal


def test_kappa_min_eq_proj():
    out = prove_kappa_min_eq_proj(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert C_diag(p) + 24 * rho_min2(p) == proj_kappa2(p)


def test_orth_delta_budgets():
    out = prove_orth_delta_identity(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13):
        n = n_of(p)
        d = d_of(p)
        assert delta2_budget_96n(p) == room_orth(p) / 24
        assert delta2_budget_96n(p) == Fraction(4 * n * (p * p - 9), 3 * (p * p - 5))
        assert delta2_budget_hyp(p) == room_orth(p) / 24 * Fraction((d - 1) ** 2, d * d)
        # kappa_hyp = proj + 24 * delta2_budget_hyp
        assert proj_kappa2(p) + 24 * delta2_budget_hyp(p) == kappa_hyp(p)
        assert proj_kappa2(p) + 24 * delta2_budget_96n(p) == 96 * n


def test_invertible_and_reduction():
    assert prove_invertible_case()["proved"] is True
    red = prove_reduction()
    assert red["proved"] is True
    assert red["general_orth_inequality_proved"] is False
    assert red["kappa_bound_proved_for_all_p"] is False


def test_p5_delta_budget_exact():
    # room_hyp/24 = 1536/65 at p=5
    assert delta2_budget_hyp(5) == Fraction(1536, 65)
    assert delta2_budget_96n(5) == Fraction(416, 15)  # 27.733...


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop15102.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop15102.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.102"
    assert data["L_status"] == "OPEN"
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["general_orth_inequality_proved"] is False
    assert data["source_mu2"]["proved"] is True
    assert data["rho_min_closed_form"]["proved"] is True
    assert data["kappa_min_eq_proj"]["proved"] is True
    assert data["orth_delta_identity"]["proved"] is True
    assert data["proved"] is True
    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["ok"] is True
    assert c5["ok"] is True
    assert c3["delta2_eq_0"] is True
    assert c5["delta2_eq_hyp_budget"] is True
    assert c5["mult_E4p"] == 12
    assert c5["mult_E4p_eq_d_minus_1"] is True


def test_solution_handoff_graph():
    assert "15.102" in (ROOT / "solution.md").read_text()
    assert "15.102" in (ROOT / "HANDOFF.md").read_text()
    assert "15.102" in (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
