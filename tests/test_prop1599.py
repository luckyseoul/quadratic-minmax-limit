"""Tests for Prop 15.99: κ-budget structure; min-distance; closed forms."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1599 import (  # noqa: E402
    C_diag_formula,
    d_of,
    n_of,
    n3_of,
    prove_budget_algebra,
    prove_closed_forms,
    prove_min_distance,
    prove_master_source,
    rho2_budget,
    sum_kappa,
    sum_kappa2,
    sum_m4,
    sum_m4_kappa,
    sum_rho,
    sum_rho_kappa,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_min_distance_algebra():
    out = prove_min_distance(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        row = out["by_p"][str(p)]
        assert row["min_dH"] == p + 1
        assert row["equals_formula"] is True
        n = n_of(p)
        assert row["max_abs_dot_non_antipodal"] == n - 2 * (p + 1)


def test_budget_algebra():
    out = prove_budget_algebra(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13):
        row = out["by_p"][str(p)]
        assert row["C_diag_le_96n"] is True
        assert row["room_over_24_eq_rho2_budget"] is True
        n = n_of(p)
        cd = C_diag_formula(p)
        room = Fraction(96 * n) - cd
        assert room / 24 == rho2_budget(p)
        assert room == Fraction(4 * n * (13 * p * p + 3), p * p)


def test_closed_forms_consistency():
    out = prove_closed_forms(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11):
        # ∑ρ and ∑ρκ match definitions
        assert sum_rho(p) == sum_m4(p) - sum_kappa(p) / Fraction(p * p)
        assert sum_rho_kappa(p) == sum_m4_kappa(p) - sum_kappa2(p) / Fraction(
            p * p
        )
        n = n_of(p)
        n4 = Fraction(n * (n - 1) * (n - 2) * (n - 3), 24)
        n1 = Fraction(n * (n - 1) * (n - 2) * (n - 2), 32)
        assert n1 + n3_of(p) == n4
        # known small values
        if p == 3:
            assert sum_kappa(p) == 18
            assert sum_rho(p) == -16
            assert sum_rho_kappa(p) == 40
        if p == 5:
            assert sum_kappa(p) == 150
            assert sum_rho(p) == -96
            assert sum_rho_kappa(p) == 312


def test_master_source():
    out = prove_master_source(_primes(30))
    assert out["proved"] is True
    for p in (3, 5, 7):
        n3 = n3_of(p)
        b2 = Fraction(576) * n3 / Fraction(p**4)
        assert out["by_p"][str(p)]["rhs_Tkappa_over_p2_sq_norm"] == str(b2)


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1599.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1599.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.99"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["min_distance"]["proved"] is True
    assert data["budget_algebra"]["proved"] is True
    assert data["closed_forms"]["proved"] is True

    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["ok"] is True
    assert c5["ok"] is True
    assert c3["kappa2_le_96n"] is True
    assert c5["kappa2_le_96n"] is True
    # equality only p=3
    assert abs(c3["kappa2"] - 960.0) < 1e-4
    assert c5["kappa2"] < 96 * 26 - 1
    assert c3["min_dH_ge_p_plus_1"] is True
    assert c5["min_dH_ge_p_plus_1"] is True
    assert c5["closed_forms_match_census"] is True


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.99" in sol
    idx = sol.index("15.99")
    assert "OPEN" in sol[idx : idx + 6000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.99" in graph
