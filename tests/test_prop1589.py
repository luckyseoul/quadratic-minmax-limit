"""Tests for Prop 15.89: Wick/Q4 split, κ_C·κ_B, H equivalence."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1589 import (  # noqa: E402
    H,
    certify_kappaB2_formula,
    prove_settlement,
    prove_wick_split_algebra,
    residual_budget_sum_rho_kappa,
    wick_Q4_over_N,
)


def test_wick_Q4_formula():
    assert wick_Q4_over_N(5) == Fraction(54, 25)  # 2+4/25=54/25
    assert wick_Q4_over_N(5) == Fraction(2) + Fraction(4, 25)
    for p in (3, 5, 7, 11, 13):
        w = wick_Q4_over_N(p)
        assert w == Fraction(2 * (p * p + 2), p * p)
        # Wick alone below 2H
        assert w < 2 * H(p)


def test_H_equivalence_budget():
    # p=5: (49/13 - 1 - 2/25)/4 = 437/650
    assert residual_budget_sum_rho_kappa(5) == Fraction(437, 650)
    for p in (3, 5, 7, 11, 13, 17, 19):
        rb = residual_budget_sum_rho_kappa(p)
        assert rb > 0
        # room consistency: 8 * rb = 2H - wick
        assert 8 * rb == 2 * H(p) - wick_Q4_over_N(p)


def test_wick_split_algebra_proved():
    primes = []
    for p in range(3, 40):
        if p % 2 == 0 and p > 2:
            continue
        if all(p % d != 0 for d in range(3, int(p**0.5) + 1, 2)):
            primes.append(p)
    out = prove_wick_split_algebra(primes)
    assert out["proved"] is True


def test_kappaB2_and_settlement():
    kb2 = certify_kappaB2_formula()
    assert bool(kb2["ok"])
    s = prove_settlement()
    assert s["proved_form"] is True


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1589.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.89"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False  # honest: residual still open
    assert data["proved"] is True  # structure identities
    assert data["kappaC_kappaB_census"]["all_ok"] is True
    for p in ("3", "5", "7", "11"):
        row = data["kappaC_kappaB_census"]["by_p"][p]
        assert row["ok"] is True
        n = int(row["n"])
        assert abs(row["c"] - (n + 1) / 4) < 1e-8


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.89" in sol
    idx = sol.index("15.89")
    assert "OPEN" in sol[idx : idx + 4000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
