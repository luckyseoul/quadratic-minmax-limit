"""Tests for Prop 15.86: ε(p)·n1 sum of star·τ1; value-set AP; residual budgets."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1586 import (  # noqa: E402
    B_cand,
    d1_one,
    epsilon,
    gain_L,
    gain_cand,
    n1_formula,
    prove_epsilon_algebra,
    prove_residual_budgets,
    prove_tau1_container,
    tau1_value_set_formula,
)


def test_epsilon_mod4_law():
    assert epsilon(3) == -1
    assert epsilon(5) == 1
    assert epsilon(7) == -1
    assert epsilon(11) == -1
    assert epsilon(13) == 1
    for p in range(3, 60, 2):
        if any(p % d == 0 for d in range(3, int(p**0.5) + 1)):
            continue
        assert epsilon(p) == (1 if p % 4 == 1 else -1)
        assert epsilon(p) == (-1) ** ((p - 1) // 2)


def test_tau1_value_set_algebra():
    # Certified sets from multi-W pure-C census
    assert tau1_value_set_formula(3) == [-1]
    assert tau1_value_set_formula(5) == [-1, 5]
    assert tau1_value_set_formula(7) == [-7, -1, 5]
    assert tau1_value_set_formula(11) == [-13, -7, -1, 5, 11]
    for p in (3, 5, 7, 11, 13, 17, 19):
        vals = tau1_value_set_formula(p)
        assert len(vals) == (p - 1) // 2
        assert all(v % 2 != 0 for v in vals)
        assert all((v % 6) == 5 for v in vals)
        # AP of difference 6
        for a, b in zip(vals, vals[1:]):
            assert b - a == 6
        assert d1_one(p) % 2 == 1


def test_residual_budgets_signs_and_limits():
    assert B_cand(5) == Fraction(-16, 325)
    assert B_cand(5) < 0
    for p in (7, 11, 13, 17, 19, 23):
        assert B_cand(p) > 0
        gL, gC = gain_L(p), gain_cand(p)
        assert 0 < gC < gL
        # gap identity from Prop 15.83
        assert gL - gC == Fraction(3 * (p - 2), 48 * (2 * p + 3))
    # B_cand → 1/2
    assert abs(float(B_cand(97)) - 0.5) < 0.03


def test_prove_helpers_and_evidence_json():
    eps = prove_epsilon_algebra()
    assert eps["proved_form"] is True
    assert "ε(p)" in eps["formula"] or "epsilon" in eps["formula"].lower() or "(-1)" in eps["formula"]

    primes = [p for p in range(5, 40) if all(p % d for d in range(2, int(p**0.5) + 1))]
    cont = prove_tau1_container(primes + [3])
    assert cont["proved_container_algebra"] is True

    bud = prove_residual_budgets(primes)
    assert bud["proved_algebra"] is True

    path = ROOT / "evidence" / "e1_gmin_m4_prop1586.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1586.py to produce evidence"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.86"
    assert data["L_status"] == "OPEN"
    assert data["proved"] is True
    assert data["census"]["all_match"] is True
    for p in ("3", "5", "7", "11"):
        row = data["census"]["by_p"][p]
        assert row["sum_eq_eps_n1"] is True
        assert row["value_set_match"] is True
        assert row["sum_star_tau1"] == epsilon(int(p)) * n1_formula(int(p))


def test_solution_and_handoff_open():
    sol = (ROOT / "solution.md").read_text()
    assert "15.86" in sol
    idx = sol.index("15.86")
    chunk = sol[idx : idx + 2500]
    assert "OPEN" in chunk
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "OPEN" in hand
    assert "15.86" in hand
