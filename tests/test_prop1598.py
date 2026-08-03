"""Tests for Prop 15.98: mult≥d−1 (PSL); gap criterion with d−1."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1598 import (  # noqa: E402
    d_of,
    n_of,
    prove_gap_criterion_d_minus_1,
    prove_mult_ge_d_minus_1,
    prove_psl_min_irrep,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_psl_min_irrep():
    out = prove_psl_min_irrep(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        d = d_of(p)
        assert out["by_p"][str(p)]["min_nontrivial_irrep_dim"] == d - 1
        assert out["by_p"][str(p)]["equals_d_minus_1"] is True


def test_mult_ge_d_minus_1_theorem():
    out = prove_mult_ge_d_minus_1()
    assert out["proved"] is True
    assert "d−1" in out["theorem"] or "d-1" in out["theorem"]


def test_gap_algebra_d_minus_1():
    out = prove_gap_criterion_d_minus_1(_primes(40))
    assert out["proved"] is True
    # p=3 fails, p>=5 holds
    assert out["by_p"]["3"]["gap_algebra_holds"] is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert out["by_p"][str(p)]["gap_algebra_holds"] is True
        d = d_of(p)
        assert d * d - 9 * d - 24 > 0
        # exact Fraction check
        n = n_of(p)
        wick = 12 * n * n + 48 * n
        lhs = Fraction(wick, 16) - d * d
        rhs = Fraction((d - 1) * d * d, 4)
        assert lhs <= rhs


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1598.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1598.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.98"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["mult_ge_d_minus_1_proved_for_paley"] is True
    assert data["proved"] is True

    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["mult_ge_d_minus_1"] is True
    assert c5["mult_ge_d_minus_1"] is True
    assert c3["gap_by_mult_d_minus_1"] is False
    assert c5["gap_by_mult_d_minus_1"] is True
    assert c5["kappa2_le_96n"] is True


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.98" in sol
    idx = sol.index("15.98")
    assert "OPEN" in sol[idx : idx + 5000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "15.98" in hand
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.98" in graph
