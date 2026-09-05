"""Tests for Prop 15.96: Wick–κ calculus; ‖κ‖²≤96n ⇔ ∑M²≤Wick."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1596 import (  # noqa: E402
    C_diag_formula,
    d_of,
    n_of,
    prove_gap_via_kappa_budget,
    prove_wick_kappa_calculus,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_wick_kappa_algebra():
    out = prove_wick_kappa_calculus(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        n = n_of(p)
        row = out["by_p"][str(p)]
        assert row["Wick_hi"] == 12 * n * n + 48 * n
        assert row["Wick_lo"] == 12 * n * n - 48 * n
        assert row["inner_Wick_M"] == 12 * n * n
        assert row["inner_Wick_kappa"] == -48 * n
        assert row["kappa2_budget_96n"] == 96 * n
        assert row["C_diag_le_96n"] is True
        assert row["room_positive"] is True
        # C_diag formula
        assert C_diag_formula(p) == Fraction(4 * n * (11 * n - 14), p * p)


def test_C_diag_le_96n_universal():
    for p in _primes(60):
        n = n_of(p)
        assert C_diag_formula(p) <= 96 * n
        # explicit disc: (11n-14) <= 24 p^2
        assert 11 * n - 14 <= 24 * p * p


def test_gap_via_kappa_budget():
    out = prove_gap_via_kappa_budget()
    assert out["proved"] is True
    assert "96n" in out["theorem"] or "κ" in out["theorem"]


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1596.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1596.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.96"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["wick_kappa_calculus"]["proved"] is True

    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3.get("skipped") is not True
    assert c5.get("skipped") is not True
    # identities
    assert c3["ySigy_constant_2n"] is True
    assert c5["ySigy_constant_2n"] is True
    assert c3["inner_match"] is True
    assert c5["inner_match"] is True
    # κ budget
    assert c3["kappa2_le_96n"] is True
    assert c5["kappa2_le_96n"] is True
    assert abs(c3["kappa2"] - c3["kappa2_budget_96n"]) < 1e-4  # eq at p=3
    assert c5["kappa2"] < c5["kappa2_budget_96n"] - 1.0  # strict at p=5
    # mult
    assert c3["mult_eq_d"] is True
    assert c5["mult_eq_d"] is True
    assert c3["gap_by_mult_d"] is False
    assert c5["gap_by_mult_d"] is True

    c7 = data["cert"].get("7")
    if c7 and not c7.get("skipped"):
        assert c7["kappa2_le_96n"] is True
        assert c7["mult_eq_d"] is True
        assert c7["gap_by_mult_d"] is True


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.96" in sol
    idx = sol.index("15.96")
    assert "OPEN" in sol[idx : idx + 4000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.96" in graph
