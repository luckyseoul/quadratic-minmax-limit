"""Tests for Prop 15.103: δ-bound cert p=3,5,7."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_96n, delta2_budget_hyp  # noqa: E402
from e1_gmin_m4_prop15103 import (  # noqa: E402
    prove_budget_formulas,
    prove_reduction,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_budget_formulas():
    out = prove_budget_formulas(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        closed = Fraction(
            4 * (p * p - 9) * (p * p - 1) ** 2,
            3 * (p * p - 5) * (p * p + 1),
        )
        assert delta2_budget_hyp(p) == closed
    assert delta2_budget_hyp(5) == Fraction(1536, 65)
    assert delta2_budget_hyp(7) == Fraction(3072, 55)


def test_reduction():
    red = prove_reduction()
    assert red["proved"] is True
    assert red["general_delta_bound_proved"] is False
    assert red["kappa_bound_proved_for_all_p"] is False


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop15103.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop15103.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.103"
    assert data["L_status"] == "OPEN"
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["general_delta_bound_proved"] is False
    assert data["delta_bound_certified_primes"] == [3, 5, 7]
    assert data["proved"] is True
    for p in ("3", "5", "7"):
        c = data["cert"][p]
        assert c.get("skipped") is not True
        assert c["ok"] is True
        assert c["kappa2_le_96n"] is True
        assert c["kappa2_le_hyp"] is True
        assert c["delta2_le_hyp_budget"] is True
    # equality only p=3,5
    assert data["cert"]["3"]["kappa2_eq_hyp"] is True
    assert data["cert"]["5"]["kappa2_eq_hyp"] is True
    assert data["cert"]["7"]["kappa2_eq_hyp"] is False
    # p=7 ratio ~0.187
    r = data["cert"]["7"]["ratio_delta2_over_hyp"]
    assert 0.15 < r < 0.25
    # 16N
    assert data["sixteen_N"]["5"]["sixteen_N_holds"] is True
    assert data["sixteen_N"]["7"]["sixteen_N_holds"] is True


def test_solution_handoff_graph():
    assert "15.103" in (ROOT / "solution.md").read_text()
    hand_head = (ROOT / "HANDOFF.md").read_text()[:5000]
    assert "**Current mathematical status:** **L OPEN.**" in hand_head
    assert "15.103" in (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
