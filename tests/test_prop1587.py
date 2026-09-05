"""Tests for Prop 15.87: K4 star theorem, S1 pattern, CS weakness, E[U1²]."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1587 import (  # noqa: E402
    d1_one,
    prove_CS_weakness,
    prove_S1_pattern_and_GD_reformulation,
    prove_k4_star_theorem,
)


def test_k4_star_theorem():
    k4 = prove_k4_star_theorem()
    assert k4["proved"] is True
    assert k4["n_labelings"] == 64
    assert k4["n_kappa_abs_1"] == 48
    assert k4["n_sum_star_0"] == 48
    assert k4["n_prod_star_1"] == 48
    assert k4["n_exactly_two_star_plus"] == 48
    assert k4["n_sigma_sum_eq_4kappa"] == 64


def test_S1_pattern_and_GD_iff():
    pat = prove_S1_pattern_and_GD_reformulation()
    assert pat["proved"] is True
    assert "g≤0" in pat["GD_iff"] or "g <= 0" in pat["GD_iff"]
    assert "S1(a)=g" in pat["S1_pattern"] or "g · star" in pat["S1_pattern"]


def test_CS_too_weak():
    primes = [5, 7, 11, 13, 17, 19, 23]
    cs = prove_CS_weakness(primes)
    assert cs["proved_CS_too_weak_for_p_ge_5"] is True
    for p in primes:
        assert cs["by_p"][str(p)]["CS_beats_Wick_scale"] is True
        assert d1_one(p) == (3 * p * p - 7) // 4


def test_evidence_json_and_open():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1587.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.87"
    assert data["L_status"] == "OPEN"
    assert data["proved"] is True
    assert data["k4_star"]["proved"] is True
    assert data["EU2_census"]["near_d1_all"] is True
    for p in ("3", "5", "7"):
        row = data["EU2_census"]["by_p"][p]
        assert row["EU2_near_d1"] is True
        assert row["CS_too_weak"] is True


def test_solution_handoff_prop1587_open():
    sol = (ROOT / "solution.md").read_text()
    assert "15.87" in sol
    idx = sol.index("15.87")
    assert "OPEN" in sol[idx : idx + 3000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
