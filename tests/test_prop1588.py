"""Tests for Prop 15.88: pairwise sum, H-gap algebra, g via S3, p=5 cert."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1588 import (  # noqa: E402
    H,
    gap_diff,
    gap_diff_closed_num,
    n_half,
    prove_H_gap_algebra,
    prove_g_via_S3,
    prove_pairwise_sum_identity,
    prove_settlement_under_H,
    three_plus_H,
)


def test_pairwise_sum_identity():
    p = prove_pairwise_sum_identity()
    assert p["proved"] is True
    assert "∑_{i<j}" in p["statement"] or "sum" in p["statement"].lower()


def test_H_gap_algebra_closed_form():
    # explicit values
    assert H(5) == Fraction(49, 13)
    assert three_plus_H(5) == Fraction(88, 13)
    assert n_half(5) == Fraction(13, 1)
    assert gap_diff(5) == Fraction(81, 13)
    assert gap_diff_closed_num(5) == 324
    # p=3 gap fails
    assert gap_diff(3) < 0
    assert three_plus_H(3) > n_half(3)
    # primes p≥5
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31]
    for p in primes:
        assert gap_diff(p) > 0
        assert three_plus_H(p) < n_half(p)
        assert gap_diff(p) == Fraction(gap_diff_closed_num(p), 2 * (p * p + 1))
    out = prove_H_gap_algebra([3] + primes)
    assert out["proved"] is True


def test_settlement_and_g_via_S3():
    s = prove_settlement_under_H()
    assert s["proved_form"] is True
    assert "bi-tight" in s["chain"].lower() or "empty" in s["chain"].lower()
    g = prove_g_via_S3()
    assert g["proved"] is True
    assert "star·S3" in g["GD_iff"] or "S3" in g["GD_iff"]


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1588.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.88"
    assert data["L_status"] == "OPEN"
    assert data["proved"] is True
    assert data["H_gap_algebra"]["proved"] is True
    assert data["pairwise_sum"]["proved"] is True
    p5 = data.get("p5_g_tau")
    if p5 and not p5.get("skipped"):
        assert p5["all_g_le_0"] is True
        assert p5["g_function_of_tau1_only"] is True
        assert p5["matches_exact_fractions"] is True


def test_solution_handoff_open():
    sol = (ROOT / "solution.md").read_text()
    assert "15.88" in sol
    idx = sol.index("15.88")
    assert "OPEN" in sol[idx : idx + 3500]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
