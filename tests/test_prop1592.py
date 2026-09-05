"""Tests for Prop 15.92: constant pairing sum; H/16N spectral reductions."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1592 import (  # noqa: E402
    H,
    d_of,
    n_of,
    prove_pairing_sum_constant,
    prove_spectral_reductions,
)


def _primes(hi: int = 40) -> list[int]:
    out = []
    for p in range(3, hi):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_pairing_sum_formula_algebra():
    out = prove_pairing_sum_constant(_primes(60))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19):
        n = n_of(p)
        target = Fraction(n * (n - 1) * (n - 2), 8)
        # direct substitute
        f2 = (n * p) ** 2
        gamma = p * p * n
        cF2 = n * (n - 1)
        pred = Fraction(f2, 8) - Fraction(gamma, 2) + Fraction(cF2, 4)
        assert pred == target
        assert out["by_p"][str(p)]["match"] is True


def test_pairing_values():
    assert n_of(3) * 9 * 8 // 8 == 90  # 10*9*8/8=90
    assert n_of(5) * 25 * 24 // 8 == 1950


def test_spectral_reductions():
    out = prove_spectral_reductions(_primes())
    assert out["proved"] is True
    for p in (3, 5, 7):
        assert out["by_p"][str(p)]["H_implies_16N"] is True
        assert H(p) <= 5


def test_H_vs_16N_thresholds():
    # H thr stricter than 16N thr for p>3
    for p in (5, 7, 11):
        d = d_of(p)
        thr16 = Fraction(4, d * d)  # over N
        thrH = (3 + H(p)) / (2 * d * d)
        assert thrH < thr16


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1592.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.92"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["sixteen_N_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["pairing_sum"]["proved"] is True
    assert data["pairing_cert"]["3"]["ok"] is True
    assert data["pairing_cert"]["5"]["ok"] is True
    assert data["W_spectrum_cert"]["3"]["H_equality"] is True
    assert data["W_spectrum_cert"]["5"]["H_equality"] is True
    # p=5 mult pattern
    eigs = data["W_spectrum_cert"]["5"]["eigs"]
    mults = [e["mult"] for e in eigs if e["value"] != 0]
    assert mults[0] == 1  # λ1
    assert mults[1] == 13  # d
    assert 26 in mults


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.92" in sol
    idx = sol.index("15.92")
    assert "OPEN" in sol[idx : idx + 3500]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
