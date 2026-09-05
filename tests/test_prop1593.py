"""Tests for Prop 15.93: Gu/FFT structure; 16N ≡ λ_max(FFT|1⊥)≤8N."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1593 import (  # noqa: E402
    H,
    d_of,
    n_of,
    prove_C_direction,
    prove_equivalences,
    prove_FFT_allones,
)


def _primes(hi: int = 40) -> list[int]:
    out = []
    for p in range(3, hi):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_FFT_allones_algebra():
    out = prove_FFT_allones(_primes())
    assert out["proved"] is True
    # Nd formula: for any p, eigenvalue is N*d
    for p in (3, 5, 7, 11):
        assert out["by_p"][str(p)]["identity_check"] is True


def test_equivalences():
    out = prove_equivalences(_primes())
    assert out["proved"] is True
    assert "8N" in out["theorem"]
    assert "16N" in out["theorem"]
    for p in (5, 7, 11):
        assert out["by_p"][str(p)]["d_ge_8_for_gap_from_16N"] is True
        assert out["by_p"][str(p)]["H_le_5"] is True


def test_C_direction():
    assert prove_C_direction()["proved"] is True


def test_H_vs_8N_threshold():
    # N(3+H) ≤ 8N iff H≤5
    for p in _primes(50):
        assert (3 + H(p)) <= 8 + Fraction(1, 10**9)
        if p > 3:
            assert (3 + H(p)) < 8


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1593.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.93"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["sixteen_N_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["FFT_allones"]["proved"] is True
    c3 = data["Gu_cert"]["3"]
    c5 = data["Gu_cert"]["5"]
    assert c3["FFT_1_equals_dN_1"] is True
    assert c5["FFT_1_equals_dN_1"] is True
    assert c3["holds_16N"] is True
    assert c3["equality_16N"] is True
    assert c5["holds_16N"] is True
    assert c5["equality_H"] is True
    assert c5["lam_max_FFT_1perp"] < c5["thr_8N"]
    # rank formula
    assert c5["rank_Gu"] == c5["rank_formula"]
    assert c3["rank_Gu"] == c3["rank_formula"]


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.93" in sol
    idx = sol.index("15.93")
    assert "OPEN" in sol[idx : idx + 4000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
