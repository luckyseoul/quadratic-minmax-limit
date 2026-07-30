"""Tests for Prop 15.95: Wick≤thr algebra; strengthened gap; C_diag; mult certs."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1595 import (  # noqa: E402
    C_diag_formula,
    d_of,
    n_of,
    prove_C_diag_formula,
    prove_strengthened_criterion,
    prove_wick_le_thr,
    thr_gap_M2,
    wick_hi,
)


def _primes_to(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_wick_le_thr_algebra():
    primes = _primes_to(40)
    out = prove_wick_le_thr(primes)
    assert out["proved"] is True
    # p=3 fails, p>=5 holds
    assert out["by_p"]["3"]["Wick_le_thr"] is False
    for p in (5, 7, 11, 13, 17, 19, 23):
        row = out["by_p"][str(p)]
        assert row["Wick_le_thr"] is True
        assert row["match_expected"] is True
        n = n_of(p)
        disc = n * n - 16 * n - 96
        assert disc > 0
        # exact thr - wick
        thr = thr_gap_M2(p)
        w = wick_hi(p)
        assert thr - w == n * disc // 2


def test_wick_le_thr_universal_disc():
    """n=p^2+1 >=26 for p>=5 implies disc n^2-16n-96>=41."""
    for p in _primes_to(60):
        n = n_of(p)
        disc = n * n - 16 * n - 96
        if p >= 5:
            assert disc >= 41
            assert wick_hi(p) <= thr_gap_M2(p)
        else:
            assert disc < 0


def test_strengthened_criterion():
    out = prove_strengthened_criterion()
    assert out["proved"] is True
    assert "Wick_hi" in out["theorem"] or "12n" in out["theorem"]
    assert "mult" in out["theorem"].lower() or "mult" in out["open_inputs"]


def test_C_diag_formula():
    out = prove_C_diag_formula(_primes_to(30))
    assert out["proved_formula"] is True
    for p in (3, 5, 7, 11):
        n = n_of(p)
        assert C_diag_formula(p) == Fraction(4 * n * (11 * n - 14), p * p)
        assert out["by_p"][str(p)]["C_diag"] == str(C_diag_formula(p))


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1595.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1595.py to generate evidence"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.95"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["sixteen_N_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["wick_le_thr"]["proved"] is True
    assert data["strengthened_criterion"]["proved"] is True

    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3.get("skipped") is not True
    assert c5.get("skipped") is not True
    # mult = d
    assert c3["lam2_mult"] == c3["d"] == 5
    assert c5["lam2_mult"] == c5["d"] == 13
    # gap: fail p=3, hold p=5
    assert c3["gap_by_mult_d"] is False
    assert c5["gap_by_mult_d"] is True
    assert c5["sum_M2_le_thr"] is True
    assert c3["sum_M2_le_thr"] is False
    assert c5["sum_M2_le_Wick_hi"] is True
    assert c5["Wick_hi_le_thr"] is True
    # p=7 if present
    c7 = data["cert"].get("7")
    if c7 and not c7.get("skipped"):
        assert c7["lam2_mult"] == c7["d"] == 25
        assert c7["gap_by_mult_d"] is True
        assert c7["sum_M2_le_thr"] is True


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.95" in sol
    idx = sol.index("15.95")
    assert "OPEN" in sol[idx : idx + 5000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "15.95" in hand
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.95" in graph
