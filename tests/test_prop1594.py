"""Tests for Prop 15.94: P⊙P annihilates range(P); gap criterion."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1594 import (  # noqa: E402
    d_of,
    n_of,
    prove_annihilator,
    prove_gap_criterion,
)


def test_annihilator_proved():
    out = prove_annihilator()
    assert out["proved"] is True
    assert "range(P)" in out["theorem"] or "range(P)" in out["theorem"].replace(" ", "")


def test_gap_criterion_algebra():
    primes = []
    for p in range(3, 40):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            primes.append(p)
    out = prove_gap_criterion(primes)
    assert out["proved"] is True
    # room for kappa2 at p>=5
    for p in (5, 7, 11, 13):
        row = out["by_p"][str(p)]
        assert row["room_positive"] is True
        d = d_of(p)
        n = n_of(p)
        assert row["thr_sum_M2"] == 4 * d * d * (d + 4)
        assert row["wick_base_12n2_minus_48n"] == 12 * n * n - 48 * n


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1594.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.94"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["sixteen_N_proved_for_all_p"] is False
    assert data["proved"] is True
    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["annihilator_ok"] is True
    assert c5["annihilator_ok"] is True
    assert c3["sum_P3_ok"] is True
    assert c5["sum_P3_ok"] is True
    # p=3: mult=d but gap fails (correct)
    assert c3["lam2_mult"] == c3["d"]
    assert c3["gap_by_mult_d"] is False
    # p=5: mult=d and gap holds via criterion
    assert c5["lam2_mult"] == c5["d"]
    assert c5["mult_ge_d"] is True
    assert c5["gap_by_mult_d"] is True
    assert c5["sixteen_N_by_mult_d"] is False
    assert c5["sqrt_Q_over_d"] <= c5["gap_thr"] + 1e-12
    assert c5["sqrt_Q_over_d"] > c5["thr16"]


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.94" in sol
    idx = sol.index("15.94")
    assert "OPEN" in sol[idx : idx + 4000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
