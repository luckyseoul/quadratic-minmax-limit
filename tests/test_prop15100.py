"""Tests for Prop 15.100: dual-frame proj; flat≤Wick; κ_hyp."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import (  # noqa: E402
    ED4_flat,
    d_of,
    kappa_hyp,
    n_of,
    proj_kappa2,
    prove_dual_frame_projection,
    prove_flat_le_wick,
    prove_kappa_hyp_le_96n,
    room_orth,
    wick_hi,
    wick_lo,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_dual_frame_projection():
    out = prove_dual_frame_projection(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        row = out["by_p"][str(p)]
        assert row["proj_match"] is True
        assert row["eq_96n"] is True
        n = n_of(p)
        assert proj_kappa2(p) + room_orth(p) == 96 * n


def test_flat_le_wick():
    out = prove_flat_le_wick(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        row = out["by_p"][str(p)]
        assert row["ED4_flat_le_Wick"] is True
        assert row["eq_proj"] is True
        d = d_of(p)
        room_c = Fraction(64 * d * (d - 5), d - 3)
        assert ED4_flat(p) <= wick_hi(p)
        assert wick_hi(p) - ED4_flat(p) == room_c
        assert ED4_flat(p) - wick_lo(p) == proj_kappa2(p)
    # equality only p=3
    assert wick_hi(3) - ED4_flat(3) == 0
    assert wick_hi(5) - ED4_flat(5) > 0


def test_kappa_hyp_le_96n():
    out = prove_kappa_hyp_le_96n(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        kh = kappa_hyp(p)
        budget = 96 * n_of(p)
        assert kh <= budget
        slack = budget - kh
        slack_c = Fraction(
            128 * p * p * (p - 3) * (p + 3),
            (p * p - 5) * (p * p + 1),
        )
        assert slack == slack_c
    assert kappa_hyp(3) == 960
    assert 96 * n_of(3) - kappa_hyp(3) == 0


def test_proj_room_closed_forms():
    # p=3
    assert proj_kappa2(3) == 960
    assert room_orth(3) == 0
    # p=5
    assert proj_kappa2(5) == Fraction(9152, 5)
    assert room_orth(5) == Fraction(3328, 5)
    assert kappa_hyp(5) == Fraction(31168, 13)


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop15100.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop15100.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.100"
    assert data["L_status"] == "OPEN"
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["dual_frame_projection"]["proved"] is True
    assert data["flat_le_wick"]["proved"] is True
    assert data["kappa_hyp_le_96n"]["proved"] is True
    assert data["proved"] is True
    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["ok"] is True
    assert c5["ok"] is True
    assert c3["kappa2_eq_hyp"] is True
    assert c5["kappa2_eq_hyp"] is True
    assert abs(c3["kappa2"] - 960) < 1e-4
    assert c5["kappa2"] < 96 * 26 - 1


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.100" in sol
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "15.100" in hand
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.100" in graph
