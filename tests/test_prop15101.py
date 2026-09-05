"""Tests for Prop 15.101: Fickus Gram residual / bulk-variance orth reduction."""
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
    room_orth,
    wick_hi,
    wick_lo,
)
from e1_gmin_m4_prop15101 import (  # noqa: E402
    eps_hyp,
    eps_wick,
    m_bulk,
    max_bulk_levels,
    prove_eps_criterion,
    prove_fickus_rank,
    prove_orth_variance_identity,
    prove_reduction,
    rank_popp,
)


def _primes(m: int = 40) -> list[int]:
    out = []
    for p in range(3, m):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_fickus_rank_structure():
    out = prove_fickus_rank(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17):
        d = d_of(p)
        assert rank_popp(p) == (d - 1) * (d - 2) // 2
        assert rank_popp(p) == 1 + m_bulk(p)
        assert rank_popp(p) < d * d
        assert max_bulk_levels(p) == m_bulk(p) // (d - 1)
        assert max_bulk_levels(p) >= 1


def test_orth_variance_identity():
    out = prove_orth_variance_identity(_primes(40))
    assert out["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        assert ED4_flat(p) - wick_lo(p) == proj_kappa2(p)
        assert wick_hi(p) - ED4_flat(p) == room_orth(p)
        # ε from room formula
        d = d_of(p)
        room = room_orth(p)
        eps_from_room = room * Fraction(d - 3, 32 * d * (d - 1) ** 2)
        assert eps_wick(p) == eps_from_room
        assert eps_wick(p) == Fraction(4 * (p * p - 9), (p * p - 1) ** 2)


def test_eps_criterion():
    out = prove_eps_criterion(_primes(40))
    assert out["proved"] is True
    assert eps_wick(3) == 0
    assert eps_hyp(3) == 0
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert eps_wick(p) > 0
        assert eps_hyp(p) <= eps_wick(p)
        assert eps_hyp(p) > 0
    # p=5 exact: ε=4*16/24²=64/576=1/9
    assert eps_wick(5) == Fraction(1, 9)
    d5 = d_of(5)
    assert eps_hyp(5) == Fraction(1, 9) * Fraction((d5 - 1) ** 2, d5 * d5)


def test_reduction_statement():
    red = prove_reduction()
    assert red["proved"] is True
    assert red["general_orth_inequality_proved"] is False
    assert red["kappa_bound_proved_for_all_p"] is False


def test_room_hyp_link():
    for p in (3, 5, 7, 11, 13):
        d = d_of(p)
        n = n_of(p)
        room_h = room_orth(p) * Fraction((d - 1) ** 2, d * d)
        assert proj_kappa2(p) + room_h == kappa_hyp(p)
        assert kappa_hyp(p) <= 96 * n


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop15101.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop15101.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.101"
    assert data["L_status"] == "OPEN"
    assert data["kappa_bound_proved_for_all_p"] is False
    assert data["general_orth_inequality_proved"] is False
    assert data["fickus_rank"]["proved"] is True
    assert data["orth_variance_identity"]["proved"] is True
    assert data["eps_criterion"]["proved"] is True
    assert data["proved"] is True
    c3 = data["cert"]["3"]
    c5 = data["cert"]["5"]
    assert c3["ok"] is True
    assert c5["ok"] is True
    assert c3["kappa2_eq_hyp"] is True
    assert c5["kappa2_eq_hyp"] is True
    assert c5["orth_eq_room_hyp"] is True
    assert c5["rank_match"] is True
    assert c5["mult_ge_d_minus_1"] is True
    # p=5: sufficient λ2 test does NOT fire (not necessary)
    assert c5["lambda2_le_thr_wick"] is False
    # but 16N holds at p=5
    assert c5["lambda2_le_4_over_N"] is True


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.101" in sol
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.101" in graph
