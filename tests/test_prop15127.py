"""Tests for Prop 15.127 — closed W_{p+1}; inversive plane; W=d false."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15127 import (  # noqa: E402
    W_hoffman_closed,
    W_hoffman_int,
    chi4,
    enumerate_W_hoffman,
    hoffman_ED4_part,
    is_prime,
    main,
    prove_formula_algebra,
    prove_inversive_plane,
    prove_replication,
    r_hoffman,
)


def test_inversive_plane():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_inversive_plane(primes)
    assert r["proved"] is True
    assert r["by_p"]["5"]["lam3"] == "1"
    assert r["by_p"]["5"]["lam2"] == "6"  # p+1=6


def test_formula_algebra_and_branches():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_formula_algebra(primes)
    assert r["proved"] is True
    assert chi4(5) == 1
    assert chi4(7) == -1
    assert W_hoffman_int(5) == d_of(5) == 13
    assert W_hoffman_int(7) == 11
    assert W_hoffman_int(7) != d_of(7)
    assert W_hoffman_int(11) == 17
    assert W_hoffman_closed(7) == Fraction(3 * 7 + 1, 2)


def test_census_p3_p5_p7():
    # p=3,5,7 sufficient for counterexample + formula; p=11 optional (fast)
    for p in (3, 5, 7):
        c = enumerate_W_hoffman(p, workers=8 if p < 7 else 86)
        assert c["seed_coclique"] is True
        assert c["match_formula"] is True
        assert c["W_enum"] == W_hoffman_int(p)
    assert enumerate_W_hoffman(7)["eq_d"] is False


def test_replication_ed4_main_open():
    primes = [p for p in range(3, 40) if is_prime(p)]
    d = prove_replication(primes)
    assert d["proved"] is True
    assert r_hoffman(5) == Fraction(3)
    assert r_hoffman(7) == Fraction(11 * 8, 50)

    # corrected contrib differs from 2d form at p=7
    N7 = 11452
    assert hoffman_ED4_part(7, N7) != Fraction(
        2 * d_of(7) * (49 - 14 - 1) ** 4, N7
    )

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["inversive_plane"] is True
    assert out["proved"]["W_formula_census_p3_5_7_11"] is True
    assert out["proved"]["W_eq_d_for_all_p"] is False
    assert out["proved"]["W_eq_d_counterexample_p7"] is True
    assert out["proved"]["weight_enumerator_closed"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15127.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.127"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["cert"]["counterexample_p7"]["W"] == 11
