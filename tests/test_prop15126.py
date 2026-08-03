"""Tests for Prop 15.126 — geometric Hoffman seed; 1-design; simplex bound."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15126 import (  # noqa: E402
    certify_hoffman_layer,
    hoffman_ED4_contribution,
    is_prime,
    main,
    prove_geometric_seed,
    prove_hoffman_ed4_contrib,
    prove_one_design_algebra,
    prove_simplex_bound,
    subfield_line,
)


def test_geometric_seed():
    r = prove_geometric_seed([3, 5, 7])
    assert r["proved"] is True
    for p in (3, 5, 7):
        assert len(subfield_line(p)) == p + 1
        assert r["by_p"][str(p)]["regular_after_hs_switch"] is True
        assert r["by_p"][str(p)]["alpha_is_0"] is True


def test_one_design_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_one_design_algebra(primes)
    assert r["proved"] is True
    for p in (3, 5, 7, 11):
        assert Fraction(n_of(p) * Fraction(p + 1, 2), p + 1) == d_of(p)


def test_cert_hoffman_p3_p5():
    for p in (3, 5):
        c = certify_hoffman_layer(p)
        assert c.get("skipped") is not True
        assert c["W_eq_d"] is True
        assert c["W_hoffman"] == d_of(p)
        assert c["r_constant"] is True
        assert c["r_eq_half_p_plus_1"] is True
        assert c["rank_eq_d"] is True
        assert c["basis_of_Vplus"] is True
        assert c["match"] is True
    # p=3: constant intersection
    c3 = certify_hoffman_layer(3)
    assert c3["constant_intersection"] is True


def test_simplex_ed4_and_main_open():
    d = prove_simplex_bound([p for p in range(3, 40) if is_prime(p)])
    assert d["proved"] is True
    assert d["p3_equality_case"] is True

    e = prove_hoffman_ed4_contrib([3, 5, 7])
    assert e["proved"] is True
    # formula check
    assert hoffman_ED4_contribution(5, 260) == Fraction(
        2 * 13 * (25 - 10 - 1) ** 4, 260
    )

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["geometric_hoffman_seed"] is True
    assert out["proved"]["one_design_algebra"] is True
    assert out["proved"]["hoffman_eq_d_cert_p3"] is True
    assert out["proved"]["hoffman_eq_d_cert_p5"] is True
    assert out["proved"]["W_hoffman_eq_d_for_all_p"] is False
    assert out["proved"]["weight_enumerator_closed"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15126.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.126"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
