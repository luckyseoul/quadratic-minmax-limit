"""Tests for Prop 15.147 — inclusion-density residual / ULC near-miss."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15124 import exact_le3  # noqa: E402
from e1_gmin_m4_prop15128 import W_CENSUS, moments_from_W  # noqa: E402
from e1_gmin_m4_prop15147 import (  # noqa: E402
    P5_DELTA2,
    P7_DELTA2,
    P_gap_poly,
    d1,
    d2,
    d3,
    d4_bud,
    d4_flat,
    d4_of_R4,
    d4_suf,
    main,
    primes_ge,
    prove_d4_residual,
    prove_falling_identity,
    prove_inclusion_le3,
    prove_open,
    prove_ulc_below_suf,
    prove_ulc_fails_census,
    ulc_U,
)
from e1_gmin_m4_prop15146 import R4_flat  # noqa: E402


def test_falling_and_inclusion():
    a = prove_falling_identity(primes_ge(5, 31))
    assert a["proved"] is True
    b = prove_inclusion_le3(primes_ge(5, 31))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert d1(p) == Fraction(1, 2)
        assert d2(p) == Fraction(p * p + 1, 4 * p * p)
        assert d3(p) == Fraction(p * p + 3, 8 * p * p)
        assert ulc_U(p) == d3(p) ** 2 / d2(p)
        assert ulc_U(p) == Fraction((p * p + 3) ** 2, 16 * p * p * (p * p + 1))


def test_d4_residual_and_ulc_gap():
    c = prove_d4_residual(primes_ge(5, 31))
    assert c["proved"] is True
    # p=5 saturates bud; p=7 strict
    d4_5 = d4_of_R4(5, R4_flat(5) + Fraction(3, 2) * P5_DELTA2)
    assert d4_5 == d4_bud(5)
    assert d4_5 <= d4_suf(5)
    d4_7 = d4_of_R4(7, R4_flat(7) + Fraction(3, 2) * P7_DELTA2)
    assert d4_7 <= d4_suf(7)
    d = prove_ulc_below_suf(primes_ge(5, 97))
    assert d["proved"] is True
    for p in primes_ge(5, 47):
        assert ulc_U(p) < d4_suf(p)
        assert P_gap_poly(p * p) > 0
    # factor check at a few points
    for x in (25, 49, 121, 169):
        P = P_gap_poly(x)
        Q = x**4 - 7 * x**3 + 71 * x**2 + 67 * x + 60
        assert P == (x - 1) * Q


def test_ulc_fails_but_residual_holds():
    e = prove_ulc_fails_census()
    assert e["proved"] is True
    for p in (5, 7):
        R4 = moments_from_W(W_CENSUS[p], 4) - exact_le3(p)["exact_le3"]
        d4 = d4_of_R4(p, R4)
        assert d4 > ulc_U(p)
        assert d4 <= d4_suf(p)


def test_open_main():
    f = prove_open()
    assert f["residual_closed_general"] is False
    assert f["ulc_closes_if_true"] is True
    assert f["ulc_holds"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["ULC_U_strictly_below_d4_suf"] is True
    assert out["proved"]["ULC_fails_census_p5_p7"] is True
    assert out["proved"]["type_free_d4_le_d4_suf_all_p"] is False
    assert out["status"]["ULC_would_suffice_if_true"] is True
    assert out["status"]["ULC_viable_as_proof"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15147.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.147"
