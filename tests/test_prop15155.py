"""Tests for Prop 15.155 — Aut-line e4(s), Tχ, Q; η=c1 R4+c0."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15146 import R4_flat, R4_suf  # noqa: E402
from e1_gmin_m4_prop15153 import mu4_flat, mu4_from_R4, R4_from_W  # noqa: E402
from e1_gmin_m4_prop15154 import eta_of_mu4, eta_suf, kappa_main  # noqa: E402
from e1_gmin_m4_prop15155 import (  # noqa: E402
    Q_of_s,
    c0_eta,
    c1_eta,
    e4_binomial,
    e4_of_s,
    eta_from_R4,
    main,
    primes_ge,
    prove_Tchi_Q_structure,
    prove_census_and_crude,
    prove_e4_poly,
    prove_eta_affine,
    prove_open,
)


def test_e4_poly():
    a = prove_e4_poly(primes_ge(5, 31))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        for s in range(-n, n + 1, 2):
            assert e4_of_s(n, s) == e4_binomial(n, s)
        assert e4_of_s(n, 0) == Fraction(n * (n - 2), 8)


def test_eta_affine_closed_forms():
    e = prove_eta_affine(primes_ge(5, 97))
    assert e["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        n = n_of(p)
        x = p * p
        assert c1_eta(p) == Fraction(16, n * (n - 1) * (n - 2) * (n - 3))
        assert c0_eta(p) == Fraction(-(x * x + 4 * x - 9), x * (x - 2))
        assert eta_from_R4(p, R4_flat(p)) == mu4_flat(p) - kappa_main(p)
        assert eta_from_R4(p, R4_suf(p)) == eta_suf(p)
        # residual equivalence
        assert (eta_suf(p) - c0_eta(p)) / c1_eta(p) == R4_suf(p)


def test_census_and_Q():
    d = prove_census_and_crude()
    assert d["proved_census"] is True
    assert d["crude_Es4_dead"] is True
    for p in (5, 7):
        R4 = R4_from_W(p)
        assert eta_from_R4(p, R4) == eta_of_mu4(p, mu4_from_R4(p, R4))
        assert eta_from_R4(p, R4) <= eta_suf(p)
    # Q formula sanity: at s=0
    for p in (5, 7, 11):
        n = n_of(p)
        assert Q_of_s(p, 0) == Fraction(p, 4) * n * (6 - n)


def test_open_and_main():
    b = prove_Tchi_Q_structure()
    assert b["proved"] is True
    g = prove_open()
    assert g["residual_closed_general"] is False
    assert g["aut_line_closed_forms"] is True
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["e4_poly_in_s"] is True
    assert out["proved"]["eta_affine_R4"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15155.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.155"
    assert "c1" in data["closed_forms"]
