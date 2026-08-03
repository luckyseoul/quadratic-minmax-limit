"""Tests for Prop 15.154 — avg(χ κ) and η residual budget."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15153 import mu4_flat, mu4_from_R4, mu4_suf, R4_from_W  # noqa: E402
from e1_gmin_m4_prop15154 import (  # noqa: E402
    avg_chi_kappa,
    eta_of_mu4,
    eta_suf,
    kappa_main,
    main,
    primes_ge,
    prove_avg_chi_kappa,
    prove_census,
    prove_decomposition,
    prove_open,
    sum_chi_kappa_algebra,
)


def test_avg_chi_kappa_algebra():
    b = prove_avg_chi_kappa(primes_ge(5, 97))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        n = n_of(p)
        assert avg_chi_kappa(p) == Fraction(3, n - 3)
        assert avg_chi_kappa(p) == Fraction(3, p * p - 2)
        assert sum_chi_kappa_algebra(p) / comb(n, 4) == Fraction(3, n - 3)
        assert kappa_main(p) == Fraction(3, p * p * (p * p - 2))
        assert kappa_main(p) == avg_chi_kappa(p) / Fraction(p * p)


def test_eta_budget_decomposition():
    c = prove_decomposition(primes_ge(5, 97))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        x = p * p
        km = kappa_main(p)
        es = eta_suf(p)
        assert km < mu4_flat(p) < mu4_suf(p)
        assert es == mu4_suf(p) - km
        assert es == Fraction(4 * (13 * x + 15), x * x * (x - 2) * (x - 5))
        assert mu4_flat(p) - km == Fraction(32, x * (x - 2) * (x - 5))
        # residual equivalence on synthetic points
        assert km + es == mu4_suf(p)


def test_census_eta():
    e = prove_census()
    assert e["proved_census"] is True
    for p in (5, 7):
        mu = mu4_from_R4(p, R4_from_W(p))
        eta = eta_of_mu4(p, mu)
        assert eta == mu - kappa_main(p)
        assert eta <= eta_suf(p)
        assert mu <= mu4_suf(p)
        assert eta > 0  # residual mass is positive


def test_open_and_main():
    g = prove_open()
    assert g["residual_closed_general"] is False
    assert g["eta_le_eta_suf_all_p"] is False
    assert g["kappa_main_closed"] is True
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["avg_chi_kappa_maxplus_free"] is True
    assert out["proved"]["mu4_kappa_eta_decomposition"] is True
    assert out["proved"]["census_eta_p5_p7"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["status"]["kappa_main_closed"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15154.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.154"
    assert data["closed_forms"]["kappa_main"] == "3/(p^2 (p^2-2))"
