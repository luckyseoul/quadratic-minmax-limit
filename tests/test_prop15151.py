"""Tests for Prop 15.151 — m_e covariance formula and t-identities."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15149 import k_act_from_delta2, k_suf  # noqa: E402
from e1_gmin_m4_prop15150 import mixture_tot, pi_e, triple_counts  # noqa: E402
from e1_gmin_m4_prop15151 import (  # noqa: E402
    E_te,
    M_E_P5,
    P5_DELTA2,
    main,
    primes_ge,
    prove_E_te,
    prove_constancy,
    prove_cov_formula,
    prove_cs_weak,
    prove_identities,
    prove_m_e_p5,
    prove_open,
    regular_set_identities,
)


def test_regular_set_identities():
    a = prove_identities()
    assert a["proved"] is True
    # spot-check a few k
    assert regular_set_identities(6, 5, (20, 0, 0, 0))["ok"]
    assert regular_set_identities(14, 5, (108, 180, 72, 4))["ok"]
    assert regular_set_identities(26, 5, (520, 1170, 780, 130))["ok"]


def test_E_te_and_constancy():
    b = prove_constancy()
    assert b["proved_census"] is True
    c = prove_E_te(primes_ge(5, 31))
    assert c["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        s = sum(E_te(p, e) for e in range(4))
        n = n_of(p)
        assert s == Fraction(n * (n * n - 4), 48)
        for e in range(4):
            tc = triple_counts(p)
            assert E_te(p, e) == Fraction(tc[e] * (p + 3 - 2 * e), 8 * p)


def test_cov_formula_and_m_e_p5():
    d = prove_cov_formula()
    assert d["proved"] is True
    assert d["p5_Cov_le_budget"] is True
    e = prove_m_e_p5()
    assert e["proved"] is True
    pi = pi_e(5)
    Em = sum(pi[e] * M_E_P5[e] for e in range(4))
    assert Em == k_act_from_delta2(5, P5_DELTA2)
    assert Em <= k_suf(5)
    assert M_E_P5[0] < k_suf(5) < M_E_P5[3]
    # reconstruct m from Cov identity
    p, n = 5, n_of(5)
    tc = triple_counts(5)
    for e in range(4):
        den = tc[e] * (p + 3 - 2 * e)
        Cov = (M_E_P5[e] - Fraction(n, 2)) * den / (8 * p)
        m_rec = Fraction(n, 2) + 8 * p * Cov / den
        assert m_rec == M_E_P5[e]


def test_cs_weak_and_main():
    f = prove_cs_weak([5, 7, 11])
    assert f["proved_too_weak"] is True
    g = prove_open()
    assert g["residual_closed_general"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["m_e_covariance_formula"] is True
    assert out["proved"]["E_te_maxplus_free"] is True
    assert out["proved"]["m_e_bound_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15151.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.151"
    assert "m_e" in data["closed_forms"]
