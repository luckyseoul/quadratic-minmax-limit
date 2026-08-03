"""Tests for Prop 15.153 — switched μ₄ residual dictionary."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15146 import R4_flat, R4_suf  # noqa: E402
from e1_gmin_m4_prop15149 import k_act_from_R4, k_suf  # noqa: E402
from e1_gmin_m4_prop15150 import pi_e  # noqa: E402
from e1_gmin_m4_prop15151 import M_E_P5  # noqa: E402
from e1_gmin_m4_prop15153 import (  # noqa: E402
    M_E_P7,
    R4_from_W,
    R4_from_mu4,
    baseline_B,
    d4_from_mu4,
    e2_of,
    main,
    mu4_flat,
    mu4_from_R4,
    mu4_suf,
    primes_ge,
    prove_census,
    prove_low_moments,
    prove_mu4_dictionary,
    prove_open,
)


def test_low_moments_and_e2():
    a = prove_low_moments(primes_ge(5, 37))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert e2_of(p) == Fraction(1, n_of(p) - 1)


def test_mu4_closed_forms_match_R4():
    b = prove_mu4_dictionary(primes_ge(5, 37))
    assert b["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        mf, ms = mu4_flat(p), mu4_suf(p)
        assert R4_from_mu4(p, mf) == R4_flat(p)
        assert R4_from_mu4(p, ms) == R4_suf(p)
        assert mu4_from_R4(p, R4_flat(p)) == mf
        assert mu4_from_R4(p, R4_suf(p)) == ms
        assert mf < ms
        assert mf > 0
        # explicit formulas
        x = p * p
        assert mf == Fraction(3 * x + 17, x * (x - 2) * (x - 5))
        assert ms == Fraction(3 * x * x + 37 * x + 60, x * x * (x - 2) * (x - 5))
        n = n_of(p)
        n4 = n * (n - 1) * (n - 2) * (n - 3)
        assert baseline_B(p) == Fraction(n * (n - 2) * (n - 3) * (n + 5), 16)
        assert R4_flat(p) - baseline_B(p) == Fraction(n4, 16) * mf
        assert d4_from_mu4(p, mf) == R4_flat(p) / Fraction(n4)


def test_census_mu4_and_m_e():
    d = prove_census()
    assert d["proved_census"] is True
    for p in (5, 7):
        R4 = R4_from_W(p)
        mu = mu4_from_R4(p, R4)
        assert mu <= mu4_suf(p)
        assert k_act_from_R4(p, R4) <= k_suf(p)
    # m_e
    Em5 = sum(pi_e(5)[e] * M_E_P5[e] for e in range(4))
    Em7 = sum(pi_e(7)[e] * M_E_P7[e] for e in range(4))
    assert Em5 == k_act_from_R4(5, R4_from_W(5))
    assert Em7 == k_act_from_R4(7, R4_from_W(7))
    assert M_E_P7[0] == Fraction(1397486, 51125)
    assert M_E_P7[3] == Fraction(298724, 10225)


def test_open_and_main():
    g = prove_open()
    assert g["residual_closed_general"] is False
    assert g["mu4_le_mu4_suf_all_p"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["mu4_dictionary_closed"] is True
    assert out["proved"]["switched_low_moments"] is True
    assert out["proved"]["census_mu4_and_m_e"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15153.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.153"
    assert "mu4_flat" in data["closed_forms"]
