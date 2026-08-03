"""Tests for Prop 15.158 — Q₄ closed; Max+ scheme structure; pole bound."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15128 import KNOWN_N  # noqa: E402
from e1_gmin_m4_prop15157 import mu_G4_suf  # noqa: E402
from e1_gmin_m4_prop15158 import (  # noqa: E402
    N_star,
    Q2_of,
    Q4_from_expansion,
    Q4_of,
    main,
    primes_ge,
    prove_Q4_closed,
    prove_chebyshev_weak,
    prove_decomposition_and_conditional,
    prove_homogeneous_and_spectrum,
    prove_not_association_scheme,
    prove_open,
)


def test_Q4_closed_form():
    a = prove_Q4_closed(list(range(3, 40)) + [61, 85, 145])
    assert a["proved"] is True
    for d in (13, 25, 61, 85, 145):
        assert Q4_of(d, Fraction(0)) == Fraction(3, d * d - 1)
        assert Q4_of(d, Fraction(1)) == 1
        for ti in range(0, 11):
            t = Fraction(ti, 10)
            assert Q4_of(d, t) == Q4_from_expansion(d, t)
        # explicit formula
        t = Fraction(1, 3)
        num = (
            Fraction(d + 2) * (d + 4) * t**4
            - 6 * Fraction(d + 2) * t**2
            + 3
        )
        assert Q4_of(d, t) == num / Fraction(d * d - 1)


def test_homogeneous_spectrum_nonscheme():
    b = prove_homogeneous_and_spectrum()
    assert b["proved_census"] is True
    assert b["p5"]["one_homogeneous"] is True
    assert b["p5"]["lambda_2N_mult"] == 13
    assert b["p5"]["lambda_0_mult"] == 247
    c = prove_not_association_scheme()
    assert c["proved_not_scheme"] is True
    # at least one class non-constant
    assert any(not v["constant"] for v in c["by_s"].values())


def test_decomposition_N_star():
    d = prove_decomposition_and_conditional()
    assert d["proved"] is True
    assert N_star(5) == 255
    assert KNOWN_N[5] >= N_star(5)
    assert KNOWN_N[7] >= N_star(7)
    # 2/N <= suf at certified primes
    assert Fraction(2, KNOWN_N[5]) <= mu_G4_suf(5)
    assert Fraction(2, KNOWN_N[7]) <= mu_G4_suf(7)
    # support hypothesis fails (documented)
    assert d["by_p"]["5"]["support_only_nonpos_Q4_except_poles"] is False
    assert d["by_p"]["7"]["support_only_nonpos_Q4_except_poles"] is False
    assert d["by_p"]["5"]["pole_bound_implies_residual"] is False


def test_chebyshev_open_main():
    f = prove_chebyshev_weak(primes_ge(5, 23))
    assert f["proved_too_weak"] is True
    g = prove_open()
    assert g["residual_closed_general"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["Q4_closed_form"] is True
    assert out["proved"]["one_homogeneous_gram_spectrum"] is True
    assert out["proved"]["not_IP_association_scheme"] is True
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["status"]["IP_scheme_blocked"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15158.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.158"
    assert "Q4" in data["closed_forms"]
