"""Tests for Prop 15.166: 16N ⇔ Q₂ bound; Wick C-eigen."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15166 import (
    certify_q2_census,
    h2_of,
    lambda_max_Phi_from_Q2,
    main,
    prove_theorem_B,
    prove_theorem_C,
    q2_threshold_16N,
)
from e1_gmin_m4_prop15100 import d_of
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7
from e1_gmin_m4_prop15164 import Es4_star
from e1_gmin_m4_prop15165 import ES4_EXACT


def test_h2_formula():
    assert h2_of(13) == 90
    assert h2_of(25) == 324


def test_q2_threshold_positive():
    for p, N in ((5, 260), (7, 11452), (11, 1000)):
        thr = q2_threshold_16N(p, N)
        assert thr > 0
        d = d_of(p)
        assert thr == Fraction(4 * N, d * (d - 1))


def test_Phi_Q2_relation():
    for p, N, spec in ((5, 260, SPECTRUM_P5), (7, 11452, SPECTRUM_P7)):
        lam_phi = max(spec)
        # invert: Q2 from Phi, then back
        d = d_of(p)
        lam_q2 = Fraction(lam_phi) * Fraction(N, 4 * d * (d - 1))
        back = lambda_max_Phi_from_Q2(p, N, lam_q2)
        assert back == lam_phi


def test_wick_is_C_eigen():
    b = prove_theorem_B()
    assert b["proved"]
    for p in ("5", "7"):
        assert b["by_p"][p]["ok"]


def test_equivalence_algebra():
    c = prove_theorem_C()
    assert c["proved"]
    for p, r in c["by_p"].items():
        assert r["match_equiv"]


def test_census_sixteen_N_via_Q2():
    cert = certify_q2_census()
    assert cert["certified"]
    for p in (5, 7):
        r = cert["by_p"][str(p)]
        assert r["Q2_le_threshold"]
        assert r["Es4_le_star"]
        assert ES4_EXACT[p] <= Es4_star(p)


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["sixteen_N_iff_Q2_bound"] is True
    assert out["proved"]["Wick_is_C_eigen"] is True
    assert out["proved"]["Q2_bound_for_all_p"] is False
    assert out["proved"]["residual_for_all_p_ge_5"] is False
    assert out["proved"]["census_Q2_p5_p7"] is True
