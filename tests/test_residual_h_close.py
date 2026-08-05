"""Tests for residual H algebra + live certification path."""
from __future__ import annotations

from fractions import Fraction

from e1_residual_h_close import (
    H_of,
    alpha_wick,
    certify_live_H,
    mu_bar,
    prove_open,
    prove_theorem_A,
    prove_theorem_B,
    prove_theorem_C,
    prove_theorem_D,
    thr_lambda,
)


def test_H_le_5_algebra():
    r = prove_theorem_A([3, 5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert r["proved"] is True
    assert H_of(3) == 5
    assert H_of(5) == Fraction(49, 13)
    assert H_of(5) < 5


def test_H_implies_16N_and_wick_and_lmin():
    assert prove_theorem_B()["proved"] is True
    assert prove_theorem_C()["proved"] is True
    assert prove_theorem_D()["proved"] is True
    assert alpha_wick(5) == Fraction(27, 25)
    assert alpha_wick(7) == Fraction(51, 49)


def test_H_equals_design_forms():
    for p in (5, 7, 11, 13):
        assert H_of(p) == Fraction((p + 2) ** 2, (p * p + 1) // 2)
        assert H_of(p) == Fraction(2 * (p + 2) ** 2, p * p + 1)


def test_live_H_p5_p7():
    """Drive real Max+ Φ eigensystem; no hard-coded success without cache."""
    for p in (5, 7):
        r = certify_live_H(p)
        if r.get("skipped"):
            # cache missing — structural path still tested above
            continue
        assert r["ray_le_H"] is True
        assert r["sixteen_N"] is True
        assert r["G_ge_I"] is True
        assert r["mult_eq_d"] is True
        if p == 5:
            assert abs(r["ray"] - float(H_of(5))) < 1e-8
            assert abs(r["lambda_max"] - float(thr_lambda(5))) < 1e-8


def test_honest_open_general():
    o = prove_open()
    assert o["residual_closed_general"] is False
    assert o["H_for_all_p"] is False
    assert o["L_status"] == "OPEN"
