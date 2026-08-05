"""Tests for Prop 15.164: 16N via mult≥d−1 + Es4≤Es4_*."""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15164 import (
    Es4_star,
    bulk_b_star,
    certify_census,
    eta_star,
    kappa4_star,
    main,
    prove_theorem_B,
    sixteen_N_from_dminus1_Es4,
    R_star,
)
from e1_gmin_m4_prop15159 import d_of, SPECTRUM_P5, SPECTRUM_P7
from e1_gmin_m4_prop15162 import Es4_from_W_census


def test_b_star_ge_6_all_primes():
    b = prove_theorem_B()
    assert b["proved"]
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43, 47):
        assert bulk_b_star(p) >= 6
        assert bulk_b_star(p) <= 16
        assert eta_star(p) > 0


def test_Es4_star_consistent():
    for p in (5, 7, 11, 13):
        es = Es4_star(p)
        assert es == kappa4_star(p) + 12 * (p * p + 1) ** 2
        assert R_star(p) == es - (
            __import__(
                "e1_gmin_m4_prop15162", fromlist=["C0_of"]
            ).C0_of(p)
        )


def test_census_sixteen_N():
    c = certify_census()
    assert c["certified"]
    # Global Es4 from Φ spectrum / Gram (15.165); not single-root W at p=7.
    ES4_GLOBAL = {
        5: Fraction(120400, 13),
        7: Fraction(5218435600, 167281),
    }
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        Es4 = ES4_GLOBAL[p]
        mult = spec[max(spec)]
        pred = sixteen_N_from_dminus1_Es4(mult, d_of(p), Es4, p)
        assert pred["sixteen_N"] is True
        assert pred["Es4_le_star"] is True
        assert pred["mult_ge_d_minus_1"] is True
        assert c["by_p"][str(p)]["Es4"] == str(Es4)


def test_predicate_fails_if_Es4_too_large():
    # force Es4 above star
    bad = sixteen_N_from_dminus1_Es4(d_of(7) - 1, d_of(7), Es4_star(7) + 1, 7)
    assert bad["sixteen_N"] is False
    # mult too small
    bad2 = sixteen_N_from_dminus1_Es4(0, d_of(7), 0, 7)
    assert bad2["sixteen_N"] is False


def test_p5_p7_exact_fractions():
    assert Es4_star(5) == Fraction(492752, 53)
    assert bulk_b_star(5) == Fraction(432, 53)
    assert Es4_star(7) == Fraction(8116400, 251)
    assert bulk_b_star(7) == Fraction(2016, 251)


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["mult_ge_d_minus_1"] is True
    assert out["proved"]["Es4_star_for_all_p"] is False
    assert out["proved"]["sixteen_N_for_all_p"] is False
    assert out["proved"]["census_p5_p7_sixteen_N"] is True
