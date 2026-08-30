"""Tests for corrected Prop 15.167: valid algebra, retracted implication."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15100 import d_of, n_of
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, dim_Z
from e1_gmin_m4_prop15167 import (
    L_star_bulk6,
    L_star_closed,
    bitight_empty_for_all_primes_ge_5,
    bitight_from_majorization,
    lambda_cycle_ub_from_L,
    main,
    prove_theorem_A,
    prove_theorem_B,
    prove_theorem_C,
    two_d_minus_L_star,
)


def test_L_star_matches_bulk_formula():
    """Shipped closed form equals majorization (T−6(m−k))/k."""
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert L_star_closed(p) == L_star_bulk6(p)
        n, d, m = n_of(p), d_of(p), dim_Z(p)
        k = d - 1
        manual = Fraction(n * (n - 2) - 6 * (m - k), k)
        assert L_star_closed(p) == manual


def test_L_star_lt_2d_algebra():
    """2d − L_* = (p⁴−24p²−1)/(2(p²−1)) > 0 for primes p≥5."""
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 47, 97):
        gap = two_d_minus_L_star(p)
        assert gap > 0
        assert gap == Fraction(p**4 - 24 * p * p - 1, 2 * (p * p - 1))
        assert L_star_closed(p) < 2 * d_of(p)


def test_majorization_predicate_records_retracted_arrow():
    for p in (5, 7, 11, 13, 17):
        r = bitight_from_majorization(p)
        assert r["bitight_empty"] is False
        assert r["retracted"] is True
        assert r["spectral_gap_conditional_on_floor"] is True
        assert r["L_star_lt_2d"] is True
        assert r["lambda_cycle_ub_lt_d"] is True
        # cycle UB strictly below d from L_*/2
        assert lambda_cycle_ub_from_L(L_star_closed(p)) < d_of(p)


def test_p3_outside_claim():
    """p=3 is outside p≥5; poly p⁴−24p²−1 = 81−216−1 < 0 so L_* > 2d."""
    assert two_d_minus_L_star(3) < 0
    r = bitight_from_majorization(3)
    assert r["prime_ge_5"] is False
    assert r["bitight_empty"] is False


def test_all_primes_ge_5_old_bitight_claim_is_false():
    bt = bitight_empty_for_all_primes_ge_5()
    assert bt["proved"] is False
    assert bt["all_bitight_empty"] is False
    assert bt["n_checked"] >= 40


def test_theorems_proved():
    assert prove_theorem_A()["proved"]
    assert prove_theorem_B()["proved"]
    assert prove_theorem_C()["proved"] is False
    assert prove_theorem_C()["retracted"] is True


def test_census_actual_spectrum_below_L_star_but_no_cover_conclusion():
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        lam = max(spec)
        assert lam <= L_star_closed(p)
        assert Fraction(lam) / 2 < d_of(p)
        assert min(spec) >= 6
        assert spec[lam] >= d_of(p) - 1


def test_main_retracts_bi_tight_and_keeps_L_open():
    out = main()
    assert out["proved"]["bitight_empty_for_all_p_ge_5"] is False
    assert out["proved"]["residual_closed_general"] is False
    assert out["proved"]["sixteen_N_for_all_p"] is False
    assert out["proved"]["E1_closed"] is False
    assert out["L_status"] == "OPEN"
    assert out["proved"]["census_spectral_conditions_p5_p7"] is True
