"""Prop 15.196 — spectral Q bound path; residual (i) still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15196 import (
    Q_bound_proved_general,
    Q_ub_from_lambda2,
    certified_paley_Q_and_lambda2,
    hinge_status_196,
    lambda2_for_Q_ub_le_10,
    main,
    residual_i_closed_via_spectral_Q,
    theorem_Q_le_10_blocks_p_ge_5,
    theorem_entry_floor_minus_2,
    theorem_lambda2_threshold_for_Q10,
    theorem_spectral_Q_ub,
    worst_mass_min_under_Q,
)


def test_proved_algebra():
    assert theorem_entry_floor_minus_2()["proved"] is True
    assert theorem_spectral_Q_ub()["proved"] is True
    assert theorem_lambda2_threshold_for_Q10()["proved"] is True


def test_Q_ub_formula_at_threshold():
    for p in (5, 7, 11, 13, 17, 19, 23):
        thr = lambda2_for_Q_ub_le_10(p)
        assert Q_ub_from_lambda2(p, thr) == 10
        # thr = 6 + 5/(n-2)
        n = n_of(p)
        assert thr == 6 + Fraction(5, n - 2)


def test_Q10_blocks_small_primes():
    """Drive real worst_mass_min_under_Q — Q≤10 blocks dual-eq at p=5,7."""
    for p in (5, 7, 11, 13):
        r = worst_mass_min_under_Q(p, Fraction(10))
        assert r["blocks_dual_eq"] is True
    # finite-range theorem
    th = theorem_Q_le_10_blocks_p_ge_5()
    assert th["proved_for_primes_5_to_47"] is True
    assert th["proved_general_all_primes"] is False


def test_census_paley():
    cert = certified_paley_Q_and_lambda2()
    assert Fraction(cert["5"]["Q"]) == Fraction(34512, 4225)
    assert Fraction(cert["5"]["lambda2"]) == Fraction(88, 13)
    assert cert["5"]["blocks_under_Q10"] is True
    assert cert["5"]["blocks_under_Q_ub"] is False  # spectral ub too weak at p=5
    assert Fraction(cert["7"]["Q"]) == Fraction(55009, 8220)
    assert cert["7"]["blocks_under_Q10"] is True
    # p=5 actual Q < 10
    assert Fraction(cert["5"]["Q"]) < 10
    assert Fraction(cert["7"]["Q"]) < 10


def test_p3_spectral_tight():
    cert = certified_paley_Q_and_lambda2()
    assert Fraction(cert["3"]["Q"]) == Fraction(cert["3"]["Q_ub_lambda2"])


def test_predicates_open():
    assert Q_bound_proved_general() is False
    assert residual_i_closed_via_spectral_Q() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_196()
    assert h["entry_floor_minus_2"] is True
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
