"""Prop 15.198 — Wick_hi ≤ K₄ thr; partition; residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15197 import K4_threshold_for_Q10, Q_from_K4
from e1_gmin_m4_prop15198 import (
    C_diag,
    K4_bound_proved_general,
    Q_from_wick_hi,
    certify_K4_exact,
    eta_sumsq_budget_for_wick_hi,
    hinge_status_198,
    main,
    residual_i_closed_via_K4_wick,
    theorem_K4_partition_identity,
    theorem_Q_wick_lt_10,
    theorem_wick_hi_le_K4_thr,
    wick_hi,
    wick_lo,
)


def test_wick_hi_le_threshold_algebra():
    th = theorem_wick_hi_le_K4_thr()
    assert th["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        w = wick_hi(p)
        thr = int(K4_threshold_for_Q10(p))
        assert w == 12 * n * n + 48 * n
        assert thr == 16 * n * n - 14 * n
        assert w <= thr
        assert thr - w == 2 * n * (2 * n - 31)
    # p=3 correctly fails
    assert wick_hi(3) > int(K4_threshold_for_Q10(3))


def test_Q_wick_formula_and_lt_10():
    th = theorem_Q_wick_lt_10()
    assert th["proved"] is True
    for p in (5, 7, 11, 13):
        n = n_of(p)
        qw = Q_from_wick_hi(p)
        assert qw == Fraction(6 * n + 52, n - 1)
        assert qw == 6 + Fraction(58, n - 1)
        assert qw < 10
        assert qw == Q_from_K4(p, wick_hi(p))
    # p=3: Q_Wick = 112/9 > 10
    assert Q_from_wick_hi(3) == Fraction(112, 9)
    assert Q_from_wick_hi(3) > 10


def test_partition_C_diag_and_eta_budget():
    th = theorem_K4_partition_identity()
    assert th["proved"] is True
    for p in (3, 5, 7, 11):
        n = n_of(p)
        cd = C_diag(p)
        # C_diag = 20n + 12n(2n-3)/(n-1)
        alt = Fraction(20 * n) + Fraction(12 * n * (2 * n - 3), n - 1)
        assert cd == alt
        bud = eta_sumsq_budget_for_wick_hi(p)
        assert bud == Fraction(n * (13 * n - 10), 6 * (n - 1))
        # Wick_hi = wick_lo + 96n; equality case uses coll + 24*bud = 96n
        assert wick_lo(p) + cd + 24 * bud == wick_hi(p)


def test_census_exact_K4_drive_shipped():
    """Drive certify_K4_exact (shipped fractions / p=7 census identity)."""
    c3 = certify_K4_exact(3)
    assert Fraction(c3["K4"]) == 1680
    assert c3["K4_le_Wick"] is True
    assert c3["K4_over_Wick"] == 1.0
    assert c3["K4_le_thr"] is False
    assert c3["eta_le_budget"] is True  # equality

    c5 = certify_K4_exact(5)
    assert Fraction(c5["K4"]) == Fraction(120400, 13)
    assert c5["K4_le_Wick"] is True
    assert c5["K4_le_thr"] is True
    assert c5["eta_le_budget"] is True
    assert Fraction(c5["Q"]) == Fraction(34512, 4225)

    c7 = certify_K4_exact(7)
    assert Fraction(c7["K4"]) == Fraction(5218435600, 167281)
    assert c7["K4_le_Wick"] is True
    assert c7["K4_le_thr"] is True
    assert c7["eta_le_budget"] is True
    assert float(c7["K4_over_Wick"]) < 1.0


def test_predicates_still_open():
    assert K4_bound_proved_general() is False
    assert residual_i_closed_via_K4_wick() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_198()
    assert h["wick_hi_le_thr_p_ge_5"] is True
    assert h["Q_wick_lt_10_p_ge_5"] is True
    assert h["K4_le_Wick_hi_general"] is False
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["wick_hi_le_thr_p_ge_5"] is True
    assert out["proved"]["K4_bound_general"] is False
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
