"""Prop 15.197 — min-distance + Q(K4); residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15197 import (
    K4_bound_proved_general,
    K4_threshold_for_Q10,
    Q_from_K4,
    certify_K4_and_distance,
    hinge_status_197,
    main,
    max_dot_nonopposite,
    residual_i_closed_via_K4,
    theorem_Q_from_K4_identity,
    theorem_min_distance_p_plus_1,
)


def test_min_distance_theorem():
    th = theorem_min_distance_p_plus_1()
    assert th["proved"] is True
    for p in (3, 5, 7, 11, 13, 17, 19):
        assert max_dot_nonopposite(p) == (p - 1) ** 2 - 2
        assert max_dot_nonopposite(p) == n_of(p) - 2 * (p + 1)


def test_Q_K4_identity_threshold():
    th = theorem_Q_from_K4_identity()
    assert th["proved"] is True
    for p in (5, 7, 11, 13):
        thr = K4_threshold_for_Q10(p)
        assert thr == 16 * n_of(p) ** 2 - 14 * n_of(p)
        assert Q_from_K4(p, thr) == 10


def test_census_distance_and_K4():
    """Drive real certify_K4_and_distance on shipped path."""
    c5 = certify_K4_and_distance(5)
    assert c5.get("skipped") is not True
    assert Fraction(c5["K4"]) == Fraction(120400, 13)
    assert c5["K4_le_thr"] is True
    assert c5["max_dot_off"] == c5["max_dot_theory"] == 14
    assert Fraction(c5["Q_from_K4"]) == Fraction(34512, 4225)
    assert abs(c5["K2"] - c5["K2_theory_2n"]) < 1e-6

    c3 = certify_K4_and_distance(3)
    assert Fraction(c3["K4"]) == 1680
    assert c3["K4_le_thr"] is False  # p=3 correctly above thr
    assert c3["max_dot_off"] <= c3["max_dot_theory"]

    c7 = certify_K4_and_distance(7)
    if not c7.get("skipped"):
        assert c7["K4_le_thr"] is True
        assert c7["max_dot_off"] == 34
        assert c7["distance_bound_sharp_or_ok"] is True


def test_predicates_open():
    assert K4_bound_proved_general() is False
    assert residual_i_closed_via_K4() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_197()
    assert h["min_distance_p_plus_1"] is True
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["min_distance_p_plus_1"] is True
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
