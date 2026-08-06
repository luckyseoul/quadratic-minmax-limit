"""Prop 15.175 — Gμ reformulation; bound still open (no soft-close)."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15175 import (
    certify_Gmu_small_p,
    gmu_bound_implies_H,
    hinge_status_175,
    main,
)


def test_implication_fraction():
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        r = gmu_bound_implies_H(p)
        assert r["proved_implication"] is True
        assert r["bound_proved_general"] is False
        assert r["gmu_abs_thr"] == str(Fraction(1, p))


def test_census_gmu_small_p():
    c3 = certify_Gmu_small_p(3)
    c5 = certify_Gmu_small_p(5)
    assert c3["E2_is_I"] is True
    assert c5["E2_is_I"] is True
    assert c3["Gmu_wedge_zero"] is True
    assert c5["Gmu_wedge_zero"] is True
    assert c3["gmu_le_1_over_p"] is True
    assert c5["gmu_le_1_over_p"] is True
    assert c3["H_ge_minus_1_over_p"] is True
    assert c5["H_ge_minus_1_over_p"] is True
    assert c3["not_general_proof"] is True
    assert c5["not_general_proof"] is True


def test_predicates_still_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_175()
    assert h["bound_proved_general"] is False
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["bound_Gmu_le_1_over_p_general"] is False
    assert out["proved"]["gsum_disj_lb_proved_general"] is False
