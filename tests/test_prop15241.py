"""Prop 15.241: hull criterion; halfspace insufficient; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15186 import mu_part_formula
from e1_gmin_m4_prop15188 import two_over_n
from e1_gmin_m4_prop15241 import (
    f4_formula,
    hinge_status_241,
    residual_i_closed_via_241,
    theorem_halfspace_orbit_insufficient,
    theorem_hull_criterion,
    theorem_segment_convexity,
)


def test_hull_criterion_proved_hypothesis_open():
    t = theorem_hull_criterion()
    assert t["proved_criterion"] is True
    assert t["hypothesis_proved_general"] is False
    assert t["sufficient_for_residual_i"] is False
    for p, row in t["by_p_sample"].items():
        assert row["le_two_over_n"] is True


def test_segment_convexity_max_plus_free():
    t = theorem_segment_convexity()
    assert t["proved"] is True
    # direct spot check p=5,7
    for p in (5, 7, 11):
        M = two_over_n(p)
        for k in (-1, 1):
            for ph in (-2 * (p - 2), -2, 2, 2 * (p - 2)):
                mp = mu_part_formula(p, k, ph)
                ff = f4_formula(p, k, ph)
                for tnum in range(0, 5):
                    tval = Fraction(tnum, 4)
                    x = (1 - tval) * mp + tval * ff
                    assert abs(x) <= M


def test_hull_endpoints_in_two_over_n_ball():
    """Any hull endpoint among {0,μ_part,f4} has abs ≤ 2/n."""
    for p in (5, 7, 11, 13, 17):
        M = two_over_n(p)
        bound = 2 * (p - 2)
        for k in (-1, 1):
            for ph in range(-bound, bound + 1, 4):
                mp = mu_part_formula(p, k, ph)
                ff = f4_formula(p, k, ph)
                lo = min(Fraction(0), mp, ff)
                hi = max(Fraction(0), mp, ff)
                assert abs(lo) <= M
                assert abs(hi) <= M
                assert abs(mp) <= M  # maj ≤ 2/n on |κ|=1,|φ|≤2(p−2)
                assert abs(ff) <= M


def test_halfspace_orbit_marked_dead():
    t = theorem_halfspace_orbit_insufficient()
    assert t["proved_structure"] is True
    assert "halfspace" in t["dead_path"]
    c5 = t["certified"]["5"]
    assert c5["exceeds_thr"] is True
    assert Fraction(c5["m4_H_max_abs_lower_bound"]) > Fraction(c5["thr"])
    c7 = t["certified"]["7"]
    assert Fraction(c7["m4_H_max_abs_lower_bound"]) > Fraction(c7["thr"])


def test_predicates_stay_false():
    assert residual_i_closed_via_241() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
    h = hinge_status_241()
    assert h["hull_hypothesis"] is False
    assert h["residual_i_closed_via_241"] is False
    assert h["segment_convexity"] is True
    assert h["hull_criterion"] is True


def test_main_writes_evidence(tmp_path=None):
    from e1_gmin_m4_prop15241 import main

    out = main()
    assert out["prop"] == "15.241"
    assert out["residual_i_closed"] is False
    assert out["proved"]["hull_criterion_chain"] is True
    assert out["proved"]["hull_hypothesis"] is False
