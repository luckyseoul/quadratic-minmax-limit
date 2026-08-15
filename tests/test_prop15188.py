"""Tests for Prop 15.188 — residual-(i) hinge honesty; 2/n target; no soft-close."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part, mu_part_formula  # noqa: E402
from e1_gmin_m4_prop15187 import theorem_global_T2_kappa  # noqa: E402
from e1_gmin_m4_prop15188 import (  # noqa: E402
    certified_p5_mu_on_kappa1,
    certified_p7_mu_census_constants,
    theorem_two_over_n_le_mu4_thr,
    two_over_n,
)


def test_two_over_n_beats_thr_all_sample_primes():
    r = theorem_two_over_n_le_mu4_thr([5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    for p in (5, 7, 11):
        assert two_over_n(p) <= mu4_threshold(p)


def test_two_over_n_fails_at_p3_as_expected():
    # 2/10 = 1/5 < 1/3 = thr? Wait thr at p=3 is 1/6; 2/10=1/5 > 1/6
    # so 2/n does NOT beat thr at p=3 — correct for residual scope p≥5 only.
    assert two_over_n(3) > mu4_threshold(3)


def test_p5_certified_beats_thr_and_two_over_n():
    r = certified_p5_mu_on_kappa1()
    assert r["beats_mu4_thr"] is True
    assert r["beats_two_over_n"] is True
    assert Fraction(r["max_abs_mu"]) == Fraction(3, 65)


def test_p7_certified_beats_thr_and_two_over_n():
    r = certified_p7_mu_census_constants()
    assert r["beats_mu4_thr"] is True
    assert r["beats_two_over_n"] is True
    assert Fraction(r["max_abs_mu_kappa1"]) == Fraction(109, 2863)
    assert r["not_bounded_by_mu_part_pointwise"] is True
    # 109/2863 > majorant? no — max act > max part
    assert Fraction(109, 2863) > majorant_mu_part(7)


def test_mu_part_formula_majorant_still_ok():
    # majorant algebra unchanged
    assert majorant_mu_part(5) <= mu4_threshold(5)
    assert majorant_mu_part(7) <= mu4_threshold(7)
    assert abs(mu_part_formula(7, 1, -10)) == majorant_mu_part(7)


def test_global_T2_kappa_still_proved():
    r = theorem_global_T2_kappa()
    assert r["proved"] is True


def test_predicates_remain_open_no_soft_close():
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is True
