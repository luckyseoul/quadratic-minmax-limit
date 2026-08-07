"""Tests for Prop 15.209 — G+ spectrum / star saturation / dim; hinge open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15190 import alpha_particular, scheme_ker_max_kappa_e  # noqa: E402
from e1_gmin_m4_prop15205 import free_e_scheme_ratio_candidate  # noqa: E402
from e1_gmin_m4_prop15209 import (  # noqa: E402
    Gplus_pd_proved_general,
    certify_Gplus_spectrum,
    certify_free_e_star_saturation,
    dim_zero_diag_plusplus,
    free_e_bound_proved_general,
    hinge_status_209,
    residual_i_closed_via_209,
    theorem_dim_zero_diag_plusplus,
    theorem_implication_Gplus_to_residual_i,
    theorem_star_count,
)


def test_dim_formula_proved():
    t = theorem_dim_zero_diag_plusplus()
    assert t["proved"]
    assert dim_zero_diag_plusplus(5) == 65
    assert dim_zero_diag_plusplus(7) == 275
    assert dim_zero_diag_plusplus(3) == 5


def test_star_count_proved():
    t = theorem_star_count()
    assert t["proved"]
    assert t["by_p_sample"]["5"]["n_stars"] == 48
    assert t["by_p_sample"]["7"]["n_stars"] == 96


def test_implication_recorded():
    t = theorem_implication_Gplus_to_residual_i()
    assert t["proved"]


def test_Gplus_spectrum_cert_positive():
    c = certify_Gplus_spectrum()
    assert c["certified"]
    assert c["proved_general"] is False
    for p, row in c["by_p"].items():
        assert row["positive_definite"]
        assert row["dim_Wpp0"] == dim_zero_diag_plusplus(int(p))
    # exact fractions
    assert c["by_p"]["5"]["Gplus_min_eig"] == str(Fraction(40, 13))
    assert c["by_p"]["7"]["Gplus_min_eig"] == str(Fraction(1536, 409))
    # multiplicities sum to dim
    for p, row in c["by_p"].items():
        mult_sum = sum(row["moment_eigs"].values())
        assert mult_sum == row["dim_Wpp0"]


def test_free_e_star_saturation_cert():
    c = certify_free_e_star_saturation()
    assert c["certified"]
    assert c["proved_general"] is False
    for p in (5, 7):
        row = c["by_p"][str(p)]
        assert row["match_r_scheme"]
        assert row["stars_fixed_equals_unrestricted"]
        fe = Fraction(row["free_e_max"])
        r = free_e_scheme_ratio_candidate(p)
        sch = scheme_ker_max_kappa_e(p)
        assert fe == r * sch
        assert fe < 2 - alpha_particular(p)


def test_predicates_honest_false():
    assert free_e_bound_proved_general() is False
    assert Gplus_pd_proved_general() is False
    assert residual_i_closed_via_209() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_209()
    assert h["residual_i_closed"] is False
    assert h["Gplus_pd_general"] is False
