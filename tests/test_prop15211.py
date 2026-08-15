"""Tests for Prop 15.211 — G+ PSD / lambda_star; hinge open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
)
from e1_gmin_m4_prop15211 import (  # noqa: E402
    Gplus_pd_proved_general,
    certify_angle_majorant_p5_only,
    certify_comm_dual_structure,
    certify_rayleigh_vs_lambda_star,
    free_e_bound_proved_general,
    hinge_status_211,
    lambda_star,
    residual_i_closed_via_211,
    theorem_Gplus_psd_on_subspace,
    theorem_alpha_n_minus_1_dual_target,
    theorem_lambda_star_positive,
    theorem_pd_criterion,
)


def test_Gplus_psd_proved():
    assert theorem_Gplus_psd_on_subspace()["proved"]


def test_pd_criterion_proved():
    assert theorem_pd_criterion()["proved"]


def test_lambda_star_positive():
    t = theorem_lambda_star_positive()
    assert t["proved"]
    assert lambda_star(26) == Fraction(80, 13)
    assert lambda_star(50) == Fraction(176, 25)
    assert lambda_star(n_of(5)) > 0


def test_alpha_n_minus_1_target():
    assert theorem_alpha_n_minus_1_dual_target()["proved"]


def test_rayleigh_cert_sharp_p5():
    c = certify_rayleigh_vs_lambda_star()
    assert c["certified"]
    assert c["proved_general"] is False
    assert c["by_p"]["5"]["sharp"] is True
    assert c["by_p"]["5"]["min_rayleigh"] == str(Fraction(80, 13))
    assert c["by_p"]["7"]["min_ge_lambda_star"] is True


def test_comm_dual_cert():
    c = certify_comm_dual_structure()
    assert c["certified"]
    assert c["by_p"]["5"]["n_minus_1"] == 25
    assert abs(c["by_p"]["5"]["sum_ne_after_comm_scale"] - 25) < 0.1


def test_angle_majorant_p5_only():
    c = certify_angle_majorant_p5_only()
    assert c["by_p"]["5"]["ub_le_thr"] is True
    assert c["by_p"]["7"]["ub_le_thr"] is False


def test_predicates_honest_false():
    assert free_e_bound_proved_general() is False
    assert Gplus_pd_proved_general() is False
    assert residual_i_closed_via_211() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_211()
    assert h["residual_i_closed"] is False
    assert h["Gplus_psd_on_Wpp0_proved"] is True
    assert h["Gplus_pd_general"] is False
