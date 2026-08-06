"""Tests for Prop 15.182 — dual-equality normal form; residual i still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)
from e1_gmin_m4_prop15182 import (  # noqa: E402
    alpha_particular,
    certify_ker_box_infeasible_p5,
    theorem_general_solution_form,
    theorem_need_matches_half,
    theorem_particular_solution,
)


def test_alpha_and_mass_identity():
    for p in [5, 7, 11, 13, 17, 19, 23]:
        if not is_prime(p):
            continue
        n = n_of(p)
        alpha = alpha_particular(p)
        assert alpha == Fraction(6, p * n)
        E = n * (n - 1) // 2
        mass = alpha * E - 2
        assert mass == k_type_I_3p_minus_2(p)
        assert alpha - 2 < 0


def test_particular_solution_theorem():
    r = theorem_particular_solution([p for p in range(5, 60) if is_prime(p)])
    assert r["proved"] is True


def test_general_form_proved():
    assert theorem_general_solution_form()["proved"] is True


def test_need_matches():
    assert theorem_need_matches_half()["proved"] is True


def test_p5_ker_box_infeasible_census():
    r = certify_ker_box_infeasible_p5()
    assert bool(r["ker_box_infeasible"]) is True
    assert r["not_general_proof"] is True


def test_hinge_still_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
