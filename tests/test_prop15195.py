"""Prop 15.195 — mass-corrected dual-eq obstruction; residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15190 import alpha_particular
from e1_gmin_m4_prop15194 import dual_eq_need_kappa_e
from e1_gmin_m4_prop15195 import (
    census_mass_min,
    hinge_status_195,
    main,
    mass_constrained_min_sum,
    mass_corrected_bound_proved_general,
    residual_i_closed_via_mass_corrected,
    theorem_mass_corrected_criterion,
    theorem_weak_floors_do_not_suffice,
)


def test_criterion_proved():
    assert theorem_mass_corrected_criterion()["proved"] is True
    assert theorem_weak_floors_do_not_suffice()["proved"] is True


def test_mass_min_algorithm_toy():
    """Drive real mass_constrained_min_sum on a tiny feasible instance."""
    alpha = Fraction(1, 2)
    M = Fraction(1)
    # lo=-1/2, n=4, all-lo sum=-2, B=-1, rem=1
    m = mass_constrained_min_sum(
        [Fraction(-2), Fraction(3)], n_wedge=2, M=M, alpha=alpha
    )
    # raise lowest (-2) by rem=1 → that κ = -1/2+1 = 1/2; contrib (-2)*(1/2)=-1
    # others at -1/2: costs 3,0,0 → contrib 3*(-1/2)+0+0 = -3/2
    # total = -1 + -3/2 = -5/2
    assert m == Fraction(-5, 2)


def test_census_p5_exact_blocks():
    c = census_mass_min(5)
    assert c.get("skipped") is not True
    assert c["exact"] is True
    msum = Fraction(c["mass_min_sum_a_kappa"])
    target = Fraction(c["minus_2M"])
    assert msum == Fraction(-30, 13)
    assert msum > target
    assert c["blocks_dual_eq"] is True
    # target matches −2(2−α)
    M = dual_eq_need_kappa_e(5)
    assert target == -2 * M


def test_census_p3_does_not_block():
    c = census_mass_min(3)
    assert c["blocks_dual_eq"] is False


def test_predicates_open():
    assert mass_corrected_bound_proved_general() is False
    assert residual_i_closed_via_mass_corrected() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_195()
    assert h["mass_corrected_criterion_proved"] is True
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["L_closed"] is False
    assert out["certified"]["blocks_p5_p7"] is True
