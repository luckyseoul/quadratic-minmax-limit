"""Tests for Prop 15.224 — T²φ closed; μ_part solves master; residual-i open."""
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
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15224 import (  # noqa: E402
    T2_phi_formula_n,
    T2_phi_formula_p,
    hinge_status_224,
    residual_i_closed_via_224,
    theorem_T2_phi_closed,
    theorem_mu_part_solves_master_general,
    theorem_particular_majorant_recall,
)


def test_T2_phi_proved():
    assert theorem_T2_phi_closed()["proved"] is True


def test_mu_part_master_proved():
    assert theorem_mu_part_solves_master_general()["proved"] is True


def test_majorant_recall():
    assert theorem_particular_majorant_recall()["proved"] is True


def test_T2_phi_n_equals_p_form():
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        assert is_prime(p)
        n = n_of(p)
        for kappa, phi in [(1, 0), (1, 2), (-1, -4), (3, 4), (1, -6)]:
            assert T2_phi_formula_n(n, phi, kappa) == T2_phi_formula_p(
                p, phi, kappa
            )


def test_T2_phi_values():
    # n=26, p=5: 4(28)(φ-2κ)=112(φ-2κ)
    assert T2_phi_formula_n(26, 0, 1) == 112 * (0 - 2)
    assert T2_phi_formula_p(5, 0, 1) == 4 * (25 + 3) * (0 - 2)


def test_majorant_le_half_over_p():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        assert is_prime(p)
        assert majorant_mu_part(p) <= Fraction(1, 2 * p)


def test_residual_i_still_open():
    assert residual_i_closed_via_224() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_224()
    assert h["T2_phi_closed"] is True
    assert h["mu_part_solves_master"] is True
    assert h["residual_i_closed_via_224"] is False
