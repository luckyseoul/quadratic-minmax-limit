"""Tests for Prop 15.186 — |φ| all p; μ_part majorant; residual i open."""
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
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import (  # noqa: E402
    certify_kappa_cr_dichotomy,
    certify_lsqr_equals_particular,
    majorant_mu_part,
    mu_part_formula,
    theorem_H_mod_4,
    theorem_hasse,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
    theorem_supersingular_double_residue,
)


def test_H_mod_4():
    assert theorem_H_mod_4()["proved"] is True


def test_hasse():
    assert theorem_hasse()["proved"] is True


def test_supersingular_double_residue_proved():
    assert theorem_supersingular_double_residue()["proved"] is True


def test_phi_bound_all_p_proved():
    assert theorem_phi_bound_general()["proved"] is True


def test_particular_majorant_beats_thr():
    r = theorem_particular_majorant_beats_thr(
        [p for p in range(5, 60) if is_prime(p)]
    )
    assert r["proved"] is True
    assert r["cubic_positive"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        maj = majorant_mu_part(p)
        thr = mu4_threshold(p)
        assert maj <= thr
        # explicit max over phi ladder
        bound_phi = 2 * (p - 2)
        phis = list(range(-bound_phi, bound_phi + 1, 4))
        mx = max(abs(mu_part_formula(p, k, ph)) for k in (-1, 1) for ph in phis)
        assert mx <= maj


def test_majorant_formula_p5():
    assert majorant_mu_part(5) == Fraction(36, 500)
    assert majorant_mu_part(5) <= Fraction(1, 10)


def test_lsqr_particular_p5():
    r = certify_lsqr_equals_particular(5)
    assert r["equals_particular"] is True
    assert r["lsqr_le_thr"] is True


def test_dichotomy_p5():
    r = certify_kappa_cr_dichotomy(5, 600)
    assert r["dichotomy_holds"] is True


def test_hinge_still_open():
    """Actual Max± |μ₄| bound not proved ⇒ residual i / E1 stay open."""
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
