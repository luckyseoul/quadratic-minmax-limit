"""Tests for Prop 15.225 — constrained P±; φ⊥E; μ_part=μ_mn; residual-i open."""
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
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15225 import (  # noqa: E402
    hinge_status_225,
    mu_avg_part_coeff,
    projector_minus_on_eigen,
    projector_plus_on_eigen,
    residual_i_closed_via_225,
    room_delta,
    room_delta_via_thr_maj,
    theorem_L2_orthogonality,
    theorem_avg_part_beats_thr,
    theorem_constrained_projectors,
    theorem_delta_split,
    theorem_mu_part_equals_mu_mn,
    theorem_phi_perp_Epm4p,
    theorem_room_delta_closed,
)


def test_projectors_proved():
    assert theorem_constrained_projectors()["proved"] is True


def test_delta_split_proved():
    assert theorem_delta_split()["proved"] is True


def test_phi_perp_proved():
    assert theorem_phi_perp_Epm4p()["proved"] is True


def test_mu_part_equals_mu_mn_proved():
    assert theorem_mu_part_equals_mu_mn()["proved"] is True


def test_avg_part_proved():
    assert theorem_avg_part_beats_thr()["proved"] is True


def test_room_delta_proved():
    assert theorem_room_delta_closed()["proved"] is True


def test_L2_orth_proved():
    assert theorem_L2_orthogonality()["proved"] is True


def test_projector_eigenvalues():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        assert is_prime(p)
        assert projector_plus_on_eigen(4 * p, p) == 1
        assert projector_minus_on_eigen(4 * p, p) == 0
        assert projector_plus_on_eigen(-4 * p, p) == 0
        assert projector_minus_on_eigen(-4 * p, p) == 1
        # Sum to identity on both eigenspaces
        assert (
            projector_plus_on_eigen(4 * p, p)
            + projector_minus_on_eigen(4 * p, p)
            == 1
        )
        assert (
            projector_plus_on_eigen(-4 * p, p)
            + projector_minus_on_eigen(-4 * p, p)
            == 1
        )


def test_room_delta_matches_thr_minus_maj():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
        assert is_prime(p)
        r = room_delta(p)
        assert r == room_delta_via_thr_maj(p)
        assert r == mu4_threshold(p) - majorant_mu_part(p)
        assert r > 0


def test_room_delta_closed_form_values():
    # p=5: thr=1/10, maj=9/125, room=7/250
    assert room_delta(5) == Fraction(7, 250)
    # p=7: thr=1/14, maj=17/539, room=43/1078
    assert room_delta(7) == Fraction(43, 1078)


def test_avg_coeff_beats_thr():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        assert mu_avg_part_coeff(p) <= mu4_threshold(p)
        assert mu_avg_part_coeff(p) < majorant_mu_part(p)


def test_avg_coeff_values():
    assert mu_avg_part_coeff(5) == Fraction(6, 125)
    assert mu_avg_part_coeff(7) == Fraction(12, 539)


def test_residual_i_still_open():
    assert residual_i_closed_via_225() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_225()
    assert h["constrained_projectors"] is True
    assert h["phi_perp_Epm4p"] is True
    assert h["mu_part_equals_mu_mn"] is True
    assert h["room_delta_closed"] is True
    assert h["residual_i_closed_via_225"] is False
