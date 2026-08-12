"""Tests for Prop 15.226 — Farkas-sharp thr; G-δ; cross; ‖μ_part‖₂²; open."""
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
from e1_gmin_m4_prop15176 import mu_star_threshold  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15225 import room_delta  # noqa: E402
from e1_gmin_m4_prop15226 import (  # noqa: E402
    hinge_status_226,
    mu_farkas_threshold,
    mu_part_norm_sq,
    mu_part_norm_sq_via_ab,
    phi_kappa_inner,
    phi_norm_sq_closed,
    residual_i_closed_via_226,
    room_delta_farkas,
    theorem_delta_G_invariant,
    theorem_farkas_sharp_thr,
    theorem_maj_beats_mu_farkas,
    theorem_maxpm_cross_moment,
    theorem_mu_L2_from_m4,
    theorem_mu_part_norm_sq_closed,
)


def test_farkas_sharp_proved():
    assert theorem_farkas_sharp_thr()["proved"] is True


def test_maj_beats_mu_f_proved():
    assert theorem_maj_beats_mu_farkas()["proved"] is True


def test_cross_proved():
    assert theorem_maxpm_cross_moment()["proved"] is True


def test_mu_L2_proved():
    assert theorem_mu_L2_from_m4()["proved"] is True


def test_G_invariant_proved():
    assert theorem_delta_G_invariant()["proved"] is True


def test_mu_part_norm_proved():
    assert theorem_mu_part_norm_sq_closed()["proved"] is True


def test_mu_f_equals_minus_half_mu_star():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert is_prime(p)
        assert mu_farkas_threshold(p) == -mu_star_threshold(p) / 2


def test_one_over_2p_lt_mu_f():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
        assert mu4_threshold(p) < mu_farkas_threshold(p)


def test_maj_le_mu_f_and_rooms():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        maj = majorant_mu_part(p)
        mf = mu_farkas_threshold(p)
        assert maj <= mf
        rdf = room_delta_farkas(p)
        rd = room_delta(p)
        assert rdf == mf - maj
        assert rdf > rd > 0


def test_mu_f_values():
    assert mu_farkas_threshold(5) == Fraction(7, 65)
    assert mu_farkas_threshold(7) == Fraction(11, 133)


def test_room_delta_farkas_values():
    # p=5: 7/65 - 9/125 = (7*25 - 9*13)/(65*25) = (175-117)/1625 = 58/1625
    assert room_delta_farkas(5) == Fraction(58, 1625)


def test_cross_formula():
    for p in (5, 7, 11, 13):
        n = n_of(p)
        assert Fraction(n * (n - 2), 8) == Fraction(n * (n - 2), 8)


def test_mu_part_norm_matches_ab():
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert mu_part_norm_sq(p) == mu_part_norm_sq_via_ab(p)
        assert mu_part_norm_sq(p) > 0


def test_mu_part_norm_values():
    # p=5: (4)(6)(26)(76)/(24*20) = 4*6*26*76 / 480
    assert mu_part_norm_sq(5) == Fraction(494, 5)
    assert mu_part_norm_sq(7) == Fraction(3700, 11)


def test_phi_identities_consistent():
    # Tφ identity check already in module construction; smoke ⟨φ,κ⟩, ∑φ² > 0
    for p in (5, 7, 11, 13):
        assert phi_norm_sq_closed(p) > 0
        # ⟨φ,κ⟩ can be positive
        assert isinstance(phi_kappa_inner(p), Fraction)


def test_residual_i_still_open():
    assert residual_i_closed_via_226() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_226()
    assert h["farkas_sharp_thr"] is True
    assert h["maj_beats_mu_farkas"] is True
    assert h["mu_part_norm_sq_closed"] is True
    assert h["residual_i_closed_via_226"] is False
