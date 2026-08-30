"""Prop 15.200 — C−2/n in ker; free-e dual-eq box; residual (i) open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15190 import alpha_particular, scheme_ker_max_kappa_e
from e1_gmin_m4_prop15200 import (
    certify_C_minus_2n_ker,
    certify_free_e_kerbox,
    hinge_status_200,
    kerbox_empty_proved_general,
    main,
    residual_i_closed_via_kerbox,
    theorem_C_minus_2n_in_ker,
    theorem_dual_eq_free_e_box,
    theorem_three_halves_scheme_blocks_if_max_le,
)


def test_theorems_proved_structure():
    assert theorem_C_minus_2n_in_ker()["proved"] is True
    assert theorem_dual_eq_free_e_box()["proved"] is True
    th = theorem_three_halves_scheme_blocks_if_max_le()
    assert th["proved_implication"] is True
    assert th["max_bound_proved"] is False
    for p in (5, 7, 11, 13, 17):
        sch = scheme_ker_max_kappa_e(p)
        assert Fraction(3, 2) * sch < 2 - alpha_particular(p)


def test_C_minus_2n_ker_drive_shipped():
    """Drive certify_C_minus_2n_ker on real Max± path."""
    for p in (3, 5):
        c = certify_C_minus_2n_ker(p)
        assert c.get("skipped") is not True
        assert c["in_ker"] is True
        assert c["||Gsum kappa||"] < 1e-8
        assert c["max_|f^T kappa|_Maxpm"] < 1e-8


def test_free_e_kerbox_census_drive_shipped():
    """Drive free-e ker-box LP: p=3 feasible scale, p=5 blocks."""
    c3 = certify_free_e_kerbox(3)
    assert c3["success"] is True
    assert c3["blocks_dual_eq"] is False  # max_ke > need at p=3
    assert c3["max_kappa_e"] > c3["need_2_minus_alpha"]

    c5 = certify_free_e_kerbox(5)
    assert c5["success"] is True
    assert c5["blocks_dual_eq"] is True
    assert c5["max_kappa_e"] < c5["need_2_minus_alpha"]
    # matches prior census scale ~0.81
    assert 0.7 < c5["max_kappa_e"] < 0.9


def test_predicates_open():
    assert kerbox_empty_proved_general() is False
    assert residual_i_closed_via_kerbox() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_200()
    assert h["C_minus_2n_in_ker"] is True
    assert h["kerbox_empty_general"] is False
    assert h["residual_i_closed"] is False
    out = main()
    assert out["proved"]["C_minus_2n_in_ker"] is True
    assert out["proved"]["kerbox_empty_general"] is False
    assert out["proved"]["L_closed"] is False
