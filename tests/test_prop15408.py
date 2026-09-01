"""Tests for Prop 15.408 — Type I Max+ 0-1 at p=5; gap-2 bad box empty."""
from __future__ import annotations

from fractions import Fraction

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15408 import (
    BAD_WIT,
    S_H_mean_minus,
    WIT,
    badcase_ub,
    main,
    n_box_zerofar_maxplus,
    n_box_zerofar_pairmin,
    n_boxt_zerofar_maxplus,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_E,
    prove_open,
)


def test_H_lift_max_minus_2():
    A = prove_A()
    assert A["proved"] is True
    assert S_H_mean_minus(5) == Fraction(-14, 5)
    assert S_H_mean_minus(5) != -2
    assert S_H_mean_minus(5) != -3
    assert (-1 + 1) == 0
    assert (-1 + -1) == -2


def test_zerofar_counts_swapped():
    B = prove_B()
    assert B["proved"] is True
    assert n_box_zerofar_maxplus(5) == 3
    assert n_boxt_zerofar_maxplus(5) == 10
    assert n_box_zerofar_pairmin(5) == 10
    assert n_box_zerofar_maxplus(5) != n_box_zerofar_pairmin(5)
    assert n_box_zerofar_maxplus(5) - n_boxt_zerofar_maxplus(5) == 3 - 10


def test_p5_typeI_01_exists():
    C = prove_C()
    assert C["proved"] is True
    assert C["plus"] == [1, 5]
    assert C["n_S1"] == 156
    assert C["n_S5"] == 104
    assert C["n_m1"] == 36
    assert C["n_m1"] < 53
    assert C["bad"] is False
    assert C["minus_ext"] == [-11, 11]
    assert len(WIT) == 13
    assert (0, 1) not in WIT


def test_p5_unrestricted_fe_not_empty():
    D = prove_D()
    assert D["proved"] is True
    assert D["bad"] is True
    assert D["n_m1"] == 24
    assert D["n_m1"] < 53
    assert D["s_max"] == 5
    assert D["s_max"] > -1
    assert len(BAD_WIT) == 13
    assert BAD_WIT != WIT


def test_p5_gap2_badcase_box_empty():
    E = prove_E()
    assert E["proved"] is True
    assert E["gap2_lp_feasible"] is True
    assert E["badcase_lp_feasible"] is False
    assert E["gap2_Y_fplus"] > -3


def test_badcase_rows_encode_S_plus_3fe_not_S_plus_3kfe():
    Fm = np.array([[1.0, -1.0, 1.0], [-1.0, 1.0, 1.0]])
    fe = np.array([1.0, -1.0])
    x = np.array([1.0, 1.0, 0.0])
    A, b = badcase_ub(Fm, fe)
    S = Fm @ x
    np.testing.assert_allclose(A[:2] @ x - b[:2], S + 1.0)
    np.testing.assert_allclose(A[2:] @ x - b[2:], S + 3.0 * fe)
    # The old buggy rows would contain the extra factor sum(x)=2.
    old_residual = (Fm + 3.0 * fe[:, None]) @ x
    assert not np.allclose(old_residual, S + 3.0 * fe)


def test_type_I_flags_still_open():
    F = prove_open()
    assert F["proved"] is False
    assert F["type_I_multilevel_bad_case_ND_closed"] is True
    assert F["type_I_aut_e_3AB_positive_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_main():
    out = main()
    assert out["proved"]["H_lift_max_minus_2"] is True
    assert out["proved"]["zerofar_counts_swapped"] is True
    assert out["proved"]["p5_typeI_01_exists"] is True
    assert out["proved"]["p5_typeI_unrestricted_fe_not_empty"] is True
    assert out["proved"]["p5_gap2_badcase_box_empty"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
    assert out["prop"] == "15.408"
