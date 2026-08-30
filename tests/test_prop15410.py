"""Tests for Prop 15.410 — zero-far gap-2; p=7 bad-case box empty."""
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
from e1_gmin_m4_prop15410 import (
    S_plus_3fe_on_fplus,
    Y_plus_zerofar,
    main,
    mu_box,
    mu_boxt,
    prove_A,
    prove_B,
    prove_open,
)
from e1_gmin_m4_prop15408 import badcase_ub


def test_zerofar_named_not_badcase():
    A = prove_A()
    assert A["proved"] is True
    assert Y_plus_zerofar(5) == Fraction(-3, 2)
    assert Y_plus_zerofar(7) == Fraction(-4, 3)
    assert S_plus_3fe_on_fplus(5) == Fraction(3, 2)
    assert S_plus_3fe_on_fplus(7) == Fraction(5, 3)
    assert S_plus_3fe_on_fplus(7) > 0
    assert Y_plus_zerofar(7) < -1
    assert mu_box(7) == Fraction(1, 12)
    assert mu_boxt(7) == Fraction(5, 16)


def test_p7_uses_shared_correct_badcase_encoding():
    Fm = np.array([[1.0, -1.0], [-1.0, -1.0]])
    fe = np.array([1.0, -1.0])
    x = np.array([1.0, 1.0])
    A, b = badcase_ub(Fm, fe)
    np.testing.assert_allclose(A[2:] @ x - b[2:], Fm @ x + 3.0 * fe)


def test_zerofar_cache_typeI():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["plus"] == [1, 5]
    assert B["by_p"]["7"]["plus"] == [1, 5]
    assert B["by_p"]["7"]["cut"] > 0


def test_type_I_flags_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["zerofar_gap2_not_badcase"] is True
    assert out["proved"]["p5_p7_zerofar_cache"] is True
    assert out["proved"]["p7_badcase_box_empty"] is True
    assert out["algebra"]["C"]["badcase_lp_feasible"] is False
    assert out["algebra"]["C"]["badcase_status"] == 2
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["prop"] == "15.410"
