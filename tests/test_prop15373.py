"""Tests for Prop 15.373 — symbolic Q++_form interval for all odd p≥5."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1
from e1_gmin_m4_prop15367 import A_num
from e1_gmin_m4_prop15371 import D_form, D_form_drop_binom
from e1_gmin_m4_prop15372 import Qpp_form
from e1_gmin_m4_prop15373 import (
    N_lo_closed,
    N_lo_from_atoms,
    N_up_closed,
    N_up_from_atoms,
    PA_PD,
    dC_lo,
    dC_up,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    sibling_gap_poly,
)


def test_A_and_D_positive():
    A = prove_A()
    assert A["proved"] is True
    assert A_num(5) == 6
    assert D_form(5) == 13
    assert A["A_lb5"] == 6


def test_closed_forms_match_atoms():
    for p in (5, 7, 11, 13, 17, 19, 23, 29):
        assert N_lo_from_atoms(p) == N_lo_closed(p)
        assert N_up_from_atoms(p) == N_up_closed(p)
        assert dC_lo(p) > 0
        assert dC_up(p) > 0
        assert N_lo_closed(p) < 0
        assert N_up_closed(p) < 0


def test_cross_multiply_is_the_interval():
    B = prove_B()
    assert B["proved"] is True
    assert N_lo_closed(5) == -1764
    for p in (5, 7, 13, 17):
        PA, PD = PA_PD(p)
        assert PA == 2 * A_num(p)
        assert PD == 2 * D_form(p)
        assert PD * (2 * p + 1) * (p - 3) <= 4 * PA * (p - 1) ** 2
        assert 16 * PA * (p * p - 2 * p - 1) <= PD * (9 * p * p - 20 * p - 21)
        Q = Qpp_form(p)
        assert Q_lo_named(p) <= Q <= Qpp_floor_ub(p)
    Qbad = Fraction(8 * A_num(5), D_form_drop_binom(5))
    assert Qbad == 16
    assert not (Q_lo_named(5) <= Qbad <= Qpp_floor_ub(5))
    # p≡3 N_up sign is proved, not waved
    q7 = 2 * 7 ** 4 + 7 * 7 ** 3 - 85 * 7 ** 2 - 23 * 7 + 291
    assert q7 == 3168
    assert N_up_closed(7) == -3 * q7
    assert N_up_closed(7) < 0


def test_remaining_sibling_identity():
    C = prove_C()
    assert C["proved"] is True
    assert sibling_gap_poly(5) == 72
    for p in (5, 13, 17, 29):
        lhs = 3 * p ** 3 - 7 * p ** 2 + 9 * p + 19
        rhs = 3 * (p - 1) ** 3
        assert lhs - rhs == 2 * p * p + 22
        assert even_S_p1(p, Q_lo_named(p))["sm_j0"] == 6
        assert even_S_p1(p, Qpp_floor_ub(p))["sm_joth"] == 6
        assert min(even_S_p1(p, Q_lo_named(p)).values()) >= 6
        assert min(even_S_p1(p, Qpp_form(p)).values()) >= 6


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["A_D_positive"] is True
    assert out["proved"]["Qform_in_interval_all_odd_p"] is True
    assert out["proved"]["p1_even_S_ge_6_on_interval"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
