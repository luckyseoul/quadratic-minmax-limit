"""Tests for Prop 15.455 — J2/R_aff alignments; residual is aff+J2."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15454 import is_Jr
from e1_gmin_m4_prop15455 import (
    D_from_p7_example,
    count_J2_alignments,
    is_axis_halfnet,
    main,
    n_CRRR_named,
    n_J2_pool_size,
    n_affJ2,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_J2_p5_is_nNL():
    A = prove_A()
    assert A["proved"] is True
    assert n_J2_pool_size(5) == 20
    assert count_J2_alignments(5, [2, 3, 4]) == 100
    assert A["n_all"] == 8000
    assert A["n_good"] == HPLUS[5] - n_1d(5) == 100
    assert A["n_all"] != A["n_good"]


def test_CRRR_p7():
    B = prove_B()
    assert B["proved"] is True
    assert n_CRRR_named(7) == 4 * 7 * 7 * 6 == 1176
    assert 4 * n_1d(7) == 560
    assert n_CRRR_named(7) != 4 * n_1d(7)
    assert n_CRRR_named(5) == 300
    assert n_CRRR_named(5) != HPLUS[5] - n_1d(5)
    D = D_from_p7_example()
    assert len(D) == 21
    assert is_axis_halfnet(D, 7, [1, 2]) is True


def test_residual_affJ2():
    C = prove_C()
    assert C["proved"] is True
    assert tuple(sorted(n_affJ2(7, 1, 0, 6, 0))) == (0, 1, 1, 2, 5, 5, 7)
    assert is_Jr(7, (0, 1, 1, 2, 5, 5, 7), rmax=2) is False


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.455"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["J2_names_nNL_p5"] is True
    assert out["proved"]["CRRR_named_p7"] is True
    assert out["proved"]["residual_is_affJ2"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
