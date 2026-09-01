"""Tests for Prop 15.577 — Type+ 1D Johnson kernel; DEAD as Type I kill."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15577 import (
    S_particular,
    S_particular_drop2,
    deg2_ambient_dim,
    johnson_ker_dim,
    johnson_rank_named,
    main,
    p5_size_shadow,
    particular_k,
    prove_A,
    prove_B,
    prove_open,
    type_I_k,
)


def test_johnson_kernel_rank():
    A = prove_A()
    assert A["proved"] is True
    assert deg2_ambient_dim(5) == 16
    assert johnson_rank_named(5) == 10
    assert johnson_ker_dim(5) == 6
    assert johnson_rank_named(7) == 21
    assert johnson_ker_dim(7) == 8
    assert A["live"]["5"]["rank"] == 10
    assert A["live"]["7"]["rank"] == 21
    # fail-when-wrong: forget kernel
    assert johnson_rank_named(5) != deg2_ambient_dim(5)
    assert johnson_rank_named(5) != 9


def test_zero_star_particular():
    B = prove_B()
    assert B["proved"] is True
    assert particular_k(5) == 17
    assert type_I_k(5) == 13
    assert particular_k(5) != type_I_k(5)
    assert particular_k(7) == 37 != type_I_k(7) == 19
    eps = [1, 1, 1, -1, -1]  # sum=1
    assert S_particular(eps) == 3 - 2 * eps[0]
    assert S_particular(eps) != 0
    assert S_particular_drop2(eps) != 3 - 2 * eps[0]


def test_p5_size13_shadow_and_flags_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["path_dead"] is True
    sh = p5_size_shadow()
    assert sh["ok"] is True
    assert sh["k"] == 13
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["johnson_ker_p_plus_1"] is True
    assert out["proved"]["zero_star_particular"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
