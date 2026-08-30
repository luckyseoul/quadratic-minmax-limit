"""Tests for Prop 15.338 — named k-vector; interpolant p=7-only."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15338 import (
    alpha_interpolant,
    k_high,
    k_of_alpha,
    main,
    n_high,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    size_from_k,
)


def test_kvec_specialises_and_sizes():
    A = prove_A()
    assert A["proved"] is True
    assert A["matches_15337"] is True
    assert n_high(7) == 1
    assert n_high(11) == 2
    assert n_high(11) != 1
    assert size_from_k(7) == 21
    assert size_from_k(11) == 55
    assert size_from_k(19) == 171


def test_t2_alpha_levels():
    B = prove_B()
    assert B["proved"] is True
    assert k_of_alpha(7, 0) == 3
    assert k_of_alpha(7, 1) == 2
    assert k_of_alpha(7, 1) != 1
    assert k_of_alpha(7, 3) == 1
    assert k_high(11) == 3
    assert k_of_alpha(11, 1) == 3


def test_interpolant_p7_only():
    C = prove_C()
    assert C["proved"] is True
    a7 = {int(k): int(v) for k, v in C["alpha7"].items()}
    assert a7 == {1: 5, 2: 6, 3: 4, 4: 6, 5: 3, 6: 6}
    assert C["by_p"]["7"]["maxplus"] is True
    assert C["by_p"]["11"]["maxplus"] is False
    assert alpha_interpolant(7, 3) == 4


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["pair_rule_Maxplus_general"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["kvec_specialises_15337"] is True
    assert out["proved"]["t2_alpha_three_levels"] is True
    assert out["proved"]["interpolant_p7_only"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
