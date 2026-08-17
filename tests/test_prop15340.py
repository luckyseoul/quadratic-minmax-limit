"""Tests for Prop 15.340 — p=11 15.338 k-vector Max+ exists."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15340 import (
    PAIR,
    build_y,
    energy,
    kvec_of,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    roles_of,
)


def test_rebuild_is_maxplus():
    y = build_y()
    assert energy(y) < 1e-10
    assert int((y[1:] < 0).sum()) == 55
    assert y[0] > 0 and y[1] < 0


def test_kvector_and_roles():
    A = prove_A()
    assert A["proved"] is True
    assert A["k_hist"] == {"0": 1, "2": 8, "3": 2, "5": 1}
    assert A["roles"]["full"] == ["s3"]
    assert A["roles"]["empty"] == ["s9"]
    assert set(A["roles"]["high"]) == {"s0", "s8"}
    # fail-when-wrong: opposite-high {0,∞}
    assert set(A["roles"]["full"] + A["roles"]["empty"]) != {"s0", "inf"}


def test_not_interpolant_has_nc():
    B = prove_B()
    assert B["proved"] is True
    assert B["eq_interpolant"] is False
    assert B["interpolant_size"] == 59
    assert B["halfspace"] is False
    assert B["has_nc"] is True


def test_no_deg2_square_predicate():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_deg2"] == 0


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["pair_rule_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_pair_table_complete():
    assert set(PAIR) == {"inf", *range(11)}
    assert PAIR[3] == frozenset({1, 2, 3, 4, 5})
    assert PAIR[9] == frozenset()


def test_main():
    out = main()
    assert out["proved"]["p11_kvec_maxplus_exists"] is True
    assert out["proved"]["not_interpolant_not_isolated"] is True
    assert out["proved"]["no_deg2_square_predicate"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
