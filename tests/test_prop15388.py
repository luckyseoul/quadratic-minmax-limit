"""Tests for Prop 15.388 — T⊙NC Hoffman meet/delta named."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15388 import (
    O_C_from_delta,
    T_cap_S,
    circle_meets,
    delta_p1,
    delta_p3_t0,
    live_pack,
    main,
    meet_even,
    meet_tangent,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_meet_tangent_not_even():
    A = prove_A()
    assert A["proved"] is True
    assert meet_tangent(7) == (0, 0, 1, 1, 2, 2, 2)
    assert meet_even(7) == (0, 0, 0, 2, 2, 2, 2)
    assert meet_tangent(7) != meet_even(7)
    pack = live_pack(7)
    meets = circle_meets(7, pack)
    assert all(m == meet_tangent(7) for m in meets)
    assert not any(m == meet_even(7) for m in meets)


def test_half_circle_and_delta():
    B = prove_B()
    assert B["proved"] is True
    assert T_cap_S(live_pack(5)) == 3
    assert T_cap_S(live_pack(7)) == 4
    assert T_cap_S(live_pack(7)) != 8
    assert delta_p1(5) == (-2, -1, 0, 1, 2)
    assert delta_p3_t0(7) == (-2, -2, 0, 0, 1, 1, 2)


def test_O_C_named():
    C = prove_C()
    assert C["proved"] is True
    assert O_C_from_delta(5, delta_p1(5)) == (0, 1, 2, 3, 4)
    assert O_C_from_delta(7, delta_p3_t0(7)) == (1, 1, 3, 3, 4, 4, 5)
    assert O_C_from_delta(7, delta_p3_t0(7)) != (0, 1, 2, 3, 4, 5, 6)


def test_phi_F_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["meet_tangent"] is True
    assert out["proved"]["delta_named"] is True
    assert out["proved"]["O_C_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
