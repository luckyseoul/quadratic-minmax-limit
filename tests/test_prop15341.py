"""Tests for Prop 15.341 — named S^{(k)} pair-selection at p=11."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15341 import (
    build_named,
    energy,
    generators_of_squares,
    main,
    named_g,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_named_rebuild_is_maxplus():
    A = prove_A()
    assert A["proved"] is True
    assert A["energy"] < 1e-10
    assert A["nD"] == 55
    assert A["eq_15340"] is True
    assert A["g"] == 3
    y = build_named(11)
    assert energy(y, 11) < 1e-10


def test_only_named_generator():
    B = prove_B()
    assert B["proved"] is True
    assert B["named_g"] == named_g(11)
    assert set(B["gens"]) == {3, 4, 5, 9}
    assert B["by_g"]["3"]["energy"] < 1e-10
    for g in (4, 5, 9):
        assert B["by_g"][str(g)]["energy"] > 1.0
        assert B["by_g"][str(g)]["nD"] == 55


def test_mutations_break():
    C = prove_C()
    assert C["proved"] is True
    assert C["e_named"] < 1e-10
    assert C["e_S0_pm2"] > 1.0
    assert C["e_A_id"] > 1.0


def test_generators_helper():
    assert generators_of_squares(11) == [3, 4, 5, 9]


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["named_rebuild_p11"] is True
    assert out["proved"]["only_named_g"] is True
    assert out["proved"]["mutations_break"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
