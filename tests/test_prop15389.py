"""Tests for Prop 15.389 — T⊙NC R-pairing (affine n + polarised k)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15389 import (
    ab_named,
    k_polarised,
    main,
    n_affine,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    tau_old,
)


def test_affine_ab():
    A = prove_A()
    assert A["proved"] is True
    assert (1, 0) not in set(ab_named(13))
    assert set(ab_named(5)) == {(4, 0), (3, 0), (3, 4)}
    assert set(ab_named(17)) == {(16, 0), (9, 0), (9, 16)}
    assert n_affine(1, 12, 0, 13) == 12


def test_polarise():
    B = prove_B()
    assert B["proved"] is True
    assert B["p5_13"] is True
    assert B["p17"] is False
    ns = [0, 4, 3, 2, 1]
    ms = [0, 2, 1, 1, 2]
    assert k_polarised(ns, ms, 5, invert=False) != k_polarised(ns, ms, 5, invert=True)


def test_tau_dies_at_17():
    C = prove_C()
    assert C["proved"] is True
    assert C["old_tau_dies_at_17"] is True
    assert C["p17_t"] == 0
    assert tau_old(9, 17) == 4
    tans = [tuple(sorted(t)) for t in C["p17_tangents"]]
    assert (6, 11) in tans
    # generic dirs are ±6, not the old ±4 guess


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
    assert out["proved"]["affine_n"] is True
    assert out["proved"]["polarised_k"] is True
    assert out["proved"]["tau_law_dies_p17"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
