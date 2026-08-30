"""Tests for Prop 15.525 — type-index Gram misses δ; no catalog D=N(a+bi)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15525 import (
    D_as_half_sumsq,
    catalog_hits_D,
    dim_Vplus,
    laws_through,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_named_gram_matches_live():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["named_live_err"] < 1e-6
    assert A["rows"]["5"]["drop_minus_err"] > 1.0
    assert A["rows"]["7"]["named_live_err"] < 1e-6


def test_live_delta_not_eigenvector():
    B = prove_B()
    assert B["proved"] is True
    assert B["p5"]["G"] == [[72, -48], [-48, 64]]
    assert B["p5"]["wedge"] == 17792
    assert B["p5"]["wedge"] != 0
    assert B["p7"]["svals"][1] > 1.0


def test_no_catalog_plaw_for_D():
    C = prove_C()
    assert C["proved"] is True
    assert catalog_hits_D() == []
    assert laws_through(2, 3) == ["(p-1)/2"]
    assert laws_through(3, 20) == []
    assert D_as_half_sumsq(7) == dim_Vplus(7) == 25
    assert D_as_half_sumsq(7) != D_named(7)
    assert D_as_half_sumsq(5) == D_named(5)


def test_jacobi_norm_is_p_not_D():
    D = prove_D()
    assert D["proved"] is True
    assert D["Jchi_k2_N"]["5"] == 25
    assert D["Jchi_k2_N"]["7"] == 49
    assert D_named(7) != 49


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert E["pairing_forced"] is False
    assert E["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["D"]["proved"] is True
    assert out["E"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
