"""Tests for Prop 15.368 — p=7 S(α); S2=40 and 3-class die."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15368 import (
    T7,
    dim_Vp,
    load_cat,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_S2_S3_die_at_p7():
    A = prove_A()
    assert A["proved"] is True
    assert A["p7_x2_occ0"] == 179539290
    assert A["p7_S2_in_Q"] is False
    assert A["p7_S3_in_Q"] is False
    assert A["p7_x3_palindromic"] is True
    cat = load_cat()
    assert cat["nT"] == T7 == 68938800


def test_deg4_not_one_class():
    B = prove_B()
    assert B["proved"] is True
    assert B["x2_eq_x4"] is True
    assert B["x2_palindromic"] is False
    assert B["n_deg4_occ0"] >= 5
    cat = load_cat()
    assert cat["by_alpha"]["x4"]["occ0"] != cat["by_alpha"]["x4+x"]["occ0"]


def test_D_dim_dies():
    C = prove_C()
    assert C["proved"] is True
    assert dim_Vp(5) == 13
    assert dim_Vp(7) == 25
    assert dim_Vp(7) != 409


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["S2_S3_p5_die_p7"] is True
    assert out["proved"]["deg4_not_one_class_p7"] is True
    assert out["proved"]["D_dimVp_dies_p7"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
