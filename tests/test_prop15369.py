"""Tests for Prop 15.369 — affine S=T only; deg-2 dies at p=7."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15368 import dim_Vp, load_cat
from e1_gmin_m4_prop15369 import main, prove_A, prove_B, prove_C, prove_open


def test_affine_S_is_T():
    A = prove_A()
    assert A["proved"] is True
    assert A["p5_S_x_is_T"] is True
    assert A["nT"] == 2040


def test_deg2_constancy_dies():
    B = prove_B()
    assert B["proved"] is True
    assert B["p5_x2_eq_x2x"] is True
    assert B["p5_S"] == 40
    assert B["p7_x2_occ0"] != B["p7_x2x_occ0"]
    cat = load_cat()
    assert cat["by_alpha"]["x2"]["occ0"] != cat["by_alpha"]["x2+x"]["occ0"]


def test_D_not_named_linear():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_hits"] == 0
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
    assert out["proved"]["affine_S_is_T"] is True
    assert out["proved"]["deg2_constancy_dies_p7"] is True
    assert out["proved"]["D_not_3term_named"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
