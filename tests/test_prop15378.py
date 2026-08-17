"""Tests for Prop 15.378 — QR0⊙NC interpolants die at p=19."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15355 import Q_AP_named
from e1_gmin_m4_prop15367 import k_named
from e1_gmin_m4_prop15377 import Qmm_interp_p7p11, qr0_nc_Q
from e1_gmin_m4_prop15378 import (
    Qn_2pt_fit,
    Qn_3pt_fit,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_Qmm_interp_dies_at_p19():
    A = prove_A()
    assert A["proved"] is True
    assert Qmm_interp_p7p11(19) == Fraction(144, 95)
    assert qr0_nc_Q(19)["--N1"] == Fraction(256, 171)
    assert qr0_nc_Q(19)["--N1"] != Qmm_interp_p7p11(19)
    assert qr0_nc_Q(11)["++"] != Q_AP_named(11)


def test_Qn_3pt_is_interpolant_2pt_dies():
    B = prove_B()
    assert B["proved"] is True
    for p, want in ((7, 76), (11, 144), (19, 272)):
        assert qr0_nc_Q(p)["N*"] * k_named(p) == want
        assert Qn_3pt_fit(p) == qr0_nc_Q(p)["N*"]
    assert Qn_2pt_fit(19) == Fraction(280, 171)
    assert Qn_2pt_fit(19) != Fraction(272, 171)


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Qmm_interp_dies_p19"] is True
    assert out["proved"]["Qn_3pt_is_interpolant"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
