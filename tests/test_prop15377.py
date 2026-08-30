"""Tests for Prop 15.377 — all p=7 atoms constructed; QR0⊙NC ≠ Q_AP."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15355 import Q_AP_named
from e1_gmin_m4_prop15376 import JOINT_P7
from e1_gmin_m4_prop15377 import (
    Qmm_interp_p7p11,
    main,
    prove_A,
    prove_B,
    prove_open,
    qr0_nc_Q,
)


def test_all_eight_atoms_constructed():
    A = prove_A()
    assert A["proved"] is True
    assert A["all_eight"] is True
    assert "40/9|8/3|32/9" in A["ystar_nc"]
    assert "40/9|80/21|24/7" in A["ystar_nc"]
    assert "32/9|24/7|256/63" in A["ystar_nc"]
    assert "40/9|24/7|76/21" in A["qr0_nc"]
    assert A["sharp"] == ("8/3", "4", "4")
    assert set(JOINT_P7) == {
        ("40/9", "0", "0"),
        ("32/3", "0", "0"),
        ("8/3", "80/21", "24/7"),
        ("40/9", "24/7", "76/21"),
        ("32/9", "24/7", "256/63"),
        ("40/9", "8/3", "32/9"),
        ("40/9", "80/21", "24/7"),
        ("8/3", "4", "4"),
    }


def test_QR0NC_not_QAP_and_interp_dies():
    B = prove_B()
    assert B["proved"] is True
    Q7 = qr0_nc_Q(7)
    assert Q7["++"] == Q_AP_named(7) == Fraction(40, 9)
    assert Q7["--N1"] == Qmm_interp_p7p11(7)
    Q11 = qr0_nc_Q(11)
    assert Q11["++"] != Q_AP_named(11)
    assert Q11["++"] == Fraction(3976, 495)
    Q19 = qr0_nc_Q(19)
    assert Q19["++"] != Q_AP_named(19)
    assert Q19["--N1"] != Qmm_interp_p7p11(19)
    assert Q19["--N1"] == Fraction(256, 171)


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["ystar_p19"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["all_eight_atoms_constructed"] is True
    assert out["proved"]["QR0NC_not_QAP"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
