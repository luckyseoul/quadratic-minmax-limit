"""Tests for Prop 15.473 — Z_ψ DFT identity."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15473 import dft_pack, main, prove_A, prove_B, prove_C, prove_open


def test_dft_identity():
    A = prove_A()
    assert A["proved"] is True
    rec = dft_pack(5)
    assert rec["err"] < 1e-10
    assert rec["err_chiG"] > 1.0
    assert dft_pack(7)["err"] < 1e-10


def test_drop_psi():
    B = prove_B()
    assert B["proved"] is True
    rec = dft_pack(5)
    assert rec["err_drop"] > 1.0
    assert abs(rec["W_drop_mean"] - rec["k_qmink"]) < 1e-6


def test_ratio_not_const():
    C = prove_C()
    assert C["proved"] is True
    r5 = dft_pack(5)["EZ2"] / float(LIVE_QPP[5])
    r7 = dft_pack(7)["EZ2"] / float(LIVE_QPP[7])
    assert abs(r5 - 68.75) < 1e-6
    assert abs(r7 - 68.75) > 1.0


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.473"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["dft_identity"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
