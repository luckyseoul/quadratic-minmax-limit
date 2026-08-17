"""Tests for Prop 15.480 — quadratic-χ kernels are Kloosterman."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15480 import (
    S_exact,
    S_named,
    S_named_wrong_sign,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_S_is_kl():
    A = prove_A()
    assert A["proved"] is True
    F = _kernel_field(5)
    err = abs(S_exact(F) - S_named(F)).max()
    assert err < 1e-10
    err_w = abs(S_exact(F)[1:] - S_named_wrong_sign(F)[1:]).max()
    assert err_w > 1.0


def test_complete_square():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["err_rm"] < 1e-10
    assert B["rows"]["5"]["err_drop_phase"] > 0.5
    assert B["rows"]["7"]["err_rp"] < 1e-10


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.480"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["S_is_kloosterman"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
