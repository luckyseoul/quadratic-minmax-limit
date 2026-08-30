"""Tests for Prop 15.447 — Q_τ as half-Gauss of split L through named C̄."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15447 import (
    SN_named,
    chsum,
    main,
    parts,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    recover_Q,
)


def test_SN_not_zero():
    assert SN_named(13) == -96
    assert SN_named(13) != 0
    F = _kernel_field(13)
    P = parts(F, 13)
    assert chsum(F, P["Sub"], P["N*"])[0] == -96
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["13"]["TN"] == 0


def test_NN_named():
    F = _kernel_field(5)
    P = parts(F, 5)
    assert chsum(F, P["N*"], P["N*"])[0] == 4
    B = prove_B()
    assert B["proved"] is True


def test_recover_and_drop_SN():
    rec = recover_Q(5, drop_ST=False)
    bad = recover_Q(5, drop_ST=True)
    assert rec["++"] == Fraction(48, 13)
    assert rec["N*"] == Fraction(32, 13)
    assert bad["++"] != Fraction(48, 13)
    C = prove_C()
    assert C["proved"] is True
    assert C["Qpp"] == "48/13"


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.447"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["SN_TN_named"] is True
    assert out["proved"]["NN_named"] is True
    assert out["proved"]["Q_from_split_L"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["C"]["Qpp"] == "48/13"
    assert out["algebra"]["C"]["Qpp_dropST"] != "48/13"
