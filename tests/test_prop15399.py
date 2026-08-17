"""Tests for Prop 15.399 — J++(ψ_2)=−4 all p; Gram cone misses floor."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named, Qn_from_Qpp
from e1_gmin_m4_prop15397 import J_type
from e1_gmin_m4_prop15398 import Q_eq_named, sbox_qn_and_four
from e1_gmin_m4_prop15399 import (
    fp_units,
    main,
    prove_A,
    prove_B,
    prove_open,
    psd_cone_forces_qn,
    psd_hi_named,
    psd_lo_named,
    psi2_sum,
    qlo_minus_psdlo_num,
    torus_U,
    usq_off,
)


def test_Jpp_all_p_via_complete_sums():
    A = prove_A()
    assert A["proved"] is True
    F = _kernel_field(5)
    assert abs(J_type(F, "++", 2) + 4) < 1e-8
    fp = fp_units(5)
    assert abs(psi2_sum(F, fp)) < 1e-8
    off = [r for r in fp if r not in (1, F["neg1"])]
    assert abs(psi2_sum(F, off) + 2) < 1e-8
    # fail-when-wrong: the full F_p^× sum is 0, not −2
    assert abs(psi2_sum(F, fp) + 2) > 0.5


def test_p3_N1_is_Usq_not_whole_U():
    F = _kernel_field(7)
    Uoff = [r for r in torus_U(F) if r not in (1, F["neg1"])]
    assert set(usq_off(F)) != set(Uoff)
    assert len(usq_off(F)) == (7 - 3) // 2


def test_I_psd_strictly_contains_floor():
    B = prove_B()
    assert B["proved"] is True
    assert psd_lo_named(5) == 2
    assert Q_lo_named(5) == Fraction(11, 4)
    assert psd_lo_named(5) < Q_lo_named(5)
    assert Q_eq_named(5) < Qpp_floor_ub(5) < psd_hi_named(5)
    assert qlo_minus_psdlo_num(5) == 12
    assert (2 - Qn_from_Qpp(5, Fraction(2))) == -3


def test_cone_does_not_force_qn_or_phi():
    C = prove_open()
    assert C["proved"] is False
    assert psd_cone_forces_qn() is False
    assert sbox_qn_and_four() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Jpp2_eq_m4_all_p"] is True
    assert out["proved"]["I_psd_contains_floor"] is True
    assert out["proved"]["psd_cone_forces_qn"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
