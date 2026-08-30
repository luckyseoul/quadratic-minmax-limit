"""Tests for Prop 15.375 — Q--_form in the p≡3 floor window."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15355 import Q_QR0_p3_named
from e1_gmin_m4_prop15367 import A_num, k_named
from e1_gmin_m4_prop15371 import D_form
from e1_gmin_m4_prop15372 import Qpp_form
from e1_gmin_m4_prop15374 import LIVE_QMM, LIVE_QN, Qm_min_sm0, even_S_p3
from e1_gmin_m4_prop15375 import (
    N_num_lo,
    Qmm_form,
    Qn_form_hit,
    dC_lo_mm,
    dC_sph,
    dC_spo,
    dC_tsh,
    f_lo,
    f_sph,
    f_spo,
    f_tsh,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    t_smh,
)


def test_affine_vanish_on_minus():
    A = prove_A()
    assert A["proved"] is True
    assert A["AP_mm"] == "0"
    assert A["QR0_mm"] == "0"
    assert A["QR0_pp"] == str(Q_QR0_p3_named(7))
    assert Fraction(A["QR0_pp"]) != 0


def test_form_hits_live_dropk_dies():
    B = prove_B()
    assert B["proved"] is True
    assert Qmm_form(7) == LIVE_QMM[7] == Fraction(1376, 409)
    assert Qpp_form(7) == Fraction(1544, 409)
    assert Qpp_form(7) != LIVE_QMM[7]
    assert Qn_form_hit(7) == LIVE_QN[7] == Fraction(1496, 409)
    assert Qn_form_hit(7) != LIVE_QMM[7]
    assert A_num(3) - k_named(3) == 0


def test_Qmm_form_ge_Qlo_and_sm0():
    C = prove_C()
    assert C["proved"] is True
    assert dC_lo_mm(7) == 132
    assert N_num_lo(7) == -8784
    assert f_lo(7) == 10944
    assert N_num_lo(7) < 0
    for p in (7, 11, 19, 23):
        assert Qmm_form(p) >= Q_lo_named(p)
        assert Qmm_form(p) <= Qpp_form(p)
        assert Qm_min_sm0(p, Qpp_form(p)) <= Q_lo_named(p)
        assert even_S_p3(p, Qpp_form(p), Qmm_form(p))["sm0"] >= 6


def test_tsh_and_all_even_S():
    D = prove_D()
    assert D["proved"] is True
    assert dC_tsh(7) == 408
    assert dC_spo(7) == 348
    assert dC_sph(7) == 696
    assert f_tsh(7) == 228672
    assert f_spo(7) == 171936
    assert f_sph(7) == 537408
    p = 7
    Qf = Qpp_form(p)
    Qm = Qmm_form(p)
    assert Qf <= t_smh(p, Qf)
    S = even_S_p3(p, Qf, Qm)
    assert min(S[k] for k in S if k != "Qn") >= 6
    # drop-k still in the window (not a floor fail; it fails the live name)
    S_drop = even_S_p3(p, Qf, Qf)
    assert S_drop["smh"] >= 6


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["D_named_in_p"] is False
    assert E["Qmm_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["affine_Qmm_zero"] is True
    assert out["proved"]["form_hits_live_p7"] is True
    assert out["proved"]["Qmm_form_ge_Qlo"] is True
    assert out["proved"]["even_S_ge_6_at_form_pair"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["proved"]["Qmm_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
    assert D_form(7) == 409
