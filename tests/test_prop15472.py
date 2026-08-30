"""Tests for Prop 15.472 — Z_ψ not constant; GJ miss Λ_1."""
from __future__ import annotations

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15468 import S_chi_named
from e1_gmin_m4_prop15471 import LIVE_LAM1_SQ, Lambda_tables, mean_on_chi
from e1_gmin_m4_prop15472 import (
    Z_psi,
    gauss_psi,
    jacobi_pp,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_Zpsi_not_constant():
    A = prove_A()
    assert A["proved"] is True
    Z = Z_psi(5)
    assert abs(Z.mean()) < 1e-8
    assert np.std(np.abs(Z)) > 1.0
    assert np.max(np.abs(Z - S_chi_named(5))) > 1.0
    Z7 = Z_psi(7)
    assert abs(Z7.mean()) < 1e-8
    assert np.std(np.abs(Z7)) > 1.0


def test_Lam1_not_GJ():
    B = prove_B()
    assert B["proved"] is True
    F, _, Lam = Lambda_tables(5)
    lam = mean_on_chi(F, Lam[1], 1)
    assert abs(lam - LIVE_LAM1_SQ[5]) < 1e-6
    assert abs(lam - gauss_psi(5)) > 1.0
    assert abs(lam - jacobi_pp(5)) > 1.0
    F7, _, Lam7 = Lambda_tables(7)
    lam7 = mean_on_chi(F7, Lam7[1], 1)
    assert abs(lam7 - gauss_psi(7)) > 1.0


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.472"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["Zpsi_not_constant"] is True
    assert out["proved"]["Lam1_not_GJ"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
