"""Tests for Prop 15.390 — m from circle-line disc; k not f(n,m)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15390 import (
    chi_p,
    m_from_disc,
    m_from_disc_flip,
    main,
    prove_A,
    prove_B,
    prove_open,
    r_dir_rows,
)


def test_m_from_disc():
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["n_ok"] == A["rows"]["5"]["n"]
    assert A["rows"]["17"]["n_ok"] == A["rows"]["17"]["n"]
    assert A["rows"]["5"]["n_flip"] < A["rows"]["5"]["n"]
    pack, recs = r_dir_rows(5)
    assert pack["t"] == 0
    assert all(r["m_pred"] == r["m_live"] for r in recs)
    assert any(r["m_flip"] != r["m_live"] for r in recs)


def test_k_not_function_of_nm():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_split"] >= 1
    assert "3,2" in B["splits"]
    assert set(B["splits"]["3,2"]) == {0, 2}


def test_chi_p_and_disc():
    assert chi_p(0, 5) == 0
    assert m_from_disc(0, 5) == 1
    assert m_from_disc(1, 5) == 2
    assert m_from_disc_flip(1, 5) == 0


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["m_from_disc"] is True
    assert out["proved"]["k_not_fn_nm"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
