"""Tests for Prop 15.549 — K_λ complete sum minus Gauss collisions."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15549 import (
    K_lam_coll_wrong_omit_aabb,
    K_lam_from_complete,
    K_lam_gauss_coll,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_G_chi_rxi_is_p_not_minus_p():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        assert abs(A["by_p"][str(p)]["G_rxi"] - p) < 1e-8
        assert abs(A["by_p"][str(p)]["G_mrxi"] - p) < 1e-8
        assert abs(A["by_p"][str(p)]["G_rxi"] + p) > 1.0


def test_K_lam_is_complete_minus_2Nq_qminus2():
    B = prove_B()
    assert B["proved"] is True
    rec = B["fold"]
    assert rec["match"] is True
    assert rec["drop_fails"] is True
    assert rec["omit_aabb_fails"] is True
    N = rec["N"]
    q = 25
    assert K_lam_gauss_coll(5, N) == N * q * (3 * q - 5) == 227500
    assert abs(rec["K_lam"] - 21550) < 1e-4
    assert abs(rec["K_all"] - rec["K_lam"] - rec["coll"]) < 1e-4
    assert abs(K_lam_from_complete(rec["K_all"], 5, N) - rec["K_lam"]) < 1e-4
    assert abs(rec["K_all"] - rec["K_lam"]) > 1.0
    wrong = rec["K_all"] - K_lam_coll_wrong_omit_aabb(5, N)
    assert abs(wrong - rec["K_lam"]) > 1.0
    assert K_lam_from_complete(rec["p7_Kall_from_B"], 7, 5726) == 2206274.0


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["K_all_named_in_p"] is False
    assert C["NUM_SUM_16pA_general"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
