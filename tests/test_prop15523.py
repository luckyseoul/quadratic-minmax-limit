"""Tests for Prop 15.523 — κ-contraction is 3p χ_p(−1)χ(ξ)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15523 import (
    G_chi,
    chi_p_minus_1,
    dropped_one_pairing,
    kappa_contraction,
    main,
    pairing_one,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_G_chi_sign():
    A = prove_A()
    assert A["proved"] is True
    assert G_chi(5) == -5
    assert G_chi(7) == 7
    assert chi_p_minus_1(5) == 1
    assert chi_p_minus_1(7) == -1
    # fail-when-wrong: G≡p
    assert G_chi(5) != 5


def test_three_pairings_vs_drop_one():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7, 11, 13):
        for cxi in (1, -1):
            full = kappa_contraction(p, cxi)
            drop = dropped_one_pairing(p, cxi)
            want = 3 * p * chi_p_minus_1(p) * cxi
            assert full == want
            assert drop == p * chi_p_minus_1(p) * cxi
            assert drop == pairing_one(p, cxi)
            # fail-when-wrong: drop one pairing
            assert full != drop
            assert full == 3 * drop


def test_field_sum_matches_three_pairings():
    C = prove_C()
    assert C["proved"] is True
    for p in (5, 7, 11, 13):
        for cxi in (1, -1):
            row = C["rows"][str(p)][str(cxi)]
            assert row["ok"] is True
            assert abs(row["live"] - row["want"]) < 1e-6
            assert abs(row["live"] - row["drop"]) > 1.0


def test_type_I_flags_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["type_I_multilevel_bad_case_ND_closed"] is False
    assert D["type_I_aut_e_3AB_positive_general"] is False
    assert D["Delta_conn_named"] is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["G_chi"] is True
    assert out["proved"]["three_pairings"] is True
    assert out["proved"]["field_sum"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is False
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
