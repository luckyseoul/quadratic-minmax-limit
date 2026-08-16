"""Tests for Prop 15.278: Aut·F spans Z ⇒ spec(Φ)=spec(Φ|_F); ≥6 still OPEN."""
from __future__ import annotations

from e1_gmin_m4_prop15159 import dim_Z
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15277 import lambda_min_ge_6_proved_general
from e1_gmin_m4_prop15270 import psl_span_F_eq_Wpp0, theorem_aut_schur
from e1_gmin_m4_prop15278 import (
    aut_span_F_equals_Z_proved_general,
    certify_PF_iota_injective,
    certify_aut_span_F,
    certify_elliptic_trPC_formula,
    certify_one_elliptic_chiV_zero,
    certify_quadratic_char_sum,
    certify_sl2_disc_identity,
    chi_V_from_signed_trace,
    chi_Z_elliptic_closed,
    chi_Z_from_weil,
    chi_Z_unipotent_minus,
    chi_Z_unipotent_plus,
    chi_minus_one_Fq,
    dim_F,
    dim_formulas_ok,
    discrete_series_identity_mass,
    elliptic_torus_order,
    lambda_min_via_F_proved_general,
    phi_F_ge_6_proved_general,
    sl2_quadratic_disc_gap,
    theorem_PF_iota_injective_criterion,
    theorem_aut_span_F_equals_Z_certified,
    theorem_aut_span_F_equals_Z_general,
    theorem_aut_span_implies_spec_equality,
    theorem_jacquet_aut_span_F_equals_Z,
    theorem_no_cuspidal_in_Z,
    theorem_phi_F_meets_global_census,
    theorem_signed_inversion_preserves_C,
    theorem_squares_transitive_good_bad,
    theorem_translate_span_is_isotypic_sum,
    theorem_unipotent_cancels_discrete_series,
    theorem_weil_zero_on_elliptic,
    trPC_elliptic_closed,
    unipotent_chiZ_mass,
    unipotent_class_size,
)


def test_dim_F_formula():
    for p in (5, 7, 11, 13, 17, 19):
        assert dim_F(p) == (p * p - 5) // 4
        assert 4 * dim_F(p) + 5 == p * p
        assert dim_F(p) < dim_Z(p)


def test_signed_inversion_preserves_C_unsigned_fails():
    r = theorem_signed_inversion_preserves_C([5, 7])
    assert r["proved"] is True
    for p in (5, 7):
        assert r["by_p"][str(p)]["signed_ok"] is True
        assert r["by_p"][str(p)]["unsigned_fails"] is True


def test_spec_equality_is_algebra():
    r = theorem_aut_span_implies_spec_equality()
    assert r["proved"] is True
    assert r["aut_schur"] is False


def test_translate_span_criterion_is_algebra():
    r = theorem_translate_span_is_isotypic_sum()
    assert r["proved"] is True
    # dim F >= dim good >= dim bad for every prime p>=5
    for p in (5, 7, 11, 13, 17, 19, 23):
        dF = dim_F(p)
        dgood = (p * p - 1) // 8
        dbad = (p * p - 9) // 8
        assert dF >= dgood >= dbad
        assert dF == (p * p - 5) // 4


def test_aut_span_F_p5_is_full_Z():
    """Live SVD: Aut·F has rank dim Z at p=5.  Fails if inversion/span is wrong."""
    c = certify_aut_span_F(5)
    assert c["dimZ"] == dim_Z(5) == 65
    assert c["dimF"] == dim_F(5) == 5
    assert c["aut_span_rank"] == 65
    assert c["full_span"] is True
    t = theorem_aut_span_F_equals_Z_certified()
    assert t["proved_census"] is True
    # census object is not the general proof; Jacquet is
    assert t["proved_general"] is False


def test_unipotent_mass_cancels_discrete_series_identity():
    """Fraction: dim Z·(q−1) = ∑_unip χ_Z for every prime p≥5."""
    from fractions import Fraction

    r = theorem_unipotent_cancels_discrete_series()
    assert r["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert discrete_series_identity_mass(p) == unipotent_chiZ_mass(p)
        assert chi_Z_unipotent_plus(p) + chi_Z_unipotent_minus(p) == Fraction(
            p * p - 5, 4
        )
        # Weil → χ_Z: Tr(PC)=±p², fix=1, g² same class
        xp = chi_V_from_signed_trace(1, p * p, p)
        xm = chi_V_from_signed_trace(1, -(p * p), p)
        assert xp == Fraction(1 + p, 2)
        assert xm == Fraction(1 - p, 2)
        assert chi_Z_from_weil(xp, xp, 1) == chi_Z_unipotent_plus(p)
        assert chi_Z_from_weil(xm, xm, 1) == chi_Z_unipotent_minus(p)
        assert unipotent_class_size(p) == Fraction(p**4 - 1, 2)
    # pinned p=5,7 (wrong closed forms fail)
    assert chi_Z_unipotent_plus(5) == 5
    assert chi_Z_unipotent_minus(5) == 0
    assert chi_Z_unipotent_plus(7) == 9
    assert chi_Z_unipotent_minus(7) == 2
    assert chi_Z_unipotent_plus(5) != 0
    assert discrete_series_identity_mass(5) == 1560
    assert unipotent_chiZ_mass(5) == 1560


def test_elliptic_trPC_and_chi_Z_vanish():
    """Gauss-sum identities: Tr(PC)=1−χ(−1)=0, χ_Z=0 on elliptics."""
    from fractions import Fraction

    g = certify_quadratic_char_sum()
    assert g["proved"] is True
    assert g["by_p"]["5"]["sum_chi"] == 0
    assert g["by_p"]["5"]["sum_off_1"] == -1
    d = certify_sl2_disc_identity()
    assert d["proved"] is True
    # polynomial: Δ(Q)−(tr²−4)=4(1−det), zero on SL(2)
    dQ, dC, gap = sl2_quadratic_disc_gap(2, 1, 1, 1)  # det=1
    assert dQ == dC
    assert gap == 0
    dQ, dC, gap = sl2_quadratic_disc_gap(2, 1, 1, 2)  # det=3
    assert dQ - dC == gap == 4 * (1 - 3)
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert chi_minus_one_Fq(p) == 1
        assert elliptic_torus_order(p) % 2 == 1
        assert trPC_elliptic_closed(p) == 0
        assert chi_Z_elliptic_closed(p) == 0
        assert chi_V_from_signed_trace(0, 0, p) == 0
        # wrong: unsigned χ(−1)=−1 would give Tr(PC)=2
        assert trPC_elliptic_closed(p) != Fraction(2)
    live = certify_one_elliptic_chiV_zero(5)
    assert live["ok"] is True
    assert live["elliptic_no_fix"] is True
    assert live["chiV_zero"] is True
    tr = certify_elliptic_trPC_formula(5, 4)
    assert tr["ok"] is True
    w = theorem_weil_zero_on_elliptic()
    assert w["proved"] is True
    assert w["gauss"] is True
    assert w["disc"] is True


def test_squares_transitive_good_bad():
    r = theorem_squares_transitive_good_bad()
    assert r["proved"] is True
    assert r["by_p"]["5"]["n_good"] == r["by_p"]["5"]["n_bad"] == 12
    assert r["by_p"]["7"]["n_good"] == 24


def test_no_cuspidal_and_jacquet_give_aut_span_general():
    """Aut·F=Z for all p≥5.  k=3 Singer span and Aut-Schur stay False."""
    n = theorem_no_cuspidal_in_Z()
    assert n["proved"] is True
    assert n["does_not_imply_k3_span"] is True
    j = theorem_jacquet_aut_span_F_equals_Z()
    assert j["proved"] is True
    assert j["aut_schur"] is False
    assert j["k3_span"] is False
    g = theorem_aut_span_F_equals_Z_general()
    assert g["proved"] is True
    assert dim_formulas_ok() is True
    # flag requires theorem + signed invert + identities + p=5 SVD
    assert aut_span_F_equals_Z_proved_general() is True
    assert theorem_aut_schur()["proved"] is False
    assert psl_span_F_eq_Wpp0() is False
    # Φ|_F floor still open, so via-F λ_min stays False
    assert phi_F_ge_6_proved_general() is False


def test_phi_F_census_not_general():
    r = theorem_phi_F_meets_global_census()
    assert r["proved_census"] is True
    assert r["proved_general"] is False
    assert r["p5_min_ge_6"] is True
    assert r["p7_min_ge_6"] is True


def test_PF_iota_criterion_and_p5_injective():
    """Aut·F=Z iff P_F∘ι injective on complex Z_μ; live at p=5."""
    c = theorem_PF_iota_injective_criterion()
    assert c["proved"] is True
    r = certify_PF_iota_injective(5)
    assert r["all_injective"] is True
    assert r["proved_general"] is False
    dims = sorted({row["dim"] for row in r["by_mu"].values()})
    assert dims == [2, 3]  # bad / good
    assert all(row["injective"] for row in r["by_mu"].values())


def test_lambda_min_via_F_and_residual_still_open():
    assert lambda_min_via_F_proved_general() is False
    assert lambda_min_ge_6_proved_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() in (True, False)
