import numpy as np
import pytest

from conference_commuting_mate_no_go import (
    commutator_parity_audit,
    exhaustive_paley_order_six_audit,
    mod_two_skew_inverse_audit,
    nomenclature_audit,
    paley_symmetric_conference,
    skew_signing_from_mask,
    symmetric_conference_commuting_no_go,
)


@pytest.mark.parametrize("n", [2, 6, 10, 14, 18])
def test_all_orders_commuting_no_go_ledger(n):
    row = symmetric_conference_commuting_no_go(n)
    assert row["number_of_diagonal_summands_is_odd"]
    assert row["commutator_diagonal_is_2_mod_4"]
    assert row["commutator_diagonal_frobenius_lower_bound"] == 4 * n
    assert row["conference_eigenspace_dimensions_are_odd"]
    assert row["skew_signing_is_invertible_mod_two"]
    assert not row["commuting_skew_signing_exists"]


@pytest.mark.parametrize("n", [2, 4, 6, 8, 10, 14])
def test_mod_two_reduction_of_even_order_skew_signing_is_an_involution(n):
    row = mod_two_skew_inverse_audit(n)
    assert row["square_is_identity"]
    assert row["all_even_order_skew_signings_have_odd_determinant"]
    assert row["all_even_order_skew_signings_are_invertible_over_R"]


@pytest.mark.parametrize("q", [5, 13, 17])
def test_small_paley_symmetric_conference_matrices(q):
    A = paley_symmetric_conference(q)
    assert np.array_equal(A, A.T)
    assert np.array_equal(A @ A, q * np.eye(q + 1, dtype=np.int64))


@pytest.mark.parametrize("mask", [0, 1, 0x1234, (1 << 15) - 1])
def test_commutator_diagonal_parity_identity(mask):
    A = paley_symmetric_conference(5)
    R = skew_signing_from_mask(6, mask)
    row = commutator_parity_audit(A, R)
    assert row["diagonal_identity_exact"]
    assert row["diagonal_nonzero"]
    assert row["diagonal_mod_four"] == [2] * 6
    assert row["diagonal_frobenius_squared"] >= 24
    assert not row["commutes"]


def test_orthogonal_design_nomenclature_translation():
    A = paley_symmetric_conference(5)
    R = skew_signing_from_mask(6, 0x1234)
    row = nomenclature_audit(A, R)
    assert row["amicable_defect_equals_negative_anticommutator"]
    assert row["anti_amicable_defect_equals_negative_commutator"]
    assert row["amicable_means"] == "AR+RA=0"
    assert row["anti_amicable_means"] == "AR-RA=0"


def test_exhaustive_order_six_audit():
    row = exhaustive_paley_order_six_audit()
    assert row["all_checks"]
    assert row["skew_signings_checked"] == 32768
    assert row["commuting_signings"] == 0
    assert row["minimum_diagonal_frobenius_squared"] >= 24
    assert row["finite_audit_is_not_the_proof"]


@pytest.mark.parametrize("n", [1, 3, 4, 8, 12])
def test_conference_order_ledger_rejects_wrong_order_class(n):
    with pytest.raises(ValueError, match="2 modulo 4"):
        symmetric_conference_commuting_no_go(n)
