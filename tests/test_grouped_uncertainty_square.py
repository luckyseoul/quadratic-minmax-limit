from __future__ import annotations

import pytest

from e1_gmin_m4_grouped_uncertainty_square import (
    canonical_square_remainder,
    even_support_grouped_uncertainty_theorem,
    exact_remainder_replay,
    grouped_uncertainty_theorem,
    theorem_record,
)


def test_symbolic_theorem_covers_all_odd_primes_and_both_support_parities():
    for p in (3, 5, 7, 11, 31, 43):
        out = grouped_uncertainty_theorem(p)
        assert out["proved"] is True
        assert out["prime_range"] == "every odd prime"
        assert out["p_congruent_3_mod_4_included"] is True
        assert out["grouped_inequality"].endswith(">=p+1")
        assert even_support_grouped_uncertainty_theorem(p)["proved"] is True


def test_square_remainder_is_nonzero_and_has_the_claimed_degree_bounds():
    supports = (
        (5, ((1, 0), (0, 1))),
        (7, ((1, 0), (0, 1), (1, 1), (1, 2))),
        (
            31,
            ((1, 11), (1, 19), (8, 11), (8, 12)),
        ),
    )
    for p, support in supports:
        out = canonical_square_remainder(p, support)
        assert out["proved"] is True
        assert out["nonzero_remainder_indices"]
        assert out["largest_remainder_degree"] == 2 * len(support)
        for index in out["nonzero_remainder_indices"]:
            assert out["homogeneous_degree_by_index"][index] == 2 * index
            assert index <= len(support)


def test_every_silent_direction_is_a_double_projective_zero():
    supports = (
        (5, ((1, 0), (0, 1))),
        (7, ((2, 1), (0, 1), (2, 2), (1, 5))),
        (31, ((1, 11), (1, 19), (8, 11), (8, 12))),
        (31, ((1, 0), (2, 0), (0, 1), (0, 2))),
    )
    for p, support in supports:
        out = exact_remainder_replay(p, support)
        assert out["proved"] is True
        assert out["grouped_bound_holds"] is True
        assert out["role"] == "fail-when-wrong replay, not theorem evidence"
        for orders in out["zero_orders"].values():
            assert all(order >= 2 for order in orders.values())


def test_record_does_not_promote_row_code_or_puncture_consequences():
    out = theorem_record(31)
    assert out["proved"]["grouped_uncertainty_all_supports"] is True
    assert out["proved"]["grouped_uncertainty_all_odd_primes"] is True
    assert out["proved"]["row_code_minimum_distance"] is False
    assert out["proved"]["minimum_word_classification"] is False
    assert out["proved"]["structured_mobius_puncture"] is False
    assert out["proved"]["symmetric_boolean_completion"] is False
    assert out["proved"]["residual_ii_closed"] is False


def test_invalid_prime_support_parity_and_duplicate_classes_fail_closed():
    with pytest.raises(ValueError, match="odd prime"):
        grouped_uncertainty_theorem(9)
    with pytest.raises(ValueError, match="even support"):
        canonical_square_remainder(7, ((1, 0), (0, 1), (1, 1)))
    with pytest.raises(ValueError, match="distinct modulo sign"):
        canonical_square_remainder(7, ((1, 2), (6, 5)))
