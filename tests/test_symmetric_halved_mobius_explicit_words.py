from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_halved_mobius_explicit_words import (
    exact_small_explicit_word_replay,
    explicit_word_intersection_theorem,
    theorem_record,
)


@pytest.mark.parametrize("p", [3, 5, 7, 31])
def test_symbolic_intersection_theorem(p: int) -> None:
    record = explicit_word_intersection_theorem(p)
    assert record["proved"]
    assert record["all_odd_prime_symbolic_proof"]
    assert record["vertical_fibre"]["one_half_intersection_at_most"] == 1
    assert record["vertical_fibre"]["cannot_be_contained"]
    assert record["scalar_graphs"]["identity_graph_one_half_intersection"] == 1
    assert record["scalar_graphs"]["nonidentity_graph_one_half_intersection"] == 0
    assert record["scalar_graphs"]["all_scalar_graphs_cannot_be_contained"]
    assert record["Mobius_half_count"] < record["delta_size"]


@pytest.mark.parametrize("p,expected", [(3, 24), (7, 336)])
def test_tiny_replay(p: int, expected: int) -> None:
    replay = exact_small_explicit_word_replay(p)
    assert replay["proved"]
    assert replay["localized_halves_checked"] == expected
    assert replay["difference_class_map_injective"]
    assert replay["vertical_hit_criterion_exact"]
    assert replay["one_identity_scalar_intersection_only"]
    assert "not theorem evidence" in replay["role"]


def test_replay_refuses_prime_census() -> None:
    with pytest.raises(ValueError, match="limited to p=3,7"):
        exact_small_explicit_word_replay(5)


@pytest.mark.parametrize("bad_p", [False, 1, 2, 9])
def test_symbolic_theorem_requires_an_odd_prime(bad_p: int) -> None:
    with pytest.raises(ValueError, match="odd prime"):
        explicit_word_intersection_theorem(bad_p)


def test_top_level_record_keeps_scope_open() -> None:
    record = theorem_record(31)
    assert record["proved_all_claimed_statements"]
    assert record["proved"]["vertical_fibre_containment_excluded"]
    assert record["proved"]["scalar_graph_containment_excluded"]
    assert not record["proved"]["all_low_weight_row_words_excluded"]
    assert not record["proved"]["structured_punctured_surjectivity"]
    assert not record["proved"]["residual_ii_closed"]
