import itertools
import json
from fractions import Fraction
from pathlib import Path

import pytest

from e1_gmin_m4_prop15764 import (
    abstract_method_barrier,
    balanced_double_star_anticommutator,
    deletion_score,
    level_five_degree_ledger,
    level_five_deeper_case_excluded_all_primes,
    minimal_gap4_shell_bridge_closed_general,
    no_bridge_size_floor,
    official_unit_entry_ledger,
    odd_small_bridge,
    paley_minus_phase_normalization,
    parity_bridge_ledger,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_odd_parity_bridge_is_an_exact_sign_sum_equivalence():
    # Under the four-gap floor, an odd H-row has score at least three.
    # Exhausting abstract sign rows here checks only the displayed identity;
    # the theorem itself is symbolic and does not rely on this finite audit.
    for h in (3, 5, 7, 9):
        for signs in itertools.product((-1, 1), repeat=h):
            total = sum(signs)
            if total < 3:
                continue
            has_level_two_deletion = any(
                deletion_score(total, sign) == 2 for sign in signs
            )
            assert has_level_two_deletion == (total == 3)
            if total == 3:
                assert all(
                    sign == 1
                    for sign in signs
                    if deletion_score(total, sign) == 2
                )
        row = parity_bridge_ledger(h)
        assert row["H_parity"] == "odd"
        assert row["deletion_parity"] == "even"
        assert row["residual_ii_level_two_possible"]


def test_even_parity_bridge_lands_at_type_one_not_residual_two():
    for h in (2, 4, 6, 8):
        for signs in itertools.product((-1, 1), repeat=h):
            total = sum(signs)
            if total < 2:
                continue
            has_deletion_at_most_two = any(
                deletion_score(total, sign) <= 2 for sign in signs
            )
            assert has_deletion_at_most_two == (total == 2)
            if has_deletion_at_most_two:
                assert any(deletion_score(total, sign) == 1 for sign in signs)
                assert all(
                    sign == 1
                    for sign in signs
                    if deletion_score(total, sign) == 1
                )
        row = parity_bridge_ledger(h)
        assert row["H_parity"] == "even"
        assert row["deletion_parity"] == "odd"
        assert not row["residual_ii_level_two_possible"]


def test_frame_averaging_and_bitight_boundaries_give_exact_open_floors():
    for p in (5, 7, 11, 13, 17):
        assert no_bridge_size_floor(p, True) == 5 * p + 2
        assert no_bridge_size_floor(p, False) == 4 * p + 2
        for h in (1, 3, 5 * p - 2, 5 * p):
            row = odd_small_bridge(p, h)
            assert row["frame_mean"] == str(Fraction(h, p))
            assert row["critical_level_two_deletion_forced"]


def test_critical_rows_enter_the_official_phase_normalized_units():
    phase = paley_minus_phase_normalization()
    assert phase["identity"] == "-C=D*P^T*C*P*D"
    assert phase["minus_shell_becomes_plus_shell"]
    assert phase["proved"]

    for p in (5, 7, 11, 13, 17):
        residual = official_unit_entry_ledger(p, True)
        assert residual["G_parity"] == "even"
        assert residual["critical_G_score"] == 2
        assert residual["active_edge_sign_on_every_G_level_two_row"] == 1
        assert residual["both_phase_G_shell_floor"] == 2
        assert residual["sharp_H_size_floor"] == 3 * p + 2
        assert residual["sharp_G_size_floor"] == 3 * p + 1
        assert residual["official_entry_proved"]

        type_i = official_unit_entry_ledger(p, False)
        assert type_i["G_parity"] == "odd"
        assert type_i["critical_G_score"] == 1
        assert type_i["active_edge_sign_on_every_G_level_one_row"] == 1
        assert type_i["sharp_H_size_floor"] == 2 * p + 2
        assert type_i["sharp_G_size_floor"] == 2 * p + 1
        assert type_i["official_entry_proved"]


def test_level_five_degree_congruence_tail_and_p7_are_empty():
    p7 = level_five_degree_ledger(7)
    assert p7["modulus"] == 24
    assert p7["common_residues"] == [11, 23]
    assert p7["minimum_possible_degree_sum"] == 550 > 70
    assert p7["arithmetic_empty"]

    for p in (11, 13, 17, 19):
        row = level_five_degree_ledger(p)
        assert row["modulus_exceeds_H_size"]
        assert 0 < Fraction(row["forced_regular_degree"]) < 1
        assert row["arithmetic_empty"]


def test_p5_level_five_profiles_and_anticommutator_are_exact():
    row = level_five_degree_ledger(5)
    assert row["modulus"] == 12
    assert row["common_residues"] == [1, 7]
    assert row["viable_common_residues"] == [1]
    assert sorted(row["profiles_as_counts_degree_25_13_1"]) == [
        [0, 2, 24],
        [1, 0, 25],
    ]
    anti = balanced_double_star_anticommutator()
    assert anti["centre_edge_forced"]
    assert anti["leaf_partition_sizes"] == [12, 12]
    assert anti["contradiction_for_all_signs"]
    assert all(entry["lhs"] == -entry["rhs"] for entry in anti["sign_audit"])
    assert level_five_deeper_case_excluded_all_primes()


def test_abstract_method_barrier_satisfies_every_claimed_scalar_identity():
    row = abstract_method_barrier()
    assert not row["paley_realizable_claimed"]
    assert row["shell_coordinate_mean"] == row["expected_frame_mean"] == "1/5"
    assert row["H_size"] == 5 * row["p"] == 25
    assert row["shell_H_score"] == 5
    assert row["deletion_shell_scores"] == [4, 6]
    assert row["deletion_shell_minimum"] == 4
    assert row["H_model_norm"] == row["Phi"] - 4
    assert row["one_deletion_model_norm"] == row["Phi"] - 2
    assert all(
        spike["score_at_H"] == row["Phi"] - 4
        for spike in row["spike_rows"]
    )
    assert all(
        witness["strictly_above_Phi_minus_4"]
        for witness in row["proper_subset_witnesses"].values()
    )
    assert row["H_inclusion_minimal_at_four_gap"]
    assert row["all_deletions_two_gap"]
    assert not row["critical_shell_deletion_exists"]
    assert row["method_barrier_proved"]


def test_parameter_guards():
    for p in (2, 3, 4, 9):
        with pytest.raises(ValueError):
            no_bridge_size_floor(p, True)
    with pytest.raises(ValueError):
        parity_bridge_ledger(0)
    with pytest.raises(ValueError):
        deletion_score(3, 0)


def test_checked_in_evidence_and_docs_preserve_the_open_scope():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15764.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["proved"]["parity_equivalence"]
    assert observed["proved"]["odd_H_at_most_5p_bridge"]
    assert not observed["proved"]["even_H_bridge_to_residual_ii"]
    assert not observed["proved"]["large_odd_H_bridge"]
    assert not observed["proved"]["residual_ii_closed"]
    assert not minimal_gap4_shell_bridge_closed_general()
    assert observed["L_status"] == "OPEN"

    for name in (
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "README.md",
        "solution.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
    ):
        text = (ROOT / name).read_text(encoding="utf-8", errors="replace")
        flat = " ".join(text.split()).lower()
        assert "15.764" in text, name
        assert "residual (ii)" in flat and "open" in flat, name
