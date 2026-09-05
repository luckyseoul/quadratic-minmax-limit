"""New mass-capacity theorem: exact scope, two carries, and failed premises."""
import json
from pathlib import Path

import pytest

import e1_gmin_m4_prop15774 as proof

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", params=[29, 31, 37, 43])
def p(request):
    return request.param


def test_quotient_floor_is_exact_for_union_mass_relaxation(p):
    m, threshold = (p + 1) // 2, 2 * p - 10
    small = {0, p - 3, p - 1, p + 1}
    for u in range(m):
        # Directly derive the first admissible mass from the four-point set.
        first = next(k for k in range(3)
                     if 2 * u + 2 * m * k in small or 2 * u + 2 * m * k >= threshold)
        assert proof.quotient_floor(p, u) == first
    assert proof.quotient_floor(p, m - 7) == 2
    assert proof.quotient_floor(p, m - 6) == 1  # strict endpoint is retained


@pytest.mark.parametrize("r", [3, 4, 5])
def test_capacity_both_types_and_all_earlier_t(p, r):
    record = proof.capacity_exclusion(p, r)
    m = record["m"]
    expected_H = r * p + 2 * (p + (r + 1) // 2 if p >= 37 else p - 6)
    assert record["maximum_H"] == expected_H
    assert record["guaranteed_isolated_vertices"] == p * p + 1 - 2 * expected_H > 0
    assert record["collisions"] == []
    assert record["both_types_use_union_spectrum"]
    assert record["even_r_phase_one_assumed"] is False
    for t in (0, 1, record["maximum_t"] - 1, record["maximum_t"]):
        residues = {u for u in range(m) if u + m * proof.quotient_floor(p, u) <= t}
        assert 0 in residues
        assert not any((r - u) % m in residues for u in residues)


@pytest.mark.parametrize("r", [3, 4, 5])
def test_first_uncovered_scalar_profile_not_a_graph(r):
    p, m = 37, 19
    record = proof.first_scalar_survivor(p, r)
    assert record["H"] == proof.capacity_exclusion(p, r)["maximum_H"] + 2
    assert [row["u"] for row in record["types"]] == [r // 2, (r + 1) // 2]
    assert record["t"] == 2 * m + (r + 1) // 2
    assert sum(sum(row["parallel_counts"]) for row in record["types"]) == record["H"]
    assert record["graph_realization_claimed"] is False
    assert record["full_local_row_realization_claimed"] is False
    for row in record["types"]:
        assert sum(row["quotients"]) == record["t"] - row["u"]
        assert set(row["quotients"]) <= {2, 3}
        assert min(row["parallel_counts"]) == r + 2
        assert min(row["masses"]) >= 2 * p - 10


@pytest.mark.parametrize("s", [1, 2])
def test_two_residual_layers_recompute_all_budgets_and_sign(p, s):
    rec = proof.residual_two_layer_exclusion(p, s)
    q, m = (p - 1) // 2, (p + 1) // 2
    assert rec["H"] == 5 * p + 2 * s
    assert rec["k"] == rec["H"] - 1
    assert rec["t"] == q + s
    assert rec["hard_sign_relative_to_transported_c_H"] == (-1) ** (q + s)
    cases = rec["residue_cases_below_q"]
    assert [row["u"] for row in cases] == list(range(q))
    for row in cases:
        u = row["u"]
        assert row["quotient_sum"] == p + s - u
        assert row["if_no_quotient_one_possible"] == (u < s)
        assert row["quotient_one_count_lower_bound"] == max(0, u + 1 - s)
    old_counts = [4, 8, 7, 7, 6, 6, 5] if p % 4 == 1 else [7, 7, 7, 7, 5, 5]
    branches = rec["carried_local_branches"].values()
    assert [row["forced_count"] for row in branches] == [c - s for c in old_counts]
    for row in branches:
        assert row["forced_count"] > 0
        assert row["hard_edges"] + row["opposite_edges"] == rec["H"]
        assert 2 * row["hard_edges"] == rec["H"] + row["hT"]
        assert row["common_low_parallel_candidates"] == [row["P"]]
        assert row["low_parallel_upper_bound"] <= 9
    fresh = rec["fresh_quotient_two_branches"]
    assert [row["u"] for row in fresh] == list(range(s))
    for row in fresh:
        u = row["u"]
        assert row["excess"] == s - u - 1
        assert row["forced_mass"] == p + 7 - 2 * u
        assert row["forced_count"] == 4 - s - u
        assert row["low_row_count_at_least"] == m - row["excess"]
        assert row["new_equality_catalog_used"] is False
        assert row["parallel_cases"][-1]["P"] == 9
        assert row["parallel_cases"][-1]["forced_Q"] == 0
    last = rec["q_nozero_branch"]
    assert last["excess"] == s
    assert last["low_row_count_at_least"] == m - s
    assert last["forced_mass"] == p + 9
    assert last["forced_count"] == 5 - s
    assert rec["q_zero_case_retains_arbitrary_quotient_heights"]


def test_minimal_bridge_bounds_keep_unproved_all_size_quantifier(p):
    rec = proof.minimal_four_gap_consequences(p)
    assert rec["odd_minimal_H_lower_bound"] == 5 * p + 6
    assert rec["even_minimal_H_lower_bound"] is None
    assert rec["even_unconditional_bound_status"] == "RETRACTED_SCOPE_MISMATCH"
    assert rec["unconditional_even_H_frame_lower_bound"] == 2 * p + 2
    assert rec["even_without_level_two_H_lower_bound"] == (6 * p + 6 if p >= 37 else 6 * p - 10)
    assert rec["odd_without_level_three_H_lower_bound"] == (7 * p + 8 if p >= 37 else 7 * p - 10)
    assert rec["restricted_Type_I_15_750_closed"] is True
    assert rec["restricted_Type_I_15_750_required_G_size"] == 3 * p - 2
    assert rec["restricted_Type_I_15_750_required_Max_plus_identity"] == "S_G=3-2*f_e on Max+"
    assert rec["even_level_two_branch_uses_proved_Type_I_15_750"] is False
    assert rec["even_level_two_branch_closed"] is False
    assert "general odd-k" in rec["missing_even_global_bridge_implication"]
    assert rec["all_size_localization_proved"] is False
    assert rec["eventual_E1_proved"] is False


@pytest.mark.parametrize("p", [True, 23, 25, 33, 49, 29.0])
def test_invalid_primes_rejected(p):
    with pytest.raises(ValueError, match="prime p>=29"):
        proof.capacity_exclusion(p, 3)


@pytest.mark.parametrize("r", [True, 2, 6, 3.0])
def test_unsupported_shell_floor_rejected(r):
    with pytest.raises(ValueError, match="shell floor"):
        proof.capacity_exclusion(37, r)


@pytest.mark.parametrize("s", [True, 0, 3, 1.0])
def test_no_accidental_claim_of_later_residual_layers(s):
    with pytest.raises(ValueError, match="two proved layers"):
        proof.residual_two_layer_exclusion(29, s)


def test_absent_spectrum_cannot_prove_capacity(monkeypatch):
    monkeypatch.setattr(proof, "affine_parity_small_mass_spectrum", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="spectrum dependency"):
        proof.capacity_exclusion(37, 4)


def test_extra_low_mass_blocks_capacity(monkeypatch):
    actual = proof.affine_parity_small_mass_spectrum(37)
    monkeypatch.setattr(proof, "affine_parity_small_mass_spectrum", lambda p: {
        **actual, "union_allowed_masses": [0, 2, 34, 36, 38]})
    with pytest.raises(ArithmeticError, match="spectrum dependency"):
        proof.capacity_exclusion(37, 3)


def test_weakened_quotient_floor_leaves_a_collision(monkeypatch):
    actual = proof.quotient_floor
    monkeypatch.setattr(proof, "quotient_floor", lambda p, u: 1 if u == 2 else actual(p, u))
    with pytest.raises(ArithmeticError, match="residue separation"):
        proof.capacity_exclusion(37, 3)


def test_missing_prior_local_classification_blocks_carry(monkeypatch):
    monkeypatch.setattr(proof, "joint_layer_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="carried layer"):
        proof.residual_two_layer_exclusion(29, 1)


def test_false_carried_branch_blocks_closure(monkeypatch):
    actual = proof._carry_branch
    monkeypatch.setattr(proof, "_carry_branch", lambda *args: {**actual(*args), "proved": False})
    with pytest.raises(ArithmeticError, match="two-layer branch"):
        proof.residual_two_layer_exclusion(29, 2)


def test_false_fresh_branch_blocks_closure(monkeypatch):
    monkeypatch.setattr(proof, "_uncatalogued_branch", lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="two-layer branch"):
        proof.residual_two_layer_exclusion(31, 2)


def test_missing_official_entry_blocks_minimal_bridge(monkeypatch):
    monkeypatch.setattr(proof, "official_unit_entry_ledger", lambda *args: {"official_entry_proved": False})
    with pytest.raises(ArithmeticError, match="official bridge"):
        proof.minimal_four_gap_consequences(37)


@pytest.mark.parametrize("dependency", ["capacity_exclusion", "residual_two_layer_exclusion"])
def test_false_minimal_consequence_dependency_is_not_discarded(monkeypatch, dependency):
    monkeypatch.setattr(proof, dependency, lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="minimal consequence dependency"):
        proof.minimal_four_gap_consequences(37)


def test_missing_restricted_type_I_theorem_blocks_its_scope_receipt(monkeypatch):
    monkeypatch.setattr(proof, "type_I_multilevel_bad_case_closed_all_primes", lambda: False)
    with pytest.raises(ArithmeticError, match="minimal consequence dependency"):
        proof.minimal_four_gap_consequences(37)


def test_true_restricted_type_I_boolean_cannot_produce_general_even_bound(monkeypatch, p):
    monkeypatch.setattr(proof, "type_I_multilevel_bad_case_closed_all_primes", lambda: True)
    rec = proof.minimal_four_gap_consequences(p)
    assert rec["restricted_Type_I_15_750_closed"] is True
    assert rec["restricted_Type_I_15_750_required_G_size"] == 3 * p - 2
    assert rec["even_minimal_H_lower_bound"] is None
    assert rec["unconditional_even_H_frame_lower_bound"] == 2 * p + 2
    assert rec["even_level_two_branch_closed"] is False
    assert rec["even_without_level_two_H_lower_bound"] > 0
    assert rec["proved"] is True


@pytest.mark.parametrize("field,value", [
    ("shell_level_entry_proved", False),
    ("official_entry_proved", True),
    ("sharp_H_size_floor", 1),
    ("sharp_G_size_floor", 1),
    ("restricted_Type_I_15_750_required_G_size", 1),
    ("restricted_Type_I_15_750_required_Max_plus_identity", "s_plus=1 only"),
    ("restricted_Type_I_15_750_required_Max_minus_inequalities", "S_G<=-1 only"),
    ("restricted_Type_I_15_750_size_forced", True),
    ("restricted_Type_I_15_750_affine_identity_forced", True),
])
def test_even_range_and_affine_scope_guards_fail_closed(monkeypatch, field, value):
    actual = proof.official_unit_entry_ledger
    def changed(p, h_is_odd):
        rec = actual(p, h_is_odd)
        return rec if h_is_odd else {**rec, field: value}
    monkeypatch.setattr(proof, "official_unit_entry_ledger", changed)
    with pytest.raises(ArithmeticError, match="official bridge dependency"):
        proof.minimal_four_gap_consequences(37)


def test_false_minimal_consequence_blocks_proposition(monkeypatch):
    monkeypatch.setattr(proof, "minimal_four_gap_consequences", lambda *args: {"proved": False})
    assert proof.proposition_15774()["proved"] is False


def test_saved_evidence_and_all_global_flags():
    result = proof.proposition_15774()
    assert result == json.loads((ROOT / "evidence/e1_gmin_m4_prop15774.json").read_text())
    assert result["proved"]
    assert result["records_are_identity_replays_not_a_prime_census"]
    assert result["new_generic_frontier"] == "p>=29,t>=q+3,k>=5p+5, q=(p-1)/2"
    for name in ("residual_ii_closed_general", "minimal_four_gap_bridge_closed_general",
                 "eventual_E1_proved", "e1_closed_general", "original_MO_limit_closed"):
        assert result[name] is False
    assert (ROOT / result["proof_note"]).is_file()
    assert (ROOT / result["local_note"]).is_file()
