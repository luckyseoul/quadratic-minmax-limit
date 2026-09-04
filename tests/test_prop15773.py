"""Joint-layer identities, complete quotient splits, and failure injection."""
import json
from pathlib import Path

import pytest

import e1_gmin_m4_prop15773 as proof

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module", params=[29, 31, 37, 43])
def record(request):
    return proof.joint_layer_exclusion(request.param)


def test_joint_all_prime_scope_and_chart(record):
    p = record["p"]
    assert record["original_k"] == 5 * p - 1
    assert record["layer_index_t"] == (p - 1) // 2
    assert record["H_edge_count"] == 5 * p
    assert record["all_boundary_sizes_excluded"]
    assert record["residual_ii_layer_excluded"]
    assert record["finite_prime_graph_or_slice_census_used"] is False
    chart = record["common_row_identity"]
    assert chart["I"] == 0
    assert chart["guaranteed_isolated_vertices"] == p * p - 10 * p + 1 > 0
    assert chart["E_z_i_z_j"] == f"-1/{p}"
    assert chart["valid_for_both_h_signs"]
    assert chart["equality_classification_or_offset_assumed"] is False
    assert chart["hard_sign_relative_to_transported_c_H"] == (1 if p % 4 == 1 else -1)


def test_residue_partition_recomputed_at_new_edge_count(record):
    p, m = record["p"], record["m"]
    q = m - 1
    ledger = record["residue_ledger"]
    assert [r["u"] for r in ledger["rows"]] == list(range(m))
    for r in ledger["rows"][:-1]:
        assert r["quotient_sum"] == p - r["u"]
        assert r["forced_low_count_at_least"] == r["u"] + 1
        assert r["forced_low_mean"] == 2 * r["u"] + p + 1
    assert ledger["surviving_positive_quotient_residues"] == (
        [0, q - 3, q - 2, q - 1] if p % 4 == 1 else [q - 2])
    zero = ledger["rows"][-1]
    assert zero["quotient_sum"] == m
    assert zero["exhaustive_alternatives"] == ["some_quotient_zero", "all_quotients_one"]
    assert zero["zero_case_low_mean"] == p - 1
    assert zero["flat_case_low_mean"] == 2 * p
    assert ledger["quotient_zero_is_not_forced_at_u_equals_q"]


def test_full_carry_forcing_counts_and_boundaries(record):
    p = record["p"]
    branches = record["carried_branch_exclusions"]
    counts = [4, 8, 7, 7, 6, 6, 5] if p % 4 == 1 else [7, 7, 7, 7, 5, 5]
    assert [r["forced_next_row_count_at_least"] for r in branches.values()] == counts
    for r in branches.values():
        assert all(r["checks"].values())
        assert 2 * r["hard_edge_count"] == 5 * p + r["hard_sign_times_global_T"]
        assert r["hard_edge_count"] + r["opposite_edge_count"] == 5 * p
        assert r["common_low_parallel_candidates"] == [r["coefficient_offset"]]
        boundary = r["opposite_local_exclusion"]
        assert [v[0] for v in boundary["nonzero_boundary_floor_excess_rows"]] == [2, p - 1]
        assert all(0 < v[2] < p - 3 for v in boundary["nonzero_boundary_floor_excess_rows"])
        assert boundary["last_pointwise_baseline"] == ("1-x_j" if p % 4 == 1 else "x_j")


def test_flat_all_ten_actual_parallel_counts_and_vacuous_negative_index(record):
    p, m = record["p"], record["m"]
    flat = record["flat_branch_exclusion"]
    assert flat["new_equality_classification_used"] is False
    assert flat["coefficient_offset_assumed"] is False
    assert flat["minimum_R"] == 9
    assert [r["P"] for r in flat["parallel_cases"]] == list(range(10))
    for r in flat["parallel_cases"]:
        P = r["P"]
        assert r["hard_sign_times_global_T"] == (p + 1) * P - 5 * p
        assert r["forced_next_Q"] == 9 - P
        assert r["forced_next_scaled_mean"] == p + 9
        assert r["surplus"] == m - 5
        assert r["forced_next_row_count_at_least"] == 5
        assert r["forbidden_Q_is_in_domain"] == (P < 9)
    last = flat["parallel_cases"][-1]
    assert last["formal_forbidden_Q"] == -1
    assert last["forced_next_Q"] == 0
    assert flat["P9_uses_Q_nonnegative_not_a_negative_index_row"]


@pytest.mark.parametrize("P", [2, 3, 4, 5])
def test_p3_two_excess_units_allow_both_elevated_profiles(P):
    p, m, u = 31, 16, 13
    single = [3] + [1] * (m - 1)
    double = [2, 2] + [1] * (m - 2)
    for quotients in (single, double):
        parallel = proof.normalized_parallel_profile(p, u, 1, P, quotients)
        assert sum(parallel) == m * P + 2
        assert parallel == [P + k - 1 for k in quotients]


@pytest.mark.parametrize("p,P", [(29, 4), (31, 3), (31, 4)])
def test_zero_case_keeps_arbitrarily_high_quotients(p, P):
    m = (p + 1) // 2
    quotients = [m] + [0] * (m - 1)
    parallel = proof.normalized_parallel_profile(p, m - 1, 0, P, quotients)
    assert parallel[0] == P + m > 9
    assert sum(parallel) == m * (P + 1)


@pytest.mark.parametrize("p", [True, False, 0, 13, 23, 25, 33, 49, 29.0])
def test_rejects_parameters_outside_claimed_family(p):
    with pytest.raises(ValueError, match="prime p>=29"):
        proof.joint_layer_exclusion(p)


@pytest.mark.parametrize("p", [29, 31])
def test_missing_mass_p9_theorem_blocks_flat_branch(monkeypatch, p):
    monkeypatch.setattr(proof, "p_plus_nine_local_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="opposite all-boundary"):
        proof.flat_mean_2p_exclusion(p)


def test_missing_p3_minus_one_lemma_blocks_residue_partition(monkeypatch):
    monkeypatch.setattr(proof, "p_minus_one_local_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="mass p-1"):
        proof.joint_residue_ledger(31)


def test_missing_p1_gap_four_offsets_cannot_be_dropped(monkeypatch):
    actual = proof.p1_hard_family_catalog(29)
    monkeypatch.setattr(proof, "p1_hard_family_catalog", lambda p: {**actual, "coefficient_offsets": [6]})
    with pytest.raises(ArithmeticError, match="hard offsets"):
        proof.joint_layer_exclusion(29)


def test_changed_phase_zero_boundary_floor_blocks_local_conclusion(monkeypatch):
    actual = proof.residual_even_floor_table(29)
    monkeypatch.setattr(proof, "residual_even_floor_table", lambda p: {
        **actual, "phase_zero_floors": {**actual["phase_zero_floors"], 4: 38}})
    with pytest.raises(ArithmeticError, match="opposite all-boundary"):
        proof.flat_mean_2p_exclusion(29)


def test_changed_carried_low_mean_cannot_be_reused(monkeypatch):
    actual = proof.p3_next_residue_ledger(31)
    actual["rows"][13]["forced_low_mean"] += 2
    monkeypatch.setattr(proof, "p3_next_residue_ledger", lambda p: actual)
    with pytest.raises(ArithmeticError, match="low-cell mean"):
        proof.joint_residue_ledger(31)


@pytest.mark.parametrize("dependency", ["common_row_identity", "joint_residue_ledger"])
def test_false_joint_setup_dependency_blocks_public_closure(monkeypatch, dependency):
    monkeypatch.setattr(proof, dependency, lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="joint identity or residue"):
        proof.joint_layer_exclusion(29)


@pytest.mark.parametrize("dependency", ["common_row_identity", "opposite_mass_exclusion"])
def test_false_flat_wrapper_dependency_blocks_closure(monkeypatch, dependency):
    monkeypatch.setattr(proof, dependency, lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="flat branch dependency"):
        proof.flat_mean_2p_exclusion(29)


def test_false_carried_wrapper_dependency_blocks_closure(monkeypatch):
    monkeypatch.setattr(proof, "opposite_mass_exclusion", lambda *args: {"proved": False})
    with pytest.raises(ArithmeticError, match="carried opposite exclusion dependency"):
        proof.joint_layer_exclusion(29)


def test_false_flat_branch_blocks_joint_aggregator(monkeypatch):
    monkeypatch.setattr(proof, "flat_mean_2p_exclusion", lambda p: {"proved": False})
    with pytest.raises(ArithmeticError, match="joint branch dependency"):
        proof.joint_layer_exclusion(29)


def test_false_carried_branch_blocks_joint_aggregator(monkeypatch):
    actual = proof._branch_ledger
    monkeypatch.setattr(proof, "_branch_ledger", lambda *args: {**actual(*args), "proved": False})
    with pytest.raises(ArithmeticError, match="joint branch dependency"):
        proof.joint_layer_exclusion(29)


def test_floor_proof_flag_is_not_discarded(monkeypatch):
    actual = proof.residual_even_floor_table(29)
    monkeypatch.setattr(proof, "residual_even_floor_table", lambda p: {**actual, "proved": False})
    with pytest.raises(ArithmeticError, match="phase-zero floor dependency"):
        proof.flat_mean_2p_exclusion(29)


def test_saved_live_evidence_and_open_global_gates():
    result = proof.proposition_15773()
    assert result == json.loads((ROOT / "evidence/e1_gmin_m4_prop15773.json").read_text())
    assert result["proved"] and result["status"] == "PROVED_INFINITE_FAMILY"
    assert result["records_are_identity_replays_not_exhaustive_prime_evidence"]
    assert result["new_mean_2p_equality_classification_used"] is False
    assert result["symbolic_range"] == {"minimum_prime": 29, "p_mod_4": [1, 3],
                                         "t": "(p-1)/2", "k": "5p-1"}
    assert result["new_generic_frontier"] == "p>=29,t>=(p+1)/2,k>=5p+1"
    assert (ROOT / result["proof_note"]).is_file()
    assert result["residual_ii_closed_general"] is False
    assert result["e1_closed_general"] is False
    assert result["original_MO_limit_closed"] is False
