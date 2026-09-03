from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from scripts.original_mo_two_half_geometry import (
    INSTANCES,
    analyze_orientation,
    boolean_states,
    directed_cut_halves,
    directed_halfcut_norms,
    quadratic_values,
    rademacher_abs_mean,
    sharp_influence_constant,
    skew_from_bits,
    solve_exhaustive,
)


def _analyze(n: int) -> dict[str, object]:
    item = INSTANCES[n]
    return analyze_orientation(
        np.asarray(item["A"], dtype=np.int64),
        int(item["m"]),
        skew_from_bits(n, str(item["r_bits"])),
    )


def test_pinned_two_half_envelopes() -> None:
    expected = {
        5: ({0: 8, 4: 12, 8: 12}, 16, 8),
        6: ({0: 8, 2: 10, 4: 12, 8: 16, 10: 14}, 18, 8),
        7: ({0: 4, 4: 8, 8: 12, 12: 16, 16: 20}, 22, 4),
        8: ({0: 4, 2: 10, 4: 12, 6: 14, 8: 16, 10: 18, 12: 16, 14: 18, 16: 16, 18: 18, 20: 16}, 28, 8),
    }
    for n, (envelope, B, excess) in expected.items():
        result = _analyze(n)
        assert result["phi_A"] == INSTANCES[n]["m"]
        assert result["B"] == B
        assert result["max_W_minus_D"] == excess
        A = np.asarray(INSTANCES[n]["A"], dtype=np.int64)
        R = skew_from_bits(n, str(INSTANCES[n]["r_bits"]))
        assert 2 * int(np.max(directed_halfcut_norms(A, R))) == B
        assert result["directed_halfcut"]["identity_twice_max_phi_equals_B"] is True
        assert {int(d): w for d, w in result["envelope_Wmax_by_D"].items()} == envelope


def test_zero_error_target_changes_at_seven() -> None:
    for n in (5, 6):
        assert not _analyze(n)["zero_error_doubling_passes"]
    for n in (7, 8):
        result = _analyze(n)
        assert result["zero_error_doubling_passes"]
        assert result["max_W_minus_D"] <= 2 * (math.sqrt(2) - 1) * result["M"]


def test_endpoint_slack_and_directed_cancellation_identities() -> None:
    """Exhaust every Boolean pair for a mixed four-vertex example."""
    A = np.array(
        [
            [0, 1, -1, 1],
            [1, 0, 1, -1],
            [-1, 1, 0, 1],
            [1, -1, 1, 0],
        ],
        dtype=np.int64,
    )
    R = np.array(
        [
            [0, 1, -1, 1],
            [-1, 0, 1, -1],
            [1, -1, 0, 1],
            [-1, 1, -1, 0],
        ],
        dtype=np.int64,
    )
    states = boolean_states(4)
    values = quadratic_values(A, states)
    M = int(np.max(np.abs(values)))
    for xi, x_state in enumerate(states):
        for yi, y_state in enumerate(states):
            qx = int(values[xi])
            qy = int(values[yi])
            epsilon = M - max(abs(qx), abs(qy))
            slack = 2 * M - abs(qx + qy)
            separation = abs(qx - qy)
            assert slack == separation + 2 * epsilon

            outward, inward = directed_cut_halves(A, R, x_state, y_state)
            cross = int(x_state @ R @ y_state)
            assert cross == 2 * (inward - outward)
            signed_overlap = 4 * min(abs(outward), abs(inward))
            if outward * inward >= 0:
                signed_overlap *= -1
            assert abs(cross) - separation == signed_overlap


def test_n5_complete_orientation_enumeration() -> None:
    item = INSTANCES[5]
    result = solve_exhaustive(np.asarray(item["A"], dtype=np.int64), int(item["m"]))
    assert result["status"] == "EXHAUSTIVE"
    assert result["orientations_checked"] == 2**10
    assert result["objective_B"] == item["expected_B"]
    assert result["optimal_orientation_count"] == item["expected_optimal_count"]


def test_persisted_solver_records_exclude_volatile_telemetry() -> None:
    """Certificate hashes must not depend on CP-SAT scheduling or runtime."""
    volatile = {"wall_time_seconds", "branches", "conflicts"}
    source = (Path(__file__).resolve().parents[1] / "scripts" / "original_mo_two_half_geometry.py").read_text(
        encoding="utf-8"
    )
    result_block = source[source.index('result = {', source.index('def solve_cpsat')):]
    result_block = result_block[:result_block.index('if status != cp_model.OPTIMAL')]
    assert all(f'"{field}"' not in result_block for field in volatile)


def test_sharp_influence_reparametrization() -> None:
    expected_mu = {1: 1, 2: 1, 3: 1.5, 4: 1.5, 5: 1.875}
    for k, value in expected_mu.items():
        assert rademacher_abs_mean(k) == value
    assert sharp_influence_constant(5, 4) == 1.875
    assert sharp_influence_constant(6, 5) == 2.25
    assert sharp_influence_constant(8, 10) == 1.75
    assert sharp_influence_constant(11, 17) == 108.28125 / 68


def test_n6_conference_is_counterexample_to_literal_degree_bound() -> None:
    A = np.asarray(INSTANCES[6]["A"], dtype=np.int64)
    assert np.array_equal(A @ A, 5 * np.eye(6, dtype=np.int64))
    values = quadratic_values(A, boolean_states(6))
    assert int(np.max(np.abs(values))) == 5
    # f=Q_A/5 has degree two, sup norm one, and total L1 influence 9/4.
    assert sharp_influence_constant(6, 5) == 9 / 4

    for n in range(2, 9):
        states = boolean_states(n)
        A = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
        A[np.triu_indices(n, 2)] *= -1
        A = np.triu(A, 1) + np.triu(A, 1).T
        q = quadratic_values(A, states)
        total_influence = 0.0
        for i in range(n):
            flipped = states.copy()
            flipped[:, i] *= -1
            derivative = (q - quadratic_values(A, flipped)) // 2
            total_influence += float(np.mean(np.abs(derivative)))
        assert total_influence == n * rademacher_abs_mean(n - 1)


def test_fixed_tie_best_response_first_moment_and_conference_output() -> None:
    # Proposition 1 requires a tie value fixed independently of the input.
    for n in (5, 6, 7):
        A = np.asarray(INSTANCES[n]["A"], dtype=np.int64)
        states = boolean_states(n)
        fields = states @ A
        response = np.where(fields >= 0, 1, -1)
        expected = rademacher_abs_mean(n - 1) / (n - 1)
        for i in range(n):
            for j in range(n):
                if i != j:
                    assert np.mean(response[:, i] * states[:, j]) == A[i, j] * expected
        action = np.einsum(
            "bi,ij,bj->b", response, A, states, dtype=np.int64
        )
        assert float(np.mean(action)) == n * rademacher_abs_mean(n - 1)

    # For the order-six symmetric conference signing, synchronous outputs
    # are pairwise unbiased and their mean quadratic energy is zero.
    A6 = np.asarray(INSTANCES[6]["A"], dtype=np.int64)
    states6 = boolean_states(6)
    response6 = np.where(states6 @ A6 >= 0, 1, -1)
    correlations = response6.T @ response6
    assert np.array_equal(correlations, len(states6) * np.eye(6, dtype=np.int64))
    assert float(np.mean(quadratic_values(A6, response6))) == 0.0
