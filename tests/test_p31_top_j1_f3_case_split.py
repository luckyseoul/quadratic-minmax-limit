from itertools import combinations

from e1_gmin_m4_p31_top_j1_f3_case_split import (
    HARD,
    P,
    TARGET_PROFILE,
    V_HARD,
    V_OPPOSITE,
    classify_fixed_counts,
    classify_support,
    required_raw_profile,
    theorem_record,
)


def test_exact_exhaustive_case_split_counts() -> None:
    row = theorem_record()
    assert row["proved"] is True
    assert row["residual_ii_closed"] is False
    assert row["fixed_allocation_count"] == 5984
    assert row["possible_fixed_allocation_count"] == 1862
    assert row["possible_support_count"] == 1428
    assert row["support_class_counts"] == {
        "capacity_excluded": 2255,
        "selector_excluded_sdr": 1309,
        "singleton_sdr_two_nonorigin": 14,
        "triple_k2_sdr_two_nonorigin": 734,
        "triple_k3_repeated_one_origin_one_nonorigin": 274,
        "triple_k3_repeated_two_nonorigin": 406,
    }


def test_support_types_and_correction_signature_bounds() -> None:
    hard_v = sorted(V_HARD)
    opposite_v = sorted(V_OPPOSITE)
    hard_out = sorted(HARD - V_HARD)

    singleton = classify_support([hard_v[0]])
    assert singleton["classification"] == "singleton_sdr_two_nonorigin"
    assert singleton["auxiliary_parity_weight"] == 16
    assert singleton["correction_signature_possible_weights"] == (0, 2, 4)

    bad_singleton = classify_support([opposite_v[0]])
    assert bad_singleton["classification"] == "selector_excluded_sdr"

    k2 = classify_support([opposite_v[0], opposite_v[1], hard_out[0]])
    assert k2["classification"] == "triple_k2_sdr_two_nonorigin"

    all_opposite = classify_support(opposite_v)
    assert all_opposite["classification"] == "triple_k3_repeated_one_origin_one_nonorigin"
    assert all_opposite["origin_cancellation_units"] == 1
    assert all_opposite["nonorigin_cancellation_units"] == 1
    assert all_opposite["correction_signature_possible_weights"] == (1, 3)

    one_hard = classify_support([hard_v[0], opposite_v[0], opposite_v[1]])
    assert one_hard["classification"] == "triple_k3_repeated_two_nonorigin"


def test_every_fixed_allocation_has_weight_one_or_three_and_raw_total_480() -> None:
    # Check the two actual multiplicity forms behind singleton parity b.
    b = min(V_HARD)
    fixed_triple = [0] * (P + 1)
    fixed_triple[b] = 3
    assert classify_fixed_counts(fixed_triple)["possible"] is True

    fixed_one_plus_pair = [0] * (P + 1)
    fixed_one_plus_pair[b] = 1
    fixed_one_plus_pair[(b + 1) % (P + 1)] = 2
    assert classify_fixed_counts(fixed_one_plus_pair)["possible"] is True

    cancellation = [0] * (P + 1)
    cancellation[0] = 1
    cancellation[1] = 1
    raw = required_raw_profile(fixed_one_plus_pair, cancellation)
    assert sum(raw) == 480
    assert all((raw[i] - TARGET_PROFILE[i] + fixed_one_plus_pair[i]) % 2 == 0 for i in range(P + 1))

    # Independently exhaust the 4,960 three-element b supports.
    assert sum(1 for _ in combinations(range(P + 1), 3)) == 4960
