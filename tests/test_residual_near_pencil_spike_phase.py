import pytest

from residual_near_pencil_spike_phase import (
    endpoint_phase_ledger,
    representative_alignment_audit,
    representative_phase_audit,
)


@pytest.mark.parametrize("p", [53, 61, 73, 89])
def test_lower_endpoint_forces_positive_boundary_product(p):
    row = endpoint_phase_ledger(p)
    assert row["proved"]
    assert row["necessary_boundary_product"] == 1
    assert row["target_feature_product"] == row["Paley_edge_product"]


@pytest.mark.parametrize(
    ("kind", "boundary_size", "translations"),
    [
        ("ordinary", 1504, [[0, 0], [0, 1]]),
        ("triple", 1500, [[0, 0], [3, 0]]),
    ],
)
def test_p53_same_spike_incidence_has_both_boundary_phases(
    kind, boundary_size, translations
):
    row = representative_alignment_audit(kind)
    assert row["proved_for_this_prime"]
    assert row["boundary_size"] == boundary_size
    assert [item["translation"] for item in row["rows"]] == translations
    assert [item["boundary_product"] for item in row["rows"]] == [1, -1]
    assert all(item["factorizations_exact"] for item in row["rows"])
    assert all(
        item["finite_spike_incidence_with_D"] == [True, True]
        for item in row["rows"]
    )
    assert row["translation_preserves_circle_mismatch"]
    assert row["circle_mismatch_does_not_determine_boundary_phase"]
    assert row["spike_incidence_does_not_determine_boundary_phase"]


def test_phase_audit_does_not_claim_residual_closure():
    row = representative_phase_audit()
    assert row["all_checks"]
    assert row["phase_compatible_alignments_exist"]
    assert row["phase_incompatible_alignments_exist"]
    assert not row["residual_ii_closed"]


@pytest.mark.parametrize("p", [3, 7, 13, 17, 29, 49])
def test_phase_ledger_rejects_wrong_prime_class(p):
    with pytest.raises(ValueError, match="prime p=1 mod 4 with p>=53"):
        endpoint_phase_ledger(p)
