import json
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15762 import (
    exact_projection_audit,
    first_possible_shell_certificate,
    integral_eigenvector_norm_floor,
    low_gap_exclusion_certificate,
    projected_norm_squared,
    spherical_ceiling,
    theorem_record,
)
from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power


ROOT = Path(__file__).resolve().parents[1]


def test_projection_normalization_and_low_gap_exclusion_are_exact():
    assert spherical_ceiling(11) == 671
    assert projected_norm_squared(11, 669) == 22
    assert projected_norm_squared(11, 667) == 44
    assert projected_norm_squared(11, 663) == 88
    assert integral_eigenvector_norm_floor(11) == 12
    for p in (5, 7, 11):
        row = low_gap_exclusion_certificate(p)
        assert row["proved"]
        assert row["excluded_positive_gaps"] == [2, 4, 6]
        assert all(row["checks"].values())


def test_first_possible_gap_has_the_forced_sparse_signed_support():
    for p in (5, 11):
        row = first_possible_shell_certificate(p)
        assert row["gap"] == 8
        assert row["projection_z_norm_squared"] == 8 * p
        assert row["v_norm_squared"] == 2 * p
        assert row["v_support"] == 2 * p
        assert row["v_positive_count"] == p - 2
        assert row["v_negative_count"] == p + 2
        assert not row["existence_claimed"]


def test_exact_projection_identities_on_a_boolean_eigenvector_and_perturbation():
    p = 3
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    eigen = np.rint(halfspace_boolean_vector(p)).astype(np.int64)
    for x in (eigen, eigen * np.array([-1] + [1] * (len(eigen) - 1))):
        row = exact_projection_audit(C.tolist(), x.tolist(), p)
        assert row["z_is_minus_p_eigenvector"]
        assert row["sum_w_identity"]
        assert row["norm_identity"]
        assert row["one_coordinate_parity"]


def test_checked_in_evidence_and_canonical_docs_preserve_nonclosure():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15762.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["proved"]["non_eigen_boolean_gap_at_least_8_for_p_at_least_5"]
    assert not observed["proved"]["nonregularizable_conference_constructed"]
    assert not observed["proved"]["residual_ii_closed"]
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
        flat = " ".join(text.split())
        assert "15.762" in text, name
        assert "conference cube gap" in flat.lower(), name
        assert "residual (ii)" in flat and "OPEN" in flat, name
