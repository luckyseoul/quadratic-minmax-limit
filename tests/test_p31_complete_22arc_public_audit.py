import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p31_complete_22arc_public_audit.py"
SPEC = importlib.util.spec_from_file_location("p31_complete_22arc_public_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

CERTIFICATE = AUDIT.build_certificate()


def test_all_eleven_public_representatives_are_exact_complete_arcs():
    assert CERTIFICATE["result_status"] == "exhaustive finite certificate"
    assert CERTIFICATE["corrected_projective_class_count"] == 12
    assert CERTIFICATE["public_representative_count"] == 11
    assert CERTIFICATE["missing_representative_count"] == 1
    assert CERTIFICATE["arc_size"] == 22
    assert CERTIFICATE["secants_per_arc"] == 231
    assert CERTIFICATE["outside_points_per_arc"] == 971
    assert CERTIFICATE["outside_secant_incidence_moment"] == 6930
    for row in CERTIFICATE["representatives"]:
        assert row["point_count"] == 22
        assert row["secant_line_count"] == 231
        assert row["outside_point_count"] == 971
        assert row["outside_secant_incidence_count"] == 6930
        assert row["minimum_outside_secant_index"] >= 1
        assert row["is_arc"] is True
        assert row["is_complete"] is True


def test_exact_c1_sequence_excludes_every_public_class_but_not_the_twelfth():
    assert CERTIFICATE["endpoint_required_c1"] == 10
    assert CERTIFICATE["index_one_point_counts_by_public_class"] == [
        0,
        2,
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        0,
        0,
    ]
    assert CERTIFICATE["maximum_c1_over_public_classes"] == 2
    assert CERTIFICATE["all_public_classes_excluded_by_c1_requirement"] is True
    assert CERTIFICATE["twelfth_class_representative_available"] is False
    assert CERTIFICATE["twelfth_class_audited"] is False
    assert CERTIFICATE["all_twelve_classes_excluded"] is False
    assert CERTIFICATE["p31_endpoint_closed"] is False


def test_positive_c1_points_share_one_secant_and_matching_number_is_one():
    assert CERTIFICATE["classes_with_positive_c1"] == [2, 8, 9]
    assert (
        CERTIFICATE[
            "positive_c1_classes_with_all_index_one_points_on_one_secant"
        ]
        == [2, 8, 9]
    )
    assert CERTIFICATE["maximum_disjoint_unique_secant_matching_by_public_class"] == [
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        0,
        0,
    ]
    rows = {row["public_class_index"]: row for row in CERTIFICATE["representatives"]}
    expected = {
        2: {
            "line": [1, 28, 8],
            "arc_pair": [[1, 13, 28], [1, 14, 9]],
            "index_one_outside_points": [[1, 24, 5], [1, 28, 22]],
        },
        8: {
            "line": [0, 1, 16],
            "arc_pair": [[1, 0, 0], [1, 21, 20]],
            "index_one_outside_points": [[1, 19, 24], [1, 20, 22]],
        },
        9: {
            "line": [1, 2, 20],
            "arc_pair": [[1, 6, 4], [1, 9, 13]],
            "index_one_outside_points": [[1, 22, 21], [1, 25, 30]],
        },
    }
    for class_index, unique_secant in expected.items():
        row = rows[class_index]
        assert row["index_one_point_count"] == 2
        assert row["unique_secant_line_count"] == 1
        assert row["index_one_points_per_unique_secant_line_histogram"] == {"2": 1}
        assert row["maximum_disjoint_unique_secant_matching"] == 1
        assert row["index_one_points_share_one_secant"] is True
        assert row["unique_secants"] == [unique_secant]


def test_checked_in_json_is_exactly_reproducible():
    checked_in = json.loads(
        (ROOT / "evidence" / "p31_complete_22arc_public_11_audit.json").read_text()
    )
    assert checked_in == CERTIFICATE
