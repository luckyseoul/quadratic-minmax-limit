import json
from pathlib import Path

from e1_gmin_m4_prop15765 import (
    P,
    boolean_shadow,
    conference_action,
    conference_quadratic,
    exceptional_complement,
    exceptional_set,
    indicator_convolution_identity,
    integral_eigenvector,
    kiss_somlai_base,
    line_profile,
    paley_character,
    paley_neighbor_counts,
    points,
    special_direction_profiles,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_cited_matrix_and_disjoint_line_have_the_claimed_sizes():
    base = kiss_somlai_base()
    assert len(base) == 33
    assert all((x, 1) not in base for x in range(P))
    assert len(exceptional_complement()) == 44
    assert len(exceptional_set()) == 77


def test_exactly_four_special_directions_and_all_are_paley_square():
    e_set = exceptional_complement()
    profiles = special_direction_profiles(e_set)
    assert set(profiles) == {"0", "2", "9", "infinity"}
    for direction in (1, 3, 4, 5, 6, 7, 8, 10):
        assert line_profile(e_set, direction) == (4,) * P
    assert paley_character((1, 0)) == 1
    assert paley_character((1, 2)) == 1
    assert paley_character((1, 9)) == 1
    assert paley_character((0, 1)) == 1


def test_nonaffineness_has_a_direct_line_profile_certificate():
    e_set = exceptional_complement()
    union_of_four_lines_profile = tuple(sorted((11,) * 4 + (0,) * 7))
    observed = [line_profile(e_set, d) for d in (*range(P), None)]
    assert all(tuple(sorted(profile)) != union_of_four_lines_profile for profile in observed)
    # If D were seven parallel lines, its complement E would be four.
    assert len(special_direction_profiles(e_set)) == 4 != 1


def test_convolution_and_equitable_partition_are_pointwise_exact():
    assert indicator_convolution_identity()
    d_set = exceptional_set()
    assert paley_neighbor_counts(d_set) == {"inside": (40,), "outside": (35,)}
    assert len(points()) == 121


def test_unique_three_integral_eigenvector_and_first_boolean_defect():
    y = integral_eigenvector()
    x = boolean_shadow()
    assert len(y) == len(x) == 122
    assert [i for i, value in enumerate(y) if abs(value) == 3] == [0]
    assert conference_action(y) == tuple(11 * value for value in y)
    assert sum(value * value for value in y) == 130
    phi = 11 * 122 // 2
    assert conference_quadratic(x) == 649
    assert phi - conference_quadratic(x) == 22 == 2 * P


def test_checked_in_evidence_and_scope_guards():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15765.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["proved"]["nonaffine_first_defect_shell_point_exists"]
    assert not observed["proved"]["all_first_shell_points_are_affine"]
    assert not observed["proved"]["common_all_deletions_H_constructed"]
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
        flat = " ".join(text.split()).lower()
        assert "15.765" in text, name
        assert "nonaffine" in flat, name
        assert "residual (ii)" in flat and "open" in flat, name
