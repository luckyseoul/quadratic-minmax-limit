import json
from fractions import Fraction
from math import comb
from pathlib import Path

import pytest

from e1_gmin_m4_prop15763 import (
    affine_alias_bound_ledger,
    affine_alias_parameters,
    critical_alias_avoidance_bound,
    least_odd_integer_at_least,
    prop15755_unsigned_bound,
    signed_affine_alias_bound,
    signed_incidence_feasibility,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_general_odd_r_alias_incidence_is_exact():
    for p, r in ((5, 1), (11, 1), (11, 3), (11, 5), (17, 7)):
        row = affine_alias_parameters(p, r)
        m = (p + 1) // 2 + r
        assert row["m"] == m
        assert row["aliases"] == comb(m, r)
        assert row["max_edge_multiplicity"] == 2 * comb(m - 2, r - 1)
        assert row["positive_to_outside_multiplicity"] <= row[
            "max_edge_multiplicity"
        ]


def test_r_one_bound_and_odd_parity_are_closed_form():
    for p in (5, 7, 11, 13, 17):
        old = prop15755_unsigned_bound(p, 1)
        signed = signed_affine_alias_bound(p, 1)
        row = affine_alias_bound_ledger(p, 1)
        assert old == Fraction((p + 1) * (p + 3), 8)
        assert signed == Fraction(p * p + 11, 4)
        assert signed.denominator == 1
        assert signed.numerator % 2 == 1
        assert row["signed_odd_integer_bound"] == signed
        assert row["signed_strictly_improves_parity_adjusted_15755"]


def test_general_signed_comparison_and_effective_bound_are_exact():
    for p, r in ((11, 3), (11, 5), (13, 5), (19, 7)):
        old = prop15755_unsigned_bound(p, r)
        signed = signed_affine_alias_bound(p, r)
        row = affine_alias_bound_ledger(p, r)
        assert signed == 2 * old - p * r * r + 2
        assert row["effective_odd_integer_bound"] == max(
            least_odd_integer_at_least(old),
            least_odd_integer_at_least(signed),
            p * r * r - 2,
        )


def test_signed_count_inequality_changes_truth_at_the_exact_odd_threshold():
    for p, r in ((5, 1), (11, 1), (11, 3)):
        threshold = least_odd_integer_at_least(signed_affine_alias_bound(p, r))
        below = signed_incidence_feasibility(p, r, threshold - 2)
        at = signed_incidence_feasibility(p, r, threshold)
        assert not below["incidence_compatible"]
        assert at["incidence_compatible"]
        assert below["proved_equivalent"] and at["proved_equivalent"]


def test_no_critical_internal_alias_forces_the_plus_three_bound():
    for p, r in ((5, 1), (7, 1), (11, 1), (11, 3)):
        base = signed_affine_alias_bound(p, r)
        avoidance = critical_alias_avoidance_bound(p, r)
        row = affine_alias_bound_ledger(p, r)
        m = row["m"]
        expected = (
            Fraction((p * r * r + 3) * m * (m - 1), 2 * r * (m - r))
            - p * r * r
            + 2
        )
        assert avoidance == expected > base
        assert row["critical_alias_alternative"]["avoidance_rational_bound"] == str(
            expected
        )
    assert critical_alias_avoidance_bound(7, 1) == 20
    assert least_odd_integer_at_least(critical_alias_avoidance_bound(7, 1)) == 21


def test_parameter_guards_reject_nonprime_or_nonadmissible_inputs():
    for p, r in ((4, 1), (9, 1), (5, 0), (5, 2), (5, 3)):
        with pytest.raises(ValueError):
            affine_alias_parameters(p, r)


def test_checked_in_evidence_and_canonical_docs_preserve_scope():
    observed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15763.json").read_text()
    )
    assert observed == theorem_record()
    assert observed["proved"]["signed_affine_alias_incidence_bound"]
    assert not observed["proved"]["all_minimum_shell_representatives_are_affine"]
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
        assert "15.763" in text, name
        assert "signed affine" in flat, name
        assert "residual (ii)" in flat and "open" in flat, name
