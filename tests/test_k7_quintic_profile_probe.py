"""Exact first-stage checks for the k=7 depressed-quintic sieve."""
from __future__ import annotations

import sys
from pathlib import Path

EVIDENCE = Path(__file__).parents[1] / "evidence"
sys.path.insert(0, str(EVIDENCE))

import k7_quintic_profile_probe as K7  # noqa: E402


def test_p13_quintic_profile_minima_and_kernel_ladder():
    report = K7.scan_prime(13)
    assert report["normal_form"] == "a*s^5+c*s^3+d*s^2+e*s+f"
    assert report["normalized_total_T"] == 21
    assert report["minimum_profile_energy"] == 1
    assert report["relevant_type_histogram"]["1"] == 12
    assert len(report["energy_partitions"]) == 105
    assert not report["empty_by_energy"]
    audit = report["rank_audit"]
    assert audit["n_direction_subsets"] == 1
    assert audit["kernel_dimension_histograms"] == {
        "1": {"5": 1},
        "2": {"4": 1},
        "3": {"3": 1},
        "4": {"2": 1},
        "5": {"1": 1},
    }
    assert audit["top_kernel_full_support"]
    assert audit["translation_removes_degree_four_kernel"]
