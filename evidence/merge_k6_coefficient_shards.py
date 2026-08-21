#!/usr/bin/env python3
"""Merge exact shard reports from k6_coefficient_sieve_fast.py."""
from __future__ import annotations

import argparse
import glob
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


def merge_reports(paths: list[Path]) -> dict:
    reports = [json.loads(path.read_text()) for path in paths]
    if not reports:
        raise ValueError("no shard reports")
    first = reports[0]
    shard_count = first["shard_count"]
    if len(reports) != shard_count:
        raise ValueError(f"expected {shard_count} shards, found {len(reports)}")
    if {report["shard_index"] for report in reports} != set(range(shard_count)):
        raise ValueError("shard indices are not exactly 0..shard_count-1")

    invariant = (
        "p",
        "algorithm",
        "shard_count",
        "n_square_directions",
        "n_direction_subsets",
        "n_cyclic_subset_orbits",
        "normalized_total_T",
        "minimum_quartic_b",
        "energy_partitions",
        "relevant_type_histogram",
        "QVAR_threshold",
    )
    for report in reports[1:]:
        for key in invariant:
            if report[key] != first[key]:
                raise ValueError(f"inconsistent shard field {key}")

    processed_orbits = sum(
        report["processed_cyclic_subset_orbits"] for report in reports
    )
    processed_subsets = sum(
        report["processed_direction_subsets"] for report in reports
    )
    if processed_orbits != first["n_cyclic_subset_orbits"]:
        raise ValueError("merged orbit count is incomplete")
    if processed_subsets != first["n_direction_subsets"]:
        raise ValueError("merged direction-subset count is incomplete")

    candidate_histogram = Counter()
    quartic_histogram = Counter()
    for report in reports:
        candidate_histogram.update(
            {int(key): value for key, value in report["coefficient_candidate_histogram"].items()}
        )
        quartic_histogram.update(
            {int(key): value for key, value in report["abs_Zpsi_sq_histogram"].items()}
        )
    representatives = sum(
        report["boolean_representatives_mod_translation"] for report in reports
    )
    moment = (
        Fraction(
            sum(value * count for value, count in quartic_histogram.items()),
            representatives,
        )
        if representatives
        else None
    )
    threshold = Fraction(first["QVAR_threshold"])
    return {
        "p": first["p"],
        "algorithm": first["algorithm"],
        "shards_merged": shard_count,
        "n_square_directions": first["n_square_directions"],
        "n_direction_subsets": first["n_direction_subsets"],
        "n_cyclic_subset_orbits": first["n_cyclic_subset_orbits"],
        "normalized_total_T": first["normalized_total_T"],
        "minimum_quartic_b": first["minimum_quartic_b"],
        "energy_partitions": first["energy_partitions"],
        "relevant_type_histogram": first["relevant_type_histogram"],
        "total_type_tuples_before_coefficient_sieve": sum(
            report["total_type_tuples_before_coefficient_sieve"]
            for report in reports
        ),
        "coefficient_candidate_histogram": {
            str(key): value for key, value in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": sum(
            report["total_coefficient_candidates"] for report in reports
        ),
        "endpoint_branches": sum(report["endpoint_branches"] for report in reports),
        "boolean_representatives_mod_translation": representatives,
        "eps_plus_count_including_translations": representatives * first["p"] ** 2,
        "abs_Zpsi_sq_histogram": {
            str(key): value for key, value in sorted(quartic_histogram.items())
        },
        "E_abs_Zpsi_sq": str(moment) if moment is not None else None,
        "QVAR_threshold": str(threshold),
        "clears_QVAR": bool(moment >= threshold if moment is not None else True),
        "k6_empty": representatives == 0,
    }


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("pattern")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    paths = sorted(Path(path) for path in glob.glob(args.pattern))
    report = merge_reports(paths)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
