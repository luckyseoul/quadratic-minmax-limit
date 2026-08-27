#!/usr/bin/env python3
"""Independent coverage audit for all exceptional mod-7 tuple certificates."""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path


ORBIT_INDICES = (0, 5, 8, 25, 26, 30, 31)


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def is_supported_unknown(batch: dict, leaf: dict) -> bool:
    if leaf.get("solver_status") == "INFEASIBLE":
        return False
    return all(
        not (int(row["b"]) == 0 and int(row["phase"]) == 0 and int(mean) > 16)
        for row, mean in zip(batch["direction_rows"], leaf["scaled_means_direction_order"])
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.evidence_dir
    projection_files = {
        profile: root / f"p7_exceptional_{profile}22_all.json"
        for profile in ("tail", "mid", "head")
    }
    projections = {profile: json.loads(path.read_text()) for profile, path in projection_files.items()}
    selected_rows = []
    for profile, payload in projections.items():
        if (
            payload.get("experiment") != "p7_exceptional_tail22_catalogs"
            or payload.get("projection_profile") != profile
            or int(payload.get("orbit_count", -1)) != 7
            or int(payload.get("projection_group_order", -1)) != 7**22
        ):
            raise ValueError(f"invalid {profile} projection summary")
        rows = tuple(int(value) for value in payload["selected_dependency_rows_7"])
        if len(rows) != 22 or len(set(rows)) != 22:
            raise AssertionError(f"invalid {profile} dependency row set")
        selected_rows.append(set(rows))
    if any(selected_rows[i] & selected_rows[j] for i in range(3) for j in range(i)):
        raise AssertionError("projection dependency row sets overlap")

    projection_by_boundary = {
        profile: {
            tuple(row["fixed_boundary"]): row
            for row in payload["orbits"]
        }
        for profile, payload in projections.items()
    }
    orbit_rows = []
    total_initial_infeasible = 0
    total_initial_unknown = 0
    total_supported = 0
    total_high_mean = 0
    for orbit_index in ORBIT_INDICES:
        mean_path = root / "mean_batches" / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"
        batch = json.loads(mean_path.read_text())
        if (
            batch.get("experiment") != "p7_fixed_boundary_mean_allocation_batch"
            or int(batch.get("allocation_count", -1)) != 180
            or len(batch.get("leaves", [])) != 180
        ):
            raise ValueError(f"invalid mean batch for orbit {orbit_index}")
        boundary = tuple(int(value) for value in batch["fixed_boundary"])
        infeasible = {int(leaf["leaf_index"]) for leaf in batch["leaves"] if leaf.get("solver_status") == "INFEASIBLE"}
        unknown = {int(leaf["leaf_index"]) for leaf in batch["leaves"] if leaf.get("solver_status") == "UNKNOWN"}
        feasible = {int(leaf["leaf_index"]) for leaf in batch["leaves"] if leaf.get("solver_status") not in ("INFEASIBLE", "UNKNOWN")}
        if feasible or infeasible & unknown or len(infeasible | unknown) != 180:
            raise AssertionError(f"mean batch status partition failed for orbit {orbit_index}")
        supported = {
            int(leaf["leaf_index"])
            for leaf in batch["leaves"]
            if is_supported_unknown(batch, leaf)
        }
        high_mean = unknown - supported
        if supported & high_mean or supported | high_mean != unknown:
            raise AssertionError("unknown leaf partition failed")

        directory = root / f"p7_exceptional_mod7triple_orbit{orbit_index:02d}"
        if not directory.exists() and orbit_index == 0:
            directory = root / "p7_exceptional_mod7triple_orbit0"
        leaf_files = sorted(directory.glob("leaf*.json"))
        certificates = {int(path.stem.removeprefix("leaf")): json.loads(path.read_text()) for path in leaf_files}
        if set(certificates) != supported:
            raise AssertionError(
                f"GPU certificate coverage mismatch for orbit {orbit_index}: "
                f"missing={sorted(supported-set(certificates))}, extra={sorted(set(certificates)-supported)}"
            )
        expected_catalog_hashes = [
            projection_by_boundary[profile][boundary]["sha256"]
            for profile in ("tail", "mid", "head")
        ]
        for leaf_index, certificate in certificates.items():
            leaf = batch["leaves"][leaf_index]
            if (
                certificate.get("status") != "complete_exact_selected_dependency_gpu_join"
                or certificate.get("projection_mode") != "injective_disjoint_mod7_22x3_tuple"
                or int(certificate.get("exact_projected_matches", -1)) != 0
                or not certificate.get("projected_modularly_infeasible")
                or not certificate.get("finite_mean_allocation_exclusion")
                or tuple(certificate.get("fixed_boundary", [])) != boundary
                or certificate.get("fixed_scaled_means") != leaf["scaled_means_direction_order"]
                or certificate.get("projection_catalog_sha256") != expected_catalog_hashes
            ):
                raise AssertionError(f"invalid GPU certificate orbit {orbit_index} leaf {leaf_index}")

        summary_candidates = (
            [root / "p7_exceptional_mod7triple_orbit0_summary_complete.json"]
            if orbit_index == 0
            else [root / f"p7_exceptional_mod7triple_orbit{orbit_index:02d}_summary.json"]
        )
        summaries = [json.loads(path.read_text()) for path in summary_candidates]
        summary = next(
            (row for row in summaries if int(row.get("eligible_leaf_count", -1)) == len(supported)
             and int(row.get("deferred_supported_leaf_count", -1)) == 0),
            None,
        )
        if summary is None or int(summary.get("excluded_leaf_count", -1)) != len(supported) or int(summary.get("unresolved_leaf_count", -1)) != 0:
            raise AssertionError(f"missing complete batch summary for orbit {orbit_index}")
        summary_indices = {int(row["leaf_index"]) for row in summary["results"]}
        if summary_indices != supported or any(not row["excluded"] or int(row["matches"]) != 0 for row in summary["results"]):
            raise AssertionError(f"batch summary coverage failed for orbit {orbit_index}")

        orbit_rows.append({
            "orbit_index": orbit_index,
            "fixed_boundary": list(boundary),
            "initial_infeasible_leaves": len(infeasible),
            "initial_unknown_leaves": len(unknown),
            "gpu_supported_unknown_leaves": len(supported),
            "gpu_excluded_leaves": len(certificates),
            "high_mean_unknown_leaves": len(high_mean),
        })
        total_initial_infeasible += len(infeasible)
        total_initial_unknown += len(unknown)
        total_supported += len(supported)
        total_high_mean += len(high_mean)

    if (total_initial_infeasible, total_initial_unknown, total_supported, total_high_mean) != (172, 1088, 662, 426):
        raise AssertionError("global exceptional coverage counts changed")
    out = {
        "experiment": "p7_exceptional_mod7_tuple_audit",
        "status": "passed_independent_all_orbit_coverage_audit",
        "p": 7,
        "c_H": -1,
        "exceptional_orbit_count": 7,
        "total_mean_allocations": 1260,
        "initial_infeasible_leaves": total_initial_infeasible,
        "initial_unknown_leaves": total_initial_unknown,
        "gpu_supported_unknown_leaves": total_supported,
        "gpu_excluded_leaves": total_supported,
        "gpu_unresolved_leaves": 0,
        "high_mean_unknown_leaves_not_claimed": total_high_mean,
        "selected_dependency_rows_pairwise_disjoint": True,
        "orbit_rows": orbit_rows,
    }
    if args.output is not None:
        atomic_json(args.output, out)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
