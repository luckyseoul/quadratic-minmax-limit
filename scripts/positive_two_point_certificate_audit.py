#!/usr/bin/env python3
"""Audit the complete small-prime positive two-point certificate bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from positive_two_point_additive_cpsat import (  # noqa: E402
    allowed_k0_values,
    exact_l1_star_profiles,
)


DIRECT_CASES = {
    5: (0, 1, 2, 3, 4, 5, 8),
    7: (2, 4, 8),
    11: (2, 8),
    13: (1, 8),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def assert_infeasible(row: dict) -> None:
    assert row["solver_status"] == "INFEASIBLE"
    assert row["feasible"] is False
    assert row["finite_infeasibility_certificate"] is True


def audit(directory: Path) -> dict:
    expected_k0 = {
        5: [0, 1, 2, 3, 4, 5, 8],
        7: [0, 2, 4, 8],
        11: [0, 2, 8],
        13: [0, 1, 8],
    }
    assert {p: allowed_k0_values(p) for p in expected_k0} == expected_k0

    direct_rows = []
    paths = []
    for p, k0_values in DIRECT_CASES.items():
        for k0 in k0_values:
            path = directory / f"positive_final_p{p}_k{k0}.json"
            row = load(path)
            assert (row["p"], row["k0"]) == (p, k0)
            assert_infeasible(row)
            direct_rows.append(row)
            paths.append(path)

    # Exact directional l1 profiles remove the populated direction required
    # at p=7,k0=4 and p=11,k0=2, while p=13,k0=1 permits only kd=1 in all
    # fourteen directions and therefore cannot have total multiplicity seven.
    assert exact_l1_star_profiles(7, 4, 1) == []
    assert exact_l1_star_profiles(11, 2, 1) == []
    assert exact_l1_star_profiles(13, 1, 0) == []
    assert len(exact_l1_star_profiles(13, 1, 1)) == 3

    classifications = {}
    for kind, name in ((-1, "minus"), (1, "plus")):
        path = directory / f"p7_positive_star_orbits_type_{name}.json"
        row = load(path)
        assert row["populated_type"] == kind
        assert row["generated_candidates"] == 238644
        assert row["survivor_count"] == 2250
        assert row["stabilizer_size"] == 48
        assert row["orbit_count"] == len(row["orbits"]) == 56
        assert sum(orbit["size"] for orbit in row["orbits"]) == 2250
        assert all(orbit["contains_zero"] for orbit in row["orbits"])
        classifications[kind] = row
        paths.append(path)

    orbit_path = directory / "p7_positive_orbit_certificate.json"
    orbit_certificate = load(orbit_path)
    assert orbit_certificate["status"] == "complete_all_infeasible"
    assert orbit_certificate["case_count"] == orbit_certificate["completed_count"] == 112
    assert orbit_certificate["infeasible_count"] == 112
    assert orbit_certificate["unknown_count"] == orbit_certificate["feasible_count"] == 0
    assert len(orbit_certificate["rows"]) == 112
    by_index = {row["case_index"]: row for row in orbit_certificate["rows"]}
    assert set(by_index) == set(range(112))
    records = [
        (kind, orbit_index, orbit)
        for kind in (-1, 1)
        for orbit_index, orbit in enumerate(classifications[kind]["orbits"])
    ]
    for index, (kind, orbit_index, orbit) in enumerate(records):
        row = by_index[index]
        assert_infeasible(row)
        assert row["populated_type"] == kind
        assert row["star_orbit_index"] == orbit_index
        assert row["star_orbit_size"] == orbit["size"]
        assert row["star"] == orbit["representative"]
    paths.append(orbit_path)

    all_one_names = (
        "p7_positive_all_kd1_star0_0_final.json",
        "p7_positive_all_kd1_square.json",
        "p7_positive_all_kd1_all_nonsquare.json",
    )
    all_one = [load(directory / name) for name in all_one_names]
    for row in all_one:
        assert_infeasible(row)
        assert row["p"] == 7 and row["k0"] == 0
        assert row["fixed_direction_multiplicities"] == [1] * 8
    assert all_one[0]["fixed_star_zero"] == 0
    assert all_one[1]["fixed_star_in"] == [0, 1]
    assert all_one[1]["fixed_star_out"] == []
    q2, _mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(7)
    squares = [u for u in range(1, q2) if chi(u) == 1]
    assert chi(1) == 1 and chi(8) == -1
    assert all_one[2]["fixed_star_in"] == [0, 8]
    assert all_one[2]["fixed_star_out"] == squares
    paths.extend(directory / name for name in all_one_names)

    # For k0=0 and p>=7, one unpopulated direction saturates the l1 bound
    # and forces every finite edge to have the opposite type.  Thus all
    # populated directions fit into one type.  At p=11,13 their required
    # total multiplicity eight exceeds the type sizes six and seven.  At
    # p=7 this leaves only the certified four-by-kd=2 type split; if there
    # is no unpopulated direction, all eight multiplicities equal one and
    # the three normalized all-one certificates cover star[0]=0 and 1.
    assert (11 + 1) // 2 < 8 and (13 + 1) // 2 < 8
    assert (7 + 1) // 2 * 2 == 8

    hashes = {path.name: sha256(path) for path in sorted(paths)}
    return {
        "experiment": "positive_two_point_certificate_audit",
        "proved": True,
        "arithmetic_k0_values": {str(p): values for p, values in expected_k0.items()},
        "direct_infeasible_cases": len(direct_rows),
        "p7_rigid_star_candidates_per_type": 238644,
        "p7_rigid_star_survivors_per_type": 2250,
        "p7_rigid_star_orbits_per_type": 56,
        "p7_rigid_infeasible_orbits": 112,
        "p7_all_one_normalized_infeasible_cases": 3,
        "unknown_count": 0,
        "feasible_count": 0,
        "files_sha256": hashes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.directory)
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
