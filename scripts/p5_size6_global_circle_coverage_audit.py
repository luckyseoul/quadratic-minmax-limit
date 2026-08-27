#!/usr/bin/env python3
"""Audit the global reduction of p=5 size-six boundaries to six circles.

The audit deliberately separates two claims:

1. the four boundary catalogs and the exact coarse SCIP batches reduce every
   p=5 size-six boundary to six residual signed-symmetry classes;
2. independent layered audits close each of those six classes.

Raw catalogs and solver batches are too large for the repository.  This
program records their hashes and compact summaries while rebuilding the
finite catalogs and all symmetry identities from source code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from p5_infinity_full_multiplier_orbits import (  # noqa: E402
    coarsen,
    field_symmetries,
)
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402
from residual_fixed_size_boundary_orbits import classify  # noqa: E402


NOINF_RANGES = [
    *(range(start, start + 250) for start in range(0, 5500, 250)),
    range(5500, 6766),
]
PLUSINF_RANGES = [range(0, 450), range(450, 905)]
MINUSINF_RANGES = [range(0, 1150)]
EXPECTED_NOINF_RESIDUALS = (881, 2529, 3032, 4731, 4939)
EXPECTED_MINUSINF_RESIDUALS = (0, 1144)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def without_elapsed(payload: dict) -> dict:
    result = dict(payload)
    result.pop("elapsed_seconds", None)
    return result


def assert_catalog_equal(recorded: dict, rebuilt: dict, label: str) -> None:
    recorded_canonical = json.dumps(
        without_elapsed(recorded), sort_keys=True, separators=(",", ":")
    )
    rebuilt_canonical = json.dumps(
        without_elapsed(rebuilt), sort_keys=True, separators=(",", ":")
    )
    if recorded_canonical != rebuilt_canonical:
        raise AssertionError(f"rebuilt {label} catalog differs from the recording")


def artifact(path: Path) -> dict:
    return {"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}


def reconstruct_orbits(
    catalog: dict, permutations: tuple[tuple[int, ...], ...]
) -> tuple[dict[tuple[int, ...], int], list[set[tuple[int, ...]]]]:
    member_to_orbit: dict[tuple[int, ...], int] = {}
    members_by_orbit = []
    for index, row in enumerate(catalog["orbits"]):
        representative = tuple(int(value) for value in row["representative_finite_field"])
        members = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in permutations
        }
        if len(members) != int(row["size"]):
            raise AssertionError(f"orbit {index} has an incorrect recorded size")
        for member in members:
            if member in member_to_orbit:
                raise AssertionError("recorded catalog orbits overlap")
            member_to_orbit[member] = index
        members_by_orbit.append(members)
    if len(member_to_orbit) != int(catalog["candidate_boundaries"]):
        raise AssertionError("catalog orbit union misses candidate boundaries")
    if sum(int(row["size"]) for row in catalog["orbits"]) != int(
        catalog["orbit_size_sum"]
    ):
        raise AssertionError("catalog orbit-size sum is internally inconsistent")
    if len(catalog["orbits"]) != int(catalog["orbit_count"]):
        raise AssertionError("catalog orbit count is internally inconsistent")
    return member_to_orbit, members_by_orbit


def signed_map_identity(
    vertex_permutation: tuple[int, ...], switching: tuple[int, ...]
) -> dict:
    """Verify a signed conference symmetry and its action on both shells."""
    if sorted(vertex_permutation) != list(range(26)):
        raise AssertionError("vertex map is not a permutation")
    if len(switching) != 26 or set(switching) - {-1, 1}:
        raise AssertionError("invalid switching vector")
    C = np.rint(paley_conference_prime_power(5)).astype(np.int8)
    factors = set()
    for a in range(26):
        for b in range(a + 1, 26):
            factors.add(
                int(switching[a])
                * int(C[vertex_permutation[a], vertex_permutation[b]])
                * int(switching[b])
                * int(C[a, b])
            )
    if len(factors) != 1:
        raise AssertionError("signed map has no uniform conference factor")
    global_factor = factors.pop()
    if global_factor not in (-1, 1):
        raise AssertionError("bad signed conference factor")

    full_shells = geometry(5, "full")["shells"]
    targets = {
        eps: {tuple(int(value) for value in row) for row in full_shells[eps]}
        for eps in (-1, 1)
    }
    shell_map = {}
    for eps in (-1, 1):
        transformed = set()
        for row in full_shells[eps]:
            image = np.empty(26, dtype=np.int8)
            image[np.asarray(vertex_permutation, dtype=np.int16)] = (
                np.asarray(switching, dtype=np.int8) * row
            )
            transformed.add(tuple(int(value) for value in image))
        matches = [target for target in (-1, 1) if transformed == targets[target]]
        if len(matches) != 1:
            raise AssertionError("signed map does not permute the two complete shells")
        shell_map[str(eps)] = matches[0]
    if set(shell_map.values()) != {-1, 1}:
        raise AssertionError("signed shell action is not bijective")
    return {
        "global_conference_factor": global_factor,
        "complete_shell_map": shell_map,
        "fixes_distinguished_edge_setwise": {
            vertex_permutation[0], vertex_permutation[1]
        }
        == {0, 1},
    }


def product_factor(
    identity: dict, switching: tuple[int, ...], boundary: tuple[int, ...]
) -> int:
    factor = int(identity["global_conference_factor"]) ** 21
    for vertex in boundary:
        factor *= int(switching[vertex])
    return factor


def nonsquare_multiplier_audit(
    noinf_plus: dict,
    noinf_minus: dict,
    square: tuple[tuple[int, ...], ...],
    nonsquare_permutation: tuple[int, ...],
) -> tuple[dict, dict[tuple[int, ...], int], list[set[tuple[int, ...]]]]:
    plus_members, plus_orbits = reconstruct_orbits(noinf_plus, square)
    minus_members, minus_orbits = reconstruct_orbits(noinf_minus, square)
    vertex = (0, *(value + 1 for value in nonsquare_permutation))
    switching = (-1, *(1 for _value in range(25)))
    identity = signed_map_identity(vertex, switching)
    if identity["global_conference_factor"] != -1:
        raise AssertionError("nonsquare multiplier is not the expected anti-isometry")

    mapped_members = {
        tuple(sorted(nonsquare_permutation[value] for value in boundary))
        for boundary in plus_members
    }
    if mapped_members != set(minus_members):
        raise AssertionError("c_H=+1 no-infinity survivors do not map onto c_H=-1")
    transfer = {}
    for index, members in enumerate(plus_orbits):
        targets = {
            minus_members[
                tuple(sorted(nonsquare_permutation[value] for value in boundary))
            ]
            for boundary in members
        }
        if len(targets) != 1:
            raise AssertionError("nonsquare map splits a square-semilinear orbit")
        transfer[index] = targets.pop()
    if set(transfer.values()) != set(range(len(minus_orbits))):
        raise AssertionError("nonsquare orbit transfer is not bijective")
    representative_factor = {
        product_factor(
            identity,
            switching,
            tuple(int(value) for value in row["representative_vertices"]),
        )
        for row in noinf_plus["orbits"]
    }
    if representative_factor != {-1}:
        raise AssertionError("nonsquare no-infinity map does not flip c_H exactly")
    return (
        {
            **identity,
            "source_orbits": len(plus_orbits),
            "target_orbits": len(minus_orbits),
            "mapped_survivors": len(mapped_members),
            "paley_product_factor": -1,
            "c_H_action": "+1 to -1",
            "orbit_transfer_is_bijective": True,
            "orbit_transfer_sha256": hashlib.sha256(
                json.dumps(transfer, sort_keys=True).encode()
            ).hexdigest(),
        },
        minus_members,
        minus_orbits,
    )


def inverse(value: int, mul) -> int:
    if value == 0:
        raise ZeroDivisionError("finite-field inverse of zero")
    result, base, exponent = 1, value, 23
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent >>= 1
    return result


def infinity_to_noinf_audit(
    minusinf: dict,
    noinf_member_to_orbit: dict[tuple[int, ...], int],
    alpha: int,
    source_index: int,
) -> dict:
    _q2, mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(5)
    sigma = int(chi(alpha))
    vertex = [1, 0]
    vertex.extend(1 + mul(alpha, inverse(value, mul)) for value in range(1, 25))
    switching = tuple(
        sigma
        if old_vertex == 1
        else int(chi(old_vertex - 1))
        if old_vertex > 1
        else 1
        for old_vertex in range(26)
    )
    identity = signed_map_identity(tuple(vertex), switching)
    if not identity["fixes_distinguished_edge_setwise"]:
        raise AssertionError("signed inversion moved the distinguished edge")
    source_boundary = tuple(
        int(value)
        for value in minusinf["orbits"][source_index]["representative_vertices"]
    )
    image_boundary = tuple(sorted(vertex[value] for value in source_boundary))
    if 0 in image_boundary:
        raise AssertionError("selected infinity orbit did not map to no-infinity")
    image_finite = tuple(value - 1 for value in image_boundary)
    target_index = noinf_member_to_orbit.get(image_finite)
    if target_index is None:
        raise AssertionError("signed-inversion image is absent from no-infinity catalog")
    factor = product_factor(identity, switching, source_boundary)
    if factor != 1:
        raise AssertionError("signed inversion failed to preserve c_H")
    return {
        **identity,
        "alpha": alpha,
        "chi_alpha": sigma,
        "source_kind": "c_H=-1, infinity present",
        "source_orbit_index": source_index,
        "source_boundary": list(source_boundary),
        "image_boundary": list(image_boundary),
        "target_kind": "c_H=-1, no infinity",
        "target_orbit_index": target_index,
        "paley_product_factor": factor,
        "c_H_is_preserved": factor == 1,
    }


def batch_filename(prefix: str, scope: range) -> str:
    return f"p5_size6_exact_{prefix}_{scope.start}_{scope.stop}.json"


def audit_batches(
    catalog_path: Path,
    catalog: dict,
    workspace: Path,
    prefix: str,
    ranges: list[range],
) -> dict:
    catalog_sha = sha256(catalog_path)
    expected_indices = set(range(int(catalog["orbit_count"])))
    covered = set()
    unknown = []
    files = []
    totals = Counter()
    for scope in ranges:
        path = workspace / batch_filename(prefix, scope)
        payload = load(path)
        if payload.get("source_sha256") != catalog_sha:
            raise AssertionError(f"source hash mismatch in {path}")
        if int(payload["start_orbit"]) != scope.start or int(
            payload["stop_orbit"]
        ) != scope.stop:
            raise AssertionError(f"range mismatch in {path}")
        rows = payload["rows"]
        row_indices = [int(row["orbit_index"]) for row in rows]
        if row_indices != list(scope):
            raise AssertionError(f"row coverage or order mismatch in {path}")
        if covered & set(row_indices):
            raise AssertionError("exact batches overlap")
        covered.update(row_indices)
        statuses = Counter(str(row["solver_status"]) for row in rows)
        if dict(statuses) != payload["status_counts"]:
            raise AssertionError(f"status histogram mismatch in {path}")
        closed_count = 0
        unknown_count = 0
        for row in rows:
            index = int(row["orbit_index"])
            orbit = catalog["orbits"][index]
            if row.get("source_sha256") != catalog_sha:
                raise AssertionError("row source hash mismatch")
            if tuple(int(value) for value in row["boundary"]) != tuple(
                int(value) for value in orbit["representative_vertices"]
            ):
                raise AssertionError("batch row has the wrong boundary")
            if int(row["c_H"]) != int(catalog["c_H"]):
                raise AssertionError("batch row has the wrong Paley-product sign")
            if row.get("feasible") is True:
                raise AssertionError("coarse exact batch contains a feasible witness")
            if (
                row.get("finite_infeasibility_certificate") is True
                and row.get("feasible") is False
            ):
                closed_count += 1
            elif (
                row.get("solver_status") == "UNKNOWN"
                and row.get("finite_infeasibility_certificate") is False
                and row.get("feasible") is False
            ):
                unknown.append(index)
                unknown_count += 1
            else:
                raise AssertionError("batch row has an unclassified solver outcome")
        if (
            int(payload["scope_orbits"]) != len(scope)
            or int(payload["completed"]) != len(scope)
            or int(payload["infeasible"]) != closed_count
            or int(payload["unknown"]) != unknown_count
            or int(payload["feasible"]) != 0
            or bool(payload["all_infeasible"]) != (unknown_count == 0)
        ):
            raise AssertionError(f"summary count mismatch in {path}")
        totals.update({"rows": len(rows), "infeasible": closed_count, "unknown": unknown_count})
        files.append(
            {
                **artifact(path),
                "start_orbit": scope.start,
                "stop_orbit": scope.stop,
                "infeasible": closed_count,
                "unknown": unknown_count,
            }
        )
    if covered != expected_indices:
        raise AssertionError("exact batches do not cover the whole catalog")
    return {
        "catalog": artifact(catalog_path),
        "catalog_orbits": len(expected_indices),
        "covered_orbits": len(covered),
        "infeasible": totals["infeasible"],
        "unknown": totals["unknown"],
        "unknown_orbit_indices": unknown,
        "batch_count": len(files),
        "batches": files,
    }


def audit_class_closures(
    paths: list[Path],
    noinf_path: Path,
    noinf: dict,
    minusinf_path: Path,
    minusinf: dict,
) -> dict:
    noinf_sha = sha256(noinf_path)
    minusinf_sha = sha256(minusinf_path)
    expected = {
        *((noinf_sha, index) for index in EXPECTED_NOINF_RESIDUALS),
        (minusinf_sha, 0),
    }
    observed = set()
    records = []
    for path in paths:
        payload = load(path)
        scope = payload.get("scope", {})
        key = (str(scope.get("source_sha256")), int(scope.get("orbit_index", -1)))
        if key not in expected:
            raise AssertionError(f"unexpected class audit {path}")
        if key in observed:
            raise AssertionError("duplicate class closure audit")
        if payload.get("proved") is not True:
            raise AssertionError(f"class audit did not prove closure: {path}")
        if [int(value) for value in payload.get("covered_internal_edge_counts", [])] != list(
            range(16)
        ):
            raise AssertionError("class audit does not cover all internal edge counts")
        catalog = noinf if key[0] == noinf_sha else minusinf
        if tuple(int(value) for value in scope["boundary"]) != tuple(
            int(value) for value in catalog["orbits"][key[1]]["representative_vertices"]
        ):
            raise AssertionError("class audit boundary does not match its catalog")
        observed.add(key)
        records.append(
            {
                **artifact(path),
                "source_sha256": key[0],
                "orbit_index": key[1],
                "boundary": scope["boundary"],
                "used_artifact_count": int(payload["used_artifact_count"]),
            }
        )
    return {
        "expected_class_count": len(expected),
        "proved_class_count": len(observed),
        "all_residual_classes_closed": observed == expected,
        "missing": [
            {"source_sha256": source_sha, "orbit_index": index}
            for source_sha, index in sorted(expected - observed)
        ],
        "audits": sorted(records, key=lambda row: (row["source_sha256"], row["orbit_index"])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--noinf-minus-source", type=Path, required=True)
    parser.add_argument("--noinf-plus-source", type=Path, required=True)
    parser.add_argument("--minusinf-square-source", type=Path, required=True)
    parser.add_argument("--plusinf-square-source", type=Path, required=True)
    parser.add_argument("--minusinf-source", type=Path, required=True)
    parser.add_argument("--plusinf-source", type=Path, required=True)
    parser.add_argument("--batch-workspace", type=Path, required=True)
    parser.add_argument("--class-audit", type=Path, action="append", default=[])
    parser.add_argument("--skip-rebuild", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    noinf_minus = load(args.noinf_minus_source)
    noinf_plus = load(args.noinf_plus_source)
    minusinf_square = load(args.minusinf_square_source)
    plusinf_square = load(args.plusinf_square_source)
    minusinf = load(args.minusinf_source)
    plusinf = load(args.plusinf_source)
    for payload, c_h, infinity, label in (
        (noinf_minus, -1, 0, "noinf-minus"),
        (noinf_plus, 1, 0, "noinf-plus"),
        (minusinf_square, -1, 1, "minusinf-square"),
        (plusinf_square, 1, 1, "plusinf-square"),
    ):
        if (
            payload.get("experiment") != "residual_fixed_size_boundary_orbits"
            or int(payload["p"]) != 5
            or int(payload["boundary_size"]) != 6
            or int(payload["c_H"]) != c_h
            or int(payload["infinity_value"]) != infinity
        ):
            raise AssertionError(f"wrong source catalog supplied for {label}")

    rebuilt = not args.skip_rebuild
    if rebuilt:
        assert_catalog_equal(noinf_minus, classify(5, -1, 6, 0), "noinf-minus")
        assert_catalog_equal(noinf_plus, classify(5, 1, 6, 0), "noinf-plus")
        assert_catalog_equal(minusinf_square, classify(5, -1, 6, 1), "minusinf-square")
        assert_catalog_equal(plusinf_square, classify(5, 1, 6, 1), "plusinf-square")
    assert_catalog_equal(minusinf, coarsen(args.minusinf_square_source), "minusinf")
    assert_catalog_equal(plusinf, coarsen(args.plusinf_square_source), "plusinf")

    square, full, nonsquare_permutation, nonsquare = field_symmetries()
    if len(square) != 24 or len(full) != 48:
        raise AssertionError("unexpected p=5 symmetry-group sizes")
    noinf_transfer, noinf_members, _noinf_orbits = nonsquare_multiplier_audit(
        noinf_plus, noinf_minus, square, nonsquare_permutation
    )
    noinf_transfer["nonsquare_multiplier"] = nonsquare
    reconstruct_orbits(minusinf, full)
    reconstruct_orbits(plusinf, full)

    noinf_batches = audit_batches(
        args.noinf_minus_source,
        noinf_minus,
        args.batch_workspace,
        "noinf",
        NOINF_RANGES,
    )
    plusinf_batches = audit_batches(
        args.plusinf_source,
        plusinf,
        args.batch_workspace,
        "plusinf",
        PLUSINF_RANGES,
    )
    minusinf_batches = audit_batches(
        args.minusinf_source,
        minusinf,
        args.batch_workspace,
        "minusinf",
        MINUSINF_RANGES,
    )
    if tuple(noinf_batches["unknown_orbit_indices"]) != EXPECTED_NOINF_RESIDUALS:
        raise AssertionError("unexpected no-infinity residual orbit set")
    if plusinf_batches["unknown_orbit_indices"]:
        raise AssertionError("c_H=+1 infinity catalog has a residual orbit")
    if tuple(minusinf_batches["unknown_orbit_indices"]) != EXPECTED_MINUSINF_RESIDUALS:
        raise AssertionError("unexpected c_H=-1 infinity residual orbit set")

    inversion = infinity_to_noinf_audit(minusinf, noinf_members, nonsquare, 1144)
    if inversion["target_orbit_index"] != 881:
        raise AssertionError("infinity orbit 1144 did not transfer to noinf orbit 881")
    residual_classes = [
        {
            "source": args.minusinf_source.name,
            "source_sha256": sha256(args.minusinf_source),
            "orbit_index": 0,
            "boundary": minusinf["orbits"][0]["representative_vertices"],
        },
        *(
            {
                "source": args.noinf_minus_source.name,
                "source_sha256": sha256(args.noinf_minus_source),
                "orbit_index": index,
                "boundary": noinf_minus["orbits"][index]["representative_vertices"],
            }
            for index in EXPECTED_NOINF_RESIDUALS
        ),
    ]
    closures = audit_class_closures(
        args.class_audit,
        args.noinf_minus_source,
        noinf_minus,
        args.minusinf_source,
        minusinf,
    )
    selection_proved = True
    proved = selection_proved and closures["all_residual_classes_closed"]
    result = {
        "experiment": "p5_size6_global_circle_coverage_audit",
        "status": "exact_global_catalog_and_layered_class_audit",
        "proved": proved,
        "selection_reduction_proved": selection_proved,
        "catalogs_rebuilt_from_definitions": rebuilt,
        "scope": {
            "p": 5,
            "boundary_size": 6,
            "distinguished_edge": [0, 1],
            "signs": [-1, 1],
            "infinity_bits": [0, 1],
        },
        "catalog_artifacts": [
            artifact(path)
            for path in (
                args.noinf_minus_source,
                args.noinf_plus_source,
                args.minusinf_square_source,
                args.plusinf_square_source,
                args.minusinf_source,
                args.plusinf_source,
            )
        ],
        "noinf_sign_transfer": noinf_transfer,
        "infinity_1144_transfer": inversion,
        "coarse_exact_batches": {
            "noinf_c_minus": noinf_batches,
            "plusinf_c_plus": plusinf_batches,
            "minusinf_c_minus": minusinf_batches,
        },
        "coarse_residual_count_before_cross_infinity_symmetry": 7,
        "residual_class_count": len(residual_classes),
        "residual_classes": residual_classes,
        "class_closures": closures,
        "solver_claim_limit": (
            "Batch and layered audits verify finite coverage and recorded exact-solver "
            "infeasibility statuses; they are not independently checkable SCIP proof logs."
        ),
        "elapsed_seconds": time.time() - started,
    }
    atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "selection_reduction_proved": selection_proved,
                "residual_class_count": len(residual_classes),
                "closed_class_count": closures["proved_class_count"],
                "proved": proved,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
