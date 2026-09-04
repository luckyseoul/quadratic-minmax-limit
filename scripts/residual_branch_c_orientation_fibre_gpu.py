#!/usr/bin/env python3
"""Exact GPU audit of the sign-gauge fibre of one paired-SDR geometry.

For a localized Mobius half, simultaneously replacing ``(alpha, sign)`` by
``(-alpha, -sign)`` leaves its physical orbit support fixed and reverses all
orbit coefficients.  Thus two paired-SDR witnesses with the same matching and
auxiliary directions can be different orientations of one geometric skeleton.

This program exhausts all ``2^m`` such orientations.  It rejects orientations
that break ternarity at a shared orbit and evaluates the exact transverse
coefficient rows of the resulting common simple graph.  The GPU calculation is
replayed with the repository's integer edge--Radon implementation for every
reported extremizer.

The input skeleton is required to have exactly one isolated orbit shared by
two halves.  Exhaustiveness refers only to sign-gauge reversals of that fixed
skeleton.  Geometries with a clean three-half ``2:1`` overlap, or any other
higher-multiplicity cancellation pattern, are outside this audit.

The compact-triangle l1 bound is necessary, not sufficient.  Consequently a
zero passer count excludes this entire gauge fibre, while a passer still needs
an exact atom decomposition.  Neither outcome alone closes residual (ii).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from e1_gmin_m4_adaptive_mobius_pairing import (  # noqa: E402
    prescribed_auxiliary_assignment_criterion,
)
from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    _functional_value,
    _negative_edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from io_atomic import write_json_atomic  # noqa: E402
from scripts.residual_branch_c_aux_sdr_cpsat import build_chart  # noqa: E402
from scripts.residual_branch_c_top_sparse_sdr_cpsat import (  # noqa: E402
    _fixed_edge,
    _fixed_word_capacity_replay,
    _graph_from_orbit_values,
    _half_orbits,
    _literal_centres,
    _projected_key,
)


def _load_profile(path: Path, profile_index: int) -> tuple[dict, int, dict]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "residual_branch_c_aux_sdr_cpsat_v1":
        raise ValueError("input is not a paired-SDR result")
    profile = payload["profiles"][profile_index]
    rows = [row for row in profile["scale_results"] if "witness" in row]
    if len(rows) != 1:
        raise ValueError("selected profile must contain exactly one witness")
    return profile, int(rows[0]["c"]), rows[0]["witness"]


def _cell_index(p: int) -> tuple[dict[tuple[int, int], int], tuple[tuple[int, int], ...]]:
    pairs = tuple((left, right) for left in range(p) for right in range(left + 1, p))
    return {pair: index for index, pair in enumerate(pairs)}, pairs


def _add_edge_projection(
    array: np.ndarray,
    p: int,
    directions: tuple[tuple[int, int], ...],
    pair_index: dict[tuple[int, int], int],
    edge,
) -> None:
    tau = paley_edge_sign(p, edge)
    for direction_index, direction in enumerate(directions):
        key = _projected_key(p, direction, edge)
        if key[0] == "P":
            continue
        coefficient = paley_direction_sign(p, direction) * tau
        array[direction_index, pair_index[(int(key[1]), int(key[2]))]] += coefficient


def _literal_direction_map(chart, alphas: tuple[int, ...]) -> tuple[dict[int, int], dict[int, int]]:
    directions = projective_functionals(chart.p)
    base_centres = _literal_centres(chart, alphas)
    hard_index_by_coordinate = {
        coordinate: index for index, coordinate in enumerate(chart.hard_coordinates)
    }
    row_to_half: dict[int, int] = {}
    for direction_index, direction in enumerate(directions):
        if paley_direction_sign(chart.p, direction) != 1:
            continue
        evaluation = _functional_value(chart.p, direction, chart.x0)
        normalized = (
            direction[0] * pow(evaluation, -1, chart.p) % chart.p,
            direction[1] * pow(evaluation, -1, chart.p) % chart.p,
        )
        coordinate = chart.direction_by_coordinate.index(normalized)
        row_to_half[direction_index] = hard_index_by_coordinate[coordinate]
    if set(row_to_half) != set(base_centres):
        raise ArithmeticError("hard direction map and literal centres disagree")
    return base_centres, row_to_half


def _orientation_orbit_values(halves: tuple[dict, ...], signs: tuple[int, ...]) -> dict:
    total = Counter()
    for orientation, half in zip(signs, halves, strict=True):
        total.update({orbit: orientation * value for orbit, value in half.items()})
    if any(abs(value) > 1 for value in total.values()):
        raise ArithmeticError("orientation is not ternary")
    return {orbit: value for orbit, value in total.items() if value}


def _cpu_replay(
    chart,
    alphas: tuple[int, ...],
    singleton_signs: tuple[int, ...],
    auxiliaries: tuple[int, ...],
    partners: tuple[int, ...],
    c: int,
    halves: tuple[dict, ...],
    orientations: tuple[int, ...],
    expected_l1: tuple[int, ...],
) -> dict[str, object]:
    p = chart.p
    directions = projective_functionals(p)
    gauged_alphas = tuple(alpha * orientation % p for alpha, orientation in zip(alphas, orientations, strict=True))
    gauged_signs = tuple(sign * orientation for sign, orientation in zip(singleton_signs, orientations, strict=True))
    orbit_values = _orientation_orbit_values(halves, orientations)
    fixed = _fixed_edge(chart, c)
    graph = _graph_from_orbit_values(p, orbit_values, fixed)
    signed_source = {edge: paley_edge_sign(p, edge) for edge in graph}
    image = edge_radon_image(p, signed_source)
    literal_centres = _literal_centres(chart, gauged_alphas)
    row_l1 = []
    parallel_counts = []
    budgets = []
    for direction_index, direction in enumerate(directions):
        direction_sign = paley_direction_sign(p, direction)
        parallel = direction_sign * image.get(("P", direction_index), 0)
        literal = literal_centres.get(direction_index)
        residual = []
        for left in range(p):
            for right in range(left + 1, p):
                value = direction_sign * image.get(
                    ("K", direction_index, left, right), 0
                )
                if literal is not None and literal in (left, right):
                    value += 1
                residual.append(value)
        row_l1.append(sum(abs(value) for value in residual))
        parallel_counts.append(parallel)
        budgets.append(3 * (parallel - 3))
    if tuple(row_l1) != expected_l1:
        raise ArithmeticError("CPU edge--Radon replay disagrees with GPU row l1")
    assignment = prescribed_auxiliary_assignment_criterion(
        p,
        chart.hard_coordinates,
        gauged_alphas,
        auxiliaries,
        partners,
        tuple(range(len(alphas))),
        c,
        gauged_signs,
    )
    fixed_word = _fixed_word_capacity_replay(
        chart, gauged_alphas, c, graph, orbit_values
    )
    graph_bytes = b"".join(
        int(coordinate).to_bytes(2, "little")
        for edge in graph
        for point in edge
        for coordinate in point
    )
    return {
        "orientations": list(orientations),
        "gauged_alphas": list(gauged_alphas),
        "gauged_singleton_signs": list(gauged_signs),
        "orbit_support": len(orbit_values),
        "graph_edge_count": len(graph),
        "graph_sha256": hashlib.sha256(graph_bytes).hexdigest(),
        "row_l1": row_l1,
        "parallel_counts": parallel_counts,
        "row_budgets": budgets,
        "maximum_l1_excess": max(value - budget for value, budget in zip(row_l1, budgets, strict=True)),
        "total_positive_l1_excess": sum(max(0, value - budget) for value, budget in zip(row_l1, budgets, strict=True)),
        "passes_all_l1_rows": all(value <= budget for value, budget in zip(row_l1, budgets, strict=True)),
        "paired_SDR_replay": bool(assignment["pair_coherent_distinct_auxiliary_assignment"]),
        "fixed_word_atom_capacity_replay": fixed_word,
        "exact_integer_edge_Radon_replay": True,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    import cupy as cp

    profile, c, witness = _load_profile(args.input, args.profile_index)
    p = int(json.loads(args.input.read_text())["p"])
    chart = build_chart(p, args.fixed_direction_index)
    alphas = tuple(int(value) % p for value in profile["alphas"])
    singleton_signs = tuple(int(value) for value in witness["signs"])
    auxiliaries = tuple(int(value) for value in witness["auxiliary_coordinates"])
    partners = tuple(int(value) for value in witness["partners"])
    m = len(alphas)
    if m > 31:
        raise ValueError("exact uint32 orientation enumeration supports at most 31 halves")

    halves = tuple(
        _half_orbits(
            chart,
            chart.hard_coordinates[index],
            alphas[index],
            auxiliaries[index],
            singleton_signs[index],
            c,
        )
        for index in range(m)
    )
    gauge_reversal = []
    for index, half in enumerate(halves):
        flipped = _half_orbits(
            chart,
            chart.hard_coordinates[index],
            -alphas[index] % p,
            auxiliaries[index],
            -singleton_signs[index],
            c,
        )
        gauge_reversal.append(
            set(flipped) == set(half)
            and all(flipped[orbit] == -value for orbit, value in half.items())
        )
    if not all(gauge_reversal):
        raise ArithmeticError("the sign gauge failed to reverse a Mobius half")

    occurrences: defaultdict[object, list[tuple[int, int]]] = defaultdict(list)
    for half_index, half in enumerate(halves):
        for orbit, value in half.items():
            occurrences[orbit].append((half_index, value))
    overlaps = {orbit: rows for orbit, rows in occurrences.items() if len(rows) > 1}
    if len(overlaps) != 1 or any(len(rows) != 2 for rows in overlaps.values()):
        raise ValueError("this exact top-fibre audit requires one two-half overlap")
    overlap_orbit, overlap_rows = next(iter(overlaps.items()))
    if sum(value for _index, value in overlap_rows) != 0:
        raise ValueError("the witness overlap is not the required cancellation")

    directions = projective_functionals(p)
    pair_index, pairs = _cell_index(p)
    width = len(pairs)
    fixed_projection = np.zeros((p + 1, width), dtype=np.int16)
    fixed = _fixed_edge(chart, c)
    _add_edge_projection(fixed_projection, p, directions, pair_index, fixed)
    base_centres, row_to_half = _literal_direction_map(chart, alphas)

    half_projection = np.zeros((m, 2, p + 1, width), dtype=np.int16)
    for half_index, half in enumerate(halves):
        for orientation_slot, orientation in enumerate((1, -1)):
            array = half_projection[half_index, orientation_slot]
            for orbit, coefficient in half.items():
                if orbit == overlap_orbit:
                    continue
                signed_difference = paley_edge_sign(p, orbit) * orientation * coefficient
                edge = orbit if signed_difference == 1 else _negative_edge(p, orbit)
                _add_edge_projection(array, p, directions, pair_index, edge)
    for direction_index, half_index in row_to_half.items():
        centre = base_centres[direction_index]
        for orientation_slot, orientation in enumerate((1, -1)):
            literal = centre if orientation == 1 else -centre % p
            for other in range(p):
                if other == literal:
                    continue
                left, right = sorted((literal, other))
                half_projection[
                    half_index, orientation_slot, direction_index, pair_index[(left, right)]
                ] += 1

    base_orbits = _orientation_orbit_values(halves, (1,) * m)
    base_graph = _graph_from_orbit_values(p, base_orbits, fixed)
    base_image = edge_radon_image(
        p, {edge: paley_edge_sign(p, edge) for edge in base_graph}
    )
    parallel_counts = np.asarray(
        [
            paley_direction_sign(p, direction)
            * base_image.get(("P", direction_index), 0)
            for direction_index, direction in enumerate(directions)
        ],
        dtype=np.int32,
    )
    budgets = 3 * (parallel_counts - 3)
    expected_edges = m * (p - 1) - 1
    if len(base_graph) != expected_edges or np.any(budgets < 0):
        raise ArithmeticError("paired-SDR geometry is not the top branch-C graph")

    state_ids = np.arange(1 << m, dtype=np.uint32)
    bits = ((state_ids[:, None] >> np.arange(m, dtype=np.uint32)) & 1).astype(np.int8)
    signs_matrix = 1 - 2 * bits
    first, first_coefficient = overlap_rows[0]
    second, second_coefficient = overlap_rows[1]
    valid = (
        signs_matrix[:, first] * first_coefficient
        + signs_matrix[:, second] * second_coefficient
        == 0
    )
    valid_ids = state_ids[valid]
    valid_bits = bits[valid]

    base = fixed_projection.astype(np.int32)
    base += half_projection[:, 0].sum(axis=0, dtype=np.int32)
    delta = (
        half_projection[:, 1].astype(np.int32)
        - half_projection[:, 0].astype(np.int32)
    ).reshape(m, -1)
    flat_base = base.reshape(-1)
    device = cp.cuda.Device(args.device)
    device.use()
    gpu_delta = cp.asarray(delta, dtype=cp.float32)
    gpu_base = cp.asarray(flat_base, dtype=cp.float32)
    gpu_budgets = cp.asarray(budgets, dtype=cp.int32)

    started = time.monotonic()
    passer_ids: list[int] = []
    best_key: tuple[int, int, int] | None = None
    best_id = -1
    best_l1: tuple[int, ...] | None = None
    evaluated = 0
    for offset in range(0, len(valid_ids), args.batch):
        batch_bits = valid_bits[offset : offset + args.batch]
        matrix = cp.asarray(batch_bits, dtype=cp.float32)
        values = gpu_base[None, :] + matrix @ gpu_delta
        values = cp.rint(values).astype(cp.int16).reshape(-1, p + 1, width)
        row_l1 = cp.sum(cp.abs(values).astype(cp.int32), axis=2, dtype=cp.int32)
        excess = row_l1 - gpu_budgets[None, :]
        maximum = cp.max(excess, axis=1)
        positive_total = cp.sum(cp.maximum(excess, 0), axis=1, dtype=cp.int32)
        maximum_host = cp.asnumpy(maximum)
        total_host = cp.asnumpy(positive_total)
        l1_host = cp.asnumpy(row_l1)
        ids_host = valid_ids[offset : offset + len(batch_bits)]
        for local_index, state_id in enumerate(ids_host):
            key = (int(maximum_host[local_index]), int(total_host[local_index]), int(state_id))
            if best_key is None or key < best_key:
                best_key = key
                best_id = int(state_id)
                best_l1 = tuple(int(value) for value in l1_host[local_index])
            if maximum_host[local_index] <= 0 and len(passer_ids) < args.max_passers:
                passer_ids.append(int(state_id))
        evaluated += len(batch_bits)
    cp.cuda.runtime.deviceSynchronize()
    elapsed = time.monotonic() - started
    properties = cp.cuda.runtime.getDeviceProperties(args.device)
    device_name = properties["name"]
    if isinstance(device_name, bytes):
        device_name = device_name.decode()
    if best_l1 is None:
        raise ArithmeticError("no valid orientation was evaluated")

    def orientation_for(state_id: int) -> tuple[int, ...]:
        return tuple(-1 if (state_id >> index) & 1 else 1 for index in range(m))

    best_replay = _cpu_replay(
        chart,
        alphas,
        singleton_signs,
        auxiliaries,
        partners,
        c,
        halves,
        orientation_for(best_id),
        best_l1,
    )
    passer_replays = []
    for state_id in passer_ids:
        orientations = orientation_for(state_id)
        orbit_values = _orientation_orbit_values(halves, orientations)
        graph = _graph_from_orbit_values(p, orbit_values, fixed)
        image = edge_radon_image(
            p, {edge: paley_edge_sign(p, edge) for edge in graph}
        )
        gauged_alphas = tuple(
            alpha * orientation % p
            for alpha, orientation in zip(alphas, orientations, strict=True)
        )
        literal_centres = _literal_centres(chart, gauged_alphas)
        exact_l1 = []
        for direction_index, direction in enumerate(directions):
            direction_sign = paley_direction_sign(p, direction)
            literal = literal_centres.get(direction_index)
            values = []
            for left, right in pairs:
                value = direction_sign * image.get(
                    ("K", direction_index, left, right), 0
                )
                if literal is not None and literal in (left, right):
                    value += 1
                values.append(value)
            exact_l1.append(sum(abs(value) for value in values))
        passer_replays.append(
            _cpu_replay(
                chart,
                alphas,
                singleton_signs,
                auxiliaries,
                partners,
                c,
                halves,
                orientations,
                tuple(exact_l1),
            )
        )

    input_hash = hashlib.sha256(args.input.read_bytes()).hexdigest()
    return {
        "schema": "residual_branch_c_orientation_fibre_gpu_v1",
        "classification": "exact exhaustive sign-gauge fibre audit of one paired-SDR geometry",
        "overlap_geometry_scope": (
            "one fixed skeleton with one isolated two-half opposite overlap; "
            "no three-half 2:1 or higher-multiplicity overlap geometry"
        ),
        "higher_multiplicity_overlap_geometries_included": False,
        "scope_warning": (
            "A zero passer count excludes only the sign-gauge fibre of this fixed "
            "two-half-overlap skeleton; a passer still needs exact atom decomposition. "
            "Residual ii is not closed by this audit."
        ),
        "host": platform.node(),
        "gpu_backend": "cupy",
        "gpu_device": str(device_name),
        "input": str(args.input),
        "input_sha256": input_hash,
        "profile_index": args.profile_index,
        "p": p,
        "c": c,
        "half_count": m,
        "raw_orientation_count": 1 << m,
        "ternary_orientation_count": int(len(valid_ids)),
        "evaluated_orientation_count": evaluated,
        "gauge_reversal_verified_for_every_half": all(gauge_reversal),
        "shared_orbit_count": len(overlaps),
        "shared_orbit_half_indices": [first, second],
        "shared_orbit_coefficients": [first_coefficient, second_coefficient],
        "fixed_graph_edge_count": expected_edges,
        "row_budgets": budgets.tolist(),
        "l1_gate_passer_count_capped": len(passer_ids),
        "l1_gate_passer_count_is_capped": len(passer_ids) == args.max_passers,
        "best_state_id": best_id,
        "best_replay": best_replay,
        "passer_replays": passer_replays,
        "gpu_seconds": elapsed,
        "exact_float32_note": "Each GEMM entry is a sum of at most m small integers and is exactly representable; every reported state is independently replayed over integers.",
        "full_residual_ii_closed": False,
    }


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--profile-index", type=int, default=0)
    parser.add_argument("--fixed-direction-index", type=int)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--batch", type=int, default=1024)
    parser.add_argument("--max-passers", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    result = run(args)
    write_json_atomic(args.output, result)
    print(
        f"wrote {args.output}: evaluated={result['evaluated_orientation_count']} "
        f"passers={result['l1_gate_passer_count_capped']} "
        f"best_max_excess={result['best_replay']['maximum_l1_excess']} "
        f"gpu_seconds={result['gpu_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
