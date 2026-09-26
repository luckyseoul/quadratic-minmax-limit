#!/usr/bin/env python3
"""Exact, bounded one-vertex certificate using deficit-aware transform sections.

This is a sufficient construction, not an extension optimizer or a convergence
proof. A stopped transform never proves that an extension row does not exist.
All geometry uses fractions. Substantial runs belong on an external compute
host; enumeration is opt-in and can be saved for subsequent certificate runs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import gcd, isqrt
from pathlib import Path
import sys


@dataclass(frozen=True)
class Slab:
    normal: tuple[int, ...]
    mass: Fraction = Fraction(1)
    reserve: Fraction = Fraction(0)
    cost: Fraction = Fraction(0)
    phases: int = 0  # 1 = positive, 2 = negative, 3 = mixed ancestry

    def width(self, buffer: Fraction) -> Fraction:
        return self.mass * buffer + self.reserve - self.cost

    def record(self, buffer: Fraction) -> dict:
        return {
            "normal": list(self.normal),
            "width": str(self.width(buffer)),
            "buffer_coefficient": str(self.mass),
            "deficit_reserve": str(self.reserve),
            "net_transform_cost": str(self.cost),
            "phases": self.phases,
        }


class ResourceLimit(Exception):
    pass


def canonical(slab: Slab) -> Slab | None:
    divisor = reduce(gcd, map(abs, slab.normal), 0)
    if not divisor:
        return None
    orientation = 1 if next(x for x in slab.normal if x) > 0 else -1
    return Slab(
        tuple(orientation * x // divisor for x in slab.normal),
        slab.mass / divisor,
        slab.reserve / divisor,
        slab.cost / divisor,
        slab.phases,
    )


def simplify(slabs, buffer: Fraction, *, prune_cube=True, max_facets=20_000):
    """Preserve every Boolean solution; dominance is ONLY at this buffer."""
    kept = {}
    for raw in slabs:
        if raw.width(buffer) <= 0:
            raise ValueError("empty open slab in a purported feasible transform")
        slab = canonical(raw)
        if slab is None:
            continue
        if prune_cube and slab.width(buffer) > sum(map(abs, slab.normal)):
            continue
        old = kept.get(slab.normal)
        if old is None or slab.width(buffer) < old.width(buffer):
            kept[slab.normal] = slab
        if len(kept) > max_facets:
            raise ResourceLimit(f"more than {max_facets} retained facets")
    return tuple(kept[key] for key in sorted(kept))


def contains(slabs, point, buffer: Fraction) -> bool:
    return all(
        abs(sum(c * z for c, z in zip(s.normal, point))) < s.width(buffer)
        for s in slabs
    )


def pivot_margin(slabs, pivot: int, buffer: Fraction) -> Fraction:
    return min(
        (s.width(buffer) - abs(s.normal[pivot]) for s in slabs),
        default=Fraction(1),
    )


def transform_section(
    slabs, pivot: int, buffer: Fraction, *, max_pairs=200_000,
    max_facets=20_000, prune_cube=True,
):
    """Section z_i=0 of ((K-e_i) intersect (K+e_i))+(-2,2)e_i.

    The retained b+|c_i| slabs are essential: the pair constraints alone
    describe a projection, which does not justify Boolean sign recovery.
    """
    if pivot_margin(slabs, pivot, buffer) <= 0:
        raise ValueError("pivot is not strictly feasible")
    oriented = []
    for s in slabs:
        p = s.normal[pivot]
        if p:
            oriented.append(replace(s, normal=tuple(
                (1 if p > 0 else -1) * x for x in s.normal
            )))
    pair_count = len(oriented) * (len(oriented) - 1) // 2
    if pair_count > max_pairs:
        raise ResourceLimit(f"{pair_count} pairs exceed limit {max_pairs}")

    def generated():
        # Original slabs expand before the zero-coordinate section is taken.
        for s in slabs:
            c = list(s.normal)
            p = abs(c[pivot])
            c[pivot] = 0
            yield replace(s, normal=tuple(c), cost=s.cost - p)
        for c, d in combinations(oriented, 2):
            p, q = c.normal[pivot], d.normal[pivot]
            yield Slab(
                tuple(q * x - p * y for x, y in zip(c.normal, d.normal)),
                q * c.mass + p * d.mass,
                q * c.reserve + p * d.reserve,
                q * c.cost + p * d.cost + 2 * p * q,
                c.phases | d.phases,
            )

    return simplify(
        generated(), buffer, prune_cube=prune_cube, max_facets=max_facets,
    ), pair_count


def construct_row(slabs, n: int, buffer: Fraction, *, order=None,
                  max_pairs=200_000, max_facets=20_000) -> dict:
    if order is not None and sorted(order) != list(range(n)):
        raise ValueError("order must be a permutation of all vertex indices")
    original = tuple(slabs)
    history, trace = [], []
    remaining = set(range(n))
    try:
        current = simplify(original, buffer, max_facets=max_facets)
        while remaining:
            margins = {i: pivot_margin(current, i, buffer) for i in remaining}
            feasible = [i for i in sorted(remaining) if margins[i] > 0]
            if order is not None:
                selected = order[len(history)]
                feasible = [selected] if margins[selected] > 0 else []
            if not feasible:
                return {
                    "status": "TRANSFORM_STOPPED", "trace": trace,
                    "scope": "this transform path only; extension feasibility unresolved",
                    "remaining": sorted(remaining),
                    "blockers": {
                        str(i): min(current, key=lambda s:
                            s.width(buffer) - abs(s.normal[i])).record(buffer)
                        for i in sorted(remaining)
                        if margins[i] <= 0
                    },
                }
            # Deterministic bounded heuristic, not a minimax ordering theorem.
            pivot = min(feasible, key=lambda i: (
                sum(s.normal[i] != 0 for s in current), -margins[i], i,
            ))
            child, pairs = transform_section(
                current, pivot, buffer, max_pairs=max_pairs,
                max_facets=max_facets,
            )
            history.append((pivot, current))
            trace.append({
                "pivot": pivot, "parent_facets": len(current),
                "child_facets": len(child), "pairs": pairs,
                "strict_margin": str(margins[pivot]),
                "positive_deficit_facets": sum(s.reserve > 0 for s in child),
                "mixed_phase_facets": sum(s.phases == 3 for s in child),
            })
            current = child
            remaining.remove(pivot)
    except ResourceLimit as exc:
        return {
            "status": "RESOURCE_LIMIT", "reason": str(exc), "trace": trace,
            "scope": "inconclusive; no extension infeasibility claim",
        }

    row = [0] * n
    for pivot, parents in reversed(history):
        for sign in (1, -1):
            row[pivot] = sign
            if contains(parents, row, buffer):
                break
        else:
            raise AssertionError("transform sign recovery failed")
    if not contains(original, row, buffer):
        raise AssertionError("recovered row violates an original slab")
    return {"status": "ROW_FOUND", "row": row, "trace": trace}


def load_matrix(path: Path) -> tuple[tuple[int, ...], ...]:
    text = path.read_text()
    try:
        obj = json.loads(text)
    except json.JSONDecodeError:
        obj = [list(map(int, line.replace(",", " ").split()))
               for line in text.splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    if isinstance(obj, dict):
        obj = next((obj[k] for k in ("A", "matrix", "M") if k in obj), None)
    if not isinstance(obj, list) or len(obj) < 2:
        raise ValueError("expected a signing of order at least two")
    n = len(obj)
    if any(not isinstance(row, list) or len(row) != n for row in obj):
        raise ValueError("matrix is not square")
    if any(type(obj[i][j]) is not int or obj[i][j] not in
           ((0,) if i == j else (-1, 1)) or obj[i][j] != obj[j][i]
           for i in range(n) for j in range(n)):
        raise ValueError("expected symmetric, zero-diagonal, integral sign matrix")
    return tuple(map(tuple, obj))


def matrix_hash(matrix) -> str:
    return sha256(json.dumps(matrix, separators=(",", ":")).encode()).hexdigest()


def enumerate_skeleton(matrix, *, max_states=1 << 20) -> dict:
    """One exact Gray-code pass; no second enumeration for row validation."""
    n = len(matrix)
    total = 1 << (n - 1)
    if total > max_states:
        raise ResourceLimit(f"{total} projective states exceed limit {max_states}")
    x = [1] * n
    fields = [sum(row) for row in matrix]
    q = sum(fields) // 2
    norm, stable = 0, []
    for step in range(total):
        if step:
            j = (step & -step).bit_length()  # coordinate zero remains +1
            old = x[j]
            q -= 2 * old * fields[j]
            x[j] = -old
            for i in range(n):
                fields[i] -= 2 * old * matrix[i][j]
        norm = max(norm, abs(q))
        local = [x[i] * fields[i] for i in range(n)]
        phase = (1 if min(local) >= 0 else 0) | (2 if max(local) <= 0 else 0)
        if phase:
            stable.append({"x": x.copy(), "energy": abs(q), "phases": phase})
    return {
        "schema": "complete_stable_skeleton_v1", "n": n, "phi": norm,
        "matrix_sha256": matrix_hash(matrix), "projective_states": total,
        "completeness": "exhaustive_projective_enumeration", "states": stable,
    }


def validate_skeleton(skeleton, matrix):
    """Check each supplied state, not the completeness claim of a cache."""
    n = len(matrix)
    if (not isinstance(skeleton, dict)
            or skeleton.get("schema") != "complete_stable_skeleton_v1"
            or skeleton.get("matrix_sha256") != matrix_hash(matrix)
            or skeleton.get("n") != n
            or skeleton.get("projective_states") != 1 << (n - 1)):
        raise ValueError("skeleton identity or enumeration metadata mismatch")
    if not isinstance(skeleton.get("states"), list):
        raise ValueError("skeleton has no state list")
    seen = set()
    for state in skeleton["states"]:
        if (not isinstance(state, dict) or not isinstance(state.get("x"), list)
                or type(state.get("energy")) is not int
                or type(state.get("phases")) is not int):
            raise ValueError("invalid stable-state record")
        x = tuple(state["x"])
        if len(x) != n or x[0] != 1 or any(type(v) is not int or v not in (-1, 1) for v in x):
            raise ValueError("invalid projective Boolean state")
        if x in seen:
            raise ValueError("duplicate projective state")
        seen.add(x)
        fields = [sum(a * v for a, v in zip(row, x)) for row in matrix]
        local = [a * v for a, v in zip(fields, x)]
        phase = (1 if min(local) >= 0 else 0) | (2 if max(local) <= 0 else 0)
        if not phase or state["phases"] != phase or state["energy"] != abs(sum(local) // 2):
            raise ValueError("incorrect stable-state energy or orientation")
    if (not seen or type(skeleton.get("phi")) is not int
            or skeleton["phi"] != max(s["energy"] for s in skeleton["states"])):
        raise ValueError("skeleton maximum does not equal its declared norm")


def target_for_neutral_extension(phi: int, n: int, excess: int) -> int:
    # Exact floor of phi*((n+1)/n)**(3/2), with no floating-point boundary.
    target = isqrt(phi * phi * (n + 1) ** 3 // n ** 3) + excess
    return target - (target - n * (n + 1) // 2) % 2


def write_new_json(path: Path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--enumerate-stable", action="store_true")
    mode.add_argument("--skeleton", type=Path)
    parser.add_argument("--save-skeleton", type=Path)
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--target-phi", type=int)
    target.add_argument("--excess", type=int, default=0,
                        help="nonnegative integer allowance above neutral target")
    parser.add_argument("--order", help="comma-separated vertex permutation")
    parser.add_argument("--max-states", type=int, default=1 << 20)
    parser.add_argument("--max-pairs", type=int, default=200_000)
    parser.add_argument("--max-facets", type=int, default=20_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.excess < 0 or min(args.max_states, args.max_pairs, args.max_facets) < 1:
        parser.error("excess must be nonnegative and resource limits positive")
    for path in (args.output, args.save_skeleton):
        if path is not None and path.exists():
            parser.error(f"refusing to overwrite {path}")
    if args.output is not None and args.output == args.save_skeleton:
        parser.error("output and skeleton paths must differ")
    matrix = load_matrix(args.matrix)
    skeleton = (enumerate_skeleton(matrix, max_states=args.max_states)
                if args.enumerate_stable else json.loads(args.skeleton.read_text()))
    validate_skeleton(skeleton, matrix)
    if args.save_skeleton:
        write_new_json(args.save_skeleton, skeleton)
    phi, n = skeleton["phi"], len(matrix)
    threshold = (args.target_phi if args.target_phi is not None else
                 target_for_neutral_extension(phi, n, args.excess))
    if threshold < phi:
        parser.error("target is below source norm; restriction already excludes it")
    # All scores are integral. The half-unit opens the body without relaxing
    # the requested integer target.
    buffer = Fraction(2 * (threshold - phi) + 1, 2)
    slabs = tuple(Slab(tuple(s["x"]), reserve=Fraction(phi - s["energy"]),
                       phases=s["phases"]) for s in skeleton["states"])
    result = construct_row(
        slabs, n, buffer,
        order=None if args.order is None else list(map(int, args.order.split(","))),
        max_pairs=args.max_pairs, max_facets=args.max_facets,
    )
    result.update({
        "schema": "deficit_transform_certificate_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "matrix_sha256": matrix_hash(matrix),
        "program_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "skeleton_sha256": sha256(json.dumps(skeleton, sort_keys=True,
            separators=(",", ":")).encode()).hexdigest(),
        "n": n, "source_phi": phi, "target_phi": threshold,
        "open_buffer": str(buffer), "stable_classes": len(slabs),
        "pivot_rule": "fewest_nonzero_pivots_then_largest_margin_then_index",
        "prescribed_order": args.order,
        "limits": {"states": args.max_states, "pairs": args.max_pairs,
                   "facets": args.max_facets},
        "completeness": ("enumerated_in_this_run" if args.enumerate_stable else
                         "cache_claim_not_independently_reverified"),
        "source_phi_scope": ("exact_in_this_run" if args.enumerate_stable else
                             "cached_claim_only"),
        "source_global_minimality": "not_established_by_this_program",
        "convergence_proved": False,
    })
    if result["status"] == "ROW_FOUND":
        scores = [s["energy"] + abs(sum(a * x for a, x in zip(result["row"], s["x"])))
                  for s in skeleton["states"]]
        maximum = max(scores)
        if maximum > threshold:
            raise AssertionError("original stable-score verification failed")
        result["stable_score"] = maximum
        result["maximizing_stable_state"] = skeleton["states"][scores.index(maximum)]
        result["certificate"] = ("exact_finite_extension" if args.enumerate_stable
                                 else "conditional_on_cached_skeleton_completeness")
        if args.enumerate_stable:
            result["extension_phi"] = maximum
    if args.output:
        write_new_json(args.output, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "ROW_FOUND" else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, ResourceLimit, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
