#!/usr/bin/env python3
r"""Exact shardable search for circle trades of the triangular completion.

This is an unnumbered diagnostic for the live residual-(ii) obstruction.  It
does not enumerate residual graphs.  For the Kiss--Somlai Boolean shadow
``x`` it scans finite square Miquelian circles ``S`` and asks whether

    C_ij x_i x_j = +1  for every i != j in S.

In that case ``-x 1_S`` is a sparse ``+p`` eigenvector, so flipping ``x`` on
``S`` gives another Boolean completion of the *same* signed three-spike
datum.  If ``S`` is disjoint from the three spikes, its intersection with
the distinguished spike circle is exactly the new circle mismatch ``mu``.

The scan is split by the integer encoding of the circle centre.  Every hit
is rechecked on all pairs, so the output is an exact finite certificate and
does not depend on a floating-point eigensolver.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from e1_gmin_m4_prop15598 import field_ctx, legendre
from residual_kiss_somlai_three_spike import (
    _add,
    _encode,
    _scalar,
    find_square_triangle_linear_map,
    transform_function,
    triangular_augmented_function,
)


def _norm_table(p: int, ia: int, ib: int, da: np.ndarray, db: np.ndarray) -> np.ndarray:
    """Norm of ``da+db*w`` when ``w^2=ia*w+ib``."""
    return (da * da + ia * da * db - ib * db * db) % p


def _triangular_data(p: int) -> tuple[np.ndarray, set[int], set[int]]:
    image_a, image_b = find_square_triangle_linear_map(p)
    f = transform_function(
        p, triangular_augmented_function(p), image_a, image_b
    )
    x = np.where(f > 0, 1, -1).astype(np.int8)
    spikes = set(int(u) for u in np.flatnonzero(f == 2))
    if len(spikes) != 2:
        raise ArithmeticError("the triangular finite spike set changed")
    base = _scalar(p, p - 2, image_b)
    gamma = {_add(p, base, _scalar(p, t, image_a)) for t in range(p)}
    if not spikes <= gamma:
        raise ArithmeticError("the finite spikes left their distinguished line")
    return x, spikes, gamma


def _all_pair_recheck(
    p: int,
    ia: int,
    ib: int,
    x: np.ndarray,
    points: np.ndarray,
) -> bool:
    for pos, u in enumerate(points):
        ua, ub = int(u % p), int(u // p)
        for v in points[pos + 1 :]:
            va, vb = int(v % p), int(v // p)
            norm = ((ua - va) ** 2 + ia * (ua - va) * (ub - vb) - ib * (ub - vb) ** 2) % p
            if int(x[u]) * int(x[v]) * legendre(norm, p) != 1:
                return False
    return True


def scan(p: int, centre_start: int, centre_stop: int) -> dict[str, object]:
    q, _mul, _add_field, _chi, _frob, _norm, ia, ib = field_ctx(p)
    if not 0 <= centre_start <= centre_stop <= q:
        raise ValueError("centre shard must lie in 0..p^2")
    x, spikes, gamma = _triangular_data(p)
    radius_type = -legendre(-1, p)
    radii = tuple(r for r in range(1, p) if legendre(r, p) == radius_type)
    leg = np.asarray([legendre(v, p) for v in range(p)], dtype=np.int8)
    coords_a = np.arange(q, dtype=np.int64) % p
    coords_b = np.arange(q, dtype=np.int64) // p
    hits: list[dict[str, object]] = []
    checked = 0
    for centre in range(centre_start, centre_stop):
        ca, cb = centre % p, centre // p
        da = (coords_a - ca) % p
        db = (coords_b - cb) % p
        norms = _norm_table(p, ia, ib, da, db)
        for radius in radii:
            points = np.flatnonzero(norms == radius)
            if len(points) != p + 1:
                raise ArithmeticError("a finite norm circle has the wrong size")
            checked += 1
            if any(int(point) in spikes for point in points):
                continue
            base = int(points[0])
            ba, bb = base % p, base // p
            nda = (coords_a[points] - ba) % p
            ndb = (coords_b[points] - bb) % p
            edge_norms = _norm_table(p, ia, ib, nda, ndb)
            relations = int(x[base]) * x[points] * leg[edge_norms]
            # Ignore the zero diagonal at the chosen base point.
            if not np.all(relations[points != base] == 1):
                continue
            if not _all_pair_recheck(p, ia, ib, x, points):
                raise ArithmeticError("base-row clique test passed but pair check failed")
            intersection = sorted(int(point) for point in points if int(point) in gamma)
            hits.append(
                {
                    "centre": centre,
                    "centre_coordinates": [ca, cb],
                    "radius": radius,
                    "mu": len(intersection),
                    "spike_circle_intersection": intersection,
                }
            )
    histogram = Counter(int(row["mu"]) for row in hits)
    return {
        "p": p,
        "centre_start": centre_start,
        "centre_stop": centre_stop,
        "circle_count_checked": checked,
        "finite_spikes": sorted(spikes),
        "switchable_circle_count": len(hits),
        "mismatch_histogram": {str(mu): count for mu, count in sorted(histogram.items())},
        "positive_mismatch_exists": any(int(row["mu"]) > 0 for row in hits),
        "hits": hits,
        "all_hits_pair_rechecked": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--centre-start", type=int, default=0)
    parser.add_argument("--centre-stop", type=int)
    parser.add_argument("--output")
    args = parser.parse_args()
    stop = args.p * args.p if args.centre_stop is None else args.centre_stop
    result = scan(args.p, args.centre_start, stop)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
