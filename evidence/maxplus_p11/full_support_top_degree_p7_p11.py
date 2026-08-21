#!/usr/bin/env python3
"""Top-profile-degree split of the exact p=7 and p=11 full-support strata.

For p=3 (mod 4), let m=(p+1)/2 and d=m-2.  The degree-d
coefficients of the m active line profiles lie in a one-dimensional
homogeneous coefficient kernel.  This script conditions the normalized
quartic imbalance B=Z_psi/(2p) on its kernel scalar.

The p=11 input must be the *deduplicated* complete Max+ array, not the raw
k6 generator output (whose translation expansion has repeated rows).
"""
from __future__ import annotations

import json
import math
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from e1_gmin_m4_prop15588 import directions, field_ctx, maxplus  # noqa: E402


def primitive_root(q: int, mul) -> int:
    for generator in range(2, q):
        x = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(x)
            x = mul(x, generator)
        if len(seen) == q - 1:
            return generator
    raise RuntimeError("no primitive root")


def quartic_direction_weights(p: int) -> np.ndarray:
    q, mul, chi, trace = field_ctx(p)
    generator = primitive_root(q, mul)
    logs = {}
    x = 1
    for exponent in range(q - 1):
        logs[x] = exponent
        x = mul(x, generator)

    raw = []
    seen = set()
    for g in range(1, q):
        if g in seen:
            continue
        seen.update(mul(t, g) for t in range(1, p))
        annihilator = next(c for c in range(1, q) if trace(mul(c, g)) == 0)
        t_of = np.array(
            [trace(mul(annihilator, value)) for value in range(q)],
            dtype=np.int64,
        )
        if chi(g) == 1:
            raw.append((t_of, 1 if logs[g] % 4 == 0 else -1))

    weights = []
    for t_of, _form in directions(p)[0]:
        matches = [weight for candidate, weight in raw if np.array_equal(candidate, t_of)]
        if len(matches) != 1:
            raise RuntimeError("direction matching failed")
        weights.append(matches[0])
    return np.asarray(weights, dtype=np.int64)


def homogeneous_matrix(forms, degree: int, p: int) -> np.ndarray:
    return np.asarray(
        [
            [
                math.comb(degree, r)
                * pow(a, degree - r, p)
                * pow(b, r, p)
                % p
                for a, b in forms
            ]
            for r in range(degree + 1)
        ],
        dtype=np.int64,
    )


def free_columns(matrix: np.ndarray, p: int) -> list[int]:
    a = matrix.copy() % p
    row = 0
    pivots = []
    for column in range(a.shape[1]):
        pivot = next((r for r in range(row, a.shape[0]) if a[r, column]), None)
        if pivot is None:
            continue
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = a[row] * pow(int(a[row, column]), p - 2, p) % p
        for r in range(a.shape[0]):
            if r != row and a[r, column]:
                a[r] = (a[r] - a[r, column] * a[row]) % p
        pivots.append(column)
        row += 1
        if row == a.shape[0]:
            break
    return [column for column in range(a.shape[1]) if column not in pivots]


def leading_coefficients(rho: np.ndarray, degree: int, p: int) -> np.ndarray:
    delta = np.zeros(rho.shape[:2], dtype=np.int64)
    for s in range(degree + 1):
        delta += (-1) ** (degree - s) * math.comb(degree, s) * rho[:, :, s]
    return delta * pow(math.factorial(degree), p - 2, p) % p


def projective_keys(coords: np.ndarray, p: int) -> np.ndarray:
    """Encode P^1(F_p) as 0..p, reserving p+1 for the zero vector."""
    keys = np.full(len(coords), p + 1, dtype=np.int16)
    finite = coords[:, 0] != 0
    for first in range(1, p):
        mask = finite & (coords[:, 0] == first)
        keys[mask] = coords[mask, 1] * pow(first, p - 2, p) % p
    infinity = (coords[:, 0] == 0) & (coords[:, 1] != 0)
    keys[infinity] = p
    return keys


def projective_label(key: int, p: int) -> str:
    if key < p:
        return f"1,{key}"
    if key == p:
        return "0,1"
    return "zero"


def add_values(record: dict, values: np.ndarray) -> None:
    record["count"] += len(values)
    record["sum_B"] += int(values.sum())
    record["sum_B2"] += int(values @ values)
    distinct, counts = np.unique(values, return_counts=True)
    record["histogram"].update(
        {int(value): int(count) for value, count in zip(distinct, counts)}
    )


def render_record(record: dict) -> dict:
    count = record["count"]
    return {
        "count": count,
        "E_B": str(Fraction(record["sum_B"], count)),
        "E_B2": str(Fraction(record["sum_B2"], count)),
        "B_histogram": {
            str(value): count for value, count in sorted(record["histogram"].items())
        },
    }


def analyze(p: int, finite_rows: np.ndarray, chunk: int = 250_000) -> dict:
    q = p * p
    square, _ = directions(p)
    m = len(square)
    degree = m - 2
    forms = [form for _t_of, form in square]
    top_matrix = homogeneous_matrix(forms, degree, p)
    top_free = free_columns(top_matrix, p)
    drop_matrix = homogeneous_matrix(forms, degree - 1, p)
    drop_free = free_columns(drop_matrix, p)
    if len(top_free) != 1 or len(drop_free) != 2:
        raise RuntimeError("unexpected homogeneous-kernel dimensions")

    incidence = np.zeros((q, m * p), dtype=np.float32)
    for j, (t_of, _form) in enumerate(square):
        incidence[np.arange(q), j * p + t_of] = 1
    weights = quartic_direction_weights(p)
    top = defaultdict(
        lambda: {"count": 0, "sum_B": 0, "sum_B2": 0, "histogram": Counter()}
    )
    drop = defaultdict(
        lambda: {"count": 0, "sum_B": 0, "sum_B2": 0, "histogram": Counter()}
    )
    n_full_support = 0
    n_drop_twice = 0

    for lo in range(0, len(finite_rows), chunk):
        hi = min(len(finite_rows), lo + chunk)
        rows = np.asarray(finite_rows[lo:hi], dtype=np.float32)
        profiles = np.rint(rows @ incidence).astype(np.int16).reshape(-1, m, p)
        full = np.any(profiles != 1, axis=2).sum(axis=1) == m
        profiles = profiles[full]
        if not len(profiles):
            continue
        n_full_support += len(profiles)
        rho = ((profiles + p - 2) // 2) % p
        coefficient = leading_coefficients(rho, degree, p)
        if np.any(top_matrix @ coefficient.T % p):
            raise RuntimeError("top coefficient kernel failed")
        scalar = coefficient[:, top_free[0]]

        h = ((profiles - 1) // 2).astype(np.int64)
        energy = np.sum(h * h, axis=2)
        if np.any(energy % (2 * p)):
            raise RuntimeError("2p energy divisibility failed")
        B = (energy // (2 * p)) @ weights
        for value in range(p):
            add_values(top[value], B[scalar == value])

        zero = scalar == 0
        if np.any(zero):
            lower = leading_coefficients(rho[zero], degree - 1, p)
            if np.any(drop_matrix @ lower.T % p):
                raise RuntimeError("degree-drop coefficient kernel failed")
            coords = lower[:, drop_free]
            n_drop_twice += int(np.all(coords == 0, axis=1).sum())
            keys = projective_keys(coords, p)
            for encoded in np.unique(keys):
                key = projective_label(int(encoded), p)
                add_values(drop[key], B[zero][keys == encoded])

    total_B2 = sum(record["sum_B2"] for record in top.values())
    target = Fraction(3 * (p * p - 1), 64)
    return {
        "p": p,
        "full_support_count": n_full_support,
        "top_profile_degree": degree,
        "quartic_direction_weights": weights.tolist(),
        "normalized_QVAR_threshold": str(target),
        "top_kernel_scalar_classes": {
            str(value): render_record(top[value])
            for value in range(p)
            if top[value]["count"]
        },
        "top_zero_count": top[0]["count"],
        "degree_drops_twice_count": n_drop_twice,
        "degree_drop_projective_classes": {
            key: render_record(record) for key, record in sorted(drop.items())
        },
        "aggregate_E_B2": str(Fraction(total_B2, n_full_support)),
    }


def main() -> dict:
    p7 = maxplus(7)
    p7_finite = p7[p7[:, 0] == 1, 1:]
    root = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    p11_all = np.load(root / "maxplus_p11_eps1.npy", mmap_mode="r")
    p11_k6_start = 2_772 + 24_200 + 58_080 + 1_306_800
    p11_finite = p11_all[p11_k6_start:, 1:]
    if p11_finite.shape != (36_065_260, 121):
        raise RuntimeError(f"unexpected deduplicated p=11 k6 shape {p11_finite.shape}")

    report = {"p7": analyze(7, p7_finite), "p11": analyze(11, p11_finite)}
    seven, eleven = report["p7"], report["p11"]
    seven_nonzero = list(seven["top_kernel_scalar_classes"].values())
    eleven_nonzero = [
        record
        for scalar, record in eleven["top_kernel_scalar_classes"].items()
        if scalar != "0"
    ]
    drop_summary = Counter(
        (record["count"], record["E_B2"])
        for record in eleven["degree_drop_projective_classes"].values()
    )
    if not (
        seven["full_support_count"] == 4_410
        and seven["top_zero_count"] == 0
        and len(seven_nonzero) == 6
        and {record["count"] for record in seven_nonzero} == {735}
        and {record["E_B2"] for record in seven_nonzero} == {"44/15"}
        and eleven["full_support_count"] == 36_065_260
        and eleven["top_zero_count"] == 2_090_880
        and eleven["top_kernel_scalar_classes"]["0"]["E_B2"] == "137/36"
        and len(eleven_nonzero) == 10
        and {record["count"] for record in eleven_nonzero} == {3_397_438}
        and {record["E_B2"] for record in eleven_nonzero} == {"111483/14039"}
        and eleven["degree_drops_twice_count"] == 0
        and drop_summary == Counter({(123_420, "151/51"): 6, (225_060, "397/93"): 6})
        and eleven["aggregate_E_B2"] == "114771/14903"
    ):
        raise RuntimeError("full-support top-degree split audit failed")

    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
