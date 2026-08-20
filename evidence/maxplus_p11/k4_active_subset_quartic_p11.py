#!/usr/bin/env python3
"""Exact p=11 k=4 quartic moment split by active direction-subset.

The input is the persisted eps=+1 finite-coordinate k=4 array.  Translation
does not change line-profile energies, so selecting the zero-linear-level
representative in each orbit reduces 58,080 rows to 480 pure parabolas.
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from e1_gmin_m4_prop15588 import directions, field_ctx, profiles_of  # noqa: E402

P = 11
Q = P * P


def primitive_root(q: int, mul) -> int:
    for g in range(2, q):
        x = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(x)
            x = mul(x, g)
        if len(seen) == q - 1:
            return g
    raise RuntimeError("no primitive root")


def quartic_direction_weights() -> np.ndarray:
    q, mul, chi, trace = field_ctx(P)
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
        seen.update(mul(t, g) for t in range(1, P))
        annihilator = next(c for c in range(1, q) if trace(mul(c, g)) == 0)
        t_of = np.array(
            [trace(mul(annihilator, value)) for value in range(q)],
            dtype=np.int64,
        )
        if chi(g) == 1:
            raw.append((t_of, 1 if logs[g] % 4 == 0 else -1))

    weights = []
    for t_of, _form in directions(P)[0]:
        matches = [weight for candidate, weight in raw if np.array_equal(candidate, t_of)]
        if len(matches) != 1:
            raise RuntimeError("direction matching failed")
        weights.append(matches[0])
    return np.asarray(weights, dtype=np.int64)


def main() -> dict:
    work = Path(os.environ.get("E1WORK_P11", "/mnt/storage/e1work/maxplus_p11"))
    finite = np.load(work / "k4_p11_full.npy", mmap_mode="r")
    if finite.shape != (58_080, Q):
        raise RuntimeError(f"unexpected k4 shape {finite.shape}")
    rows = np.column_stack((np.ones(len(finite), dtype=np.int8), finite))
    profiles = profiles_of(P, rows)
    active = np.any(profiles != 1, axis=2)

    # For rho(s)=a0+a1*s+a2*s^2, recover a1 from s=0,1,2.
    rho = ((profiles + P - 2) // 2) % P
    a0 = rho[:, :, 0]
    r1 = (rho[:, :, 1] - a0) % P
    r2 = (rho[:, :, 2] - a0) % P
    a2 = ((r2 - 2 * r1) * pow(2, P - 2, P)) % P
    a1 = (r1 - a2) % P
    pure = (active.sum(axis=1) == 4) & np.all((a1 == 0) | (~active), axis=1)
    if int(pure.sum()) != 480:
        raise RuntimeError(f"expected 480 pure reps, got {int(pure.sum())}")

    h = ((profiles[pure] - 1) // 2).astype(np.int64)
    energy = (h * h).sum(axis=2)
    if np.any(energy % (2 * P)):
        raise RuntimeError("2p divisibility failed")
    b = energy // (2 * P)
    weights = quartic_direction_weights()
    B = b @ weights

    selected_rows = np.where(pure)[0]
    grouped: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for position, row in enumerate(selected_rows):
        subset = tuple(int(j) for j in np.where(active[row])[0])
        grouped[subset].append(int(B[position]))

    subsets = {}
    n_balanced = n_unbalanced = 0
    for subset, values in sorted(grouped.items()):
        sign_sum = int(weights[list(subset)].sum())
        kind = "balanced" if sign_sum == 0 else "unbalanced"
        n_balanced += kind == "balanced"
        n_unbalanced += kind == "unbalanced"
        moment = Fraction(sum(value * value for value in values), len(values))
        subsets[",".join(map(str, subset))] = {
            "kind": kind,
            "n_pure_reps": len(values),
            "B_histogram": {
                str(value): count for value, count in sorted(Counter(values).items())
            },
            "E_B2": str(moment),
        }

    aggregate = Fraction(sum(int(value) ** 2 for value in B), len(B))
    report = {
        "p": P,
        "quartic_direction_weights": weights.tolist(),
        "n_full_vectors": len(rows),
        "translation_orbit_size": Q,
        "n_pure_reps": int(pure.sum()),
        "n_balanced_subsets": n_balanced,
        "n_unbalanced_subsets": n_unbalanced,
        "normalized_QVAR_threshold": "45/8",
        "aggregate_E_B2": str(aggregate),
        "aggregate_E_Z2": str(aggregate * (2 * P) ** 2),
        "subsets": subsets,
    }
    if (
        n_balanced != 9
        or n_unbalanced != 6
        or aggregate != Fraction(39, 2)
        or {rec["E_B2"] for rec in subsets.values()} != {"5", "63"}
    ):
        raise RuntimeError("p=11 subset split audit failed")

    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
