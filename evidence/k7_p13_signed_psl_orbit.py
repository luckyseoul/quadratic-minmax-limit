#!/usr/bin/env python3
"""Enumerate the signed-PSL orbit of the exact p=13, k=7 witness.

The signed generators are the conference automorphisms from Prop 15.588.
Vectors are represented by one Python integer whose set bits are the negative
coordinates.  Byte lookup tables make each signed permutation inexpensive.

The optional packed output has three uint64 words per orbit vector and is the
input to ``k7_p13_orbit_quartic_xpu.py``.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15588 import aut_generators, signed_matrix  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def negative_key(y: np.ndarray) -> int:
    """Pack a sign vector into an integer, with one bit per negative entry."""
    key = 0
    for index in np.flatnonzero(np.asarray(y) < 0):
        key |= 1 << int(index)
    return key


def sign_vector(key: int, n: int) -> np.ndarray:
    """Unpack a negative-bit integer into a +-1 vector."""
    return np.fromiter(
        (-1 if (key >> index) & 1 else 1 for index in range(n)),
        dtype=np.int8,
        count=n,
    )


def transform_key_direct(key: int, perm: np.ndarray, signs: np.ndarray) -> int:
    """Reference signed-permutation action: out[perm[i]]=signs[i]*in[i]."""
    out = 0
    for source, target in enumerate(perm):
        negative = ((key >> source) & 1) ^ int(signs[source] < 0)
        out |= negative << int(target)
    return out


def transform_tables(
    generators: list[tuple[np.ndarray, np.ndarray]], n: int
) -> tuple[list[int], list[list[list[int]]]]:
    """Build byte lookup tables for each signed permutation."""
    blocks = (n + 7) // 8
    sign_masks: list[int] = []
    tables: list[list[list[int]]] = []
    for perm, signs in generators:
        sign_mask = 0
        for source, target in enumerate(perm):
            if signs[source] < 0:
                sign_mask |= 1 << int(target)
        generator_tables: list[list[int]] = []
        for block in range(blocks):
            rows = [0] * 256
            for value in range(256):
                contribution = 0
                for offset in range(8):
                    source = 8 * block + offset
                    if source < n and (value >> offset) & 1:
                        contribution |= 1 << int(perm[source])
                rows[value] = contribution
            generator_tables.append(rows)
        sign_masks.append(sign_mask)
        tables.append(generator_tables)
    return sign_masks, tables


def transform_key(key: int, sign_mask: int, tables: list[list[int]]) -> int:
    out = sign_mask
    for block, rows in enumerate(tables):
        out ^= rows[(key >> (8 * block)) & 255]
    return out


def enumerate_orbit(
    seed: int,
    sign_masks: list[int],
    tables: list[list[list[int]]],
    progress: int,
    limit: int | None,
) -> set[int]:
    """Breadth-first orbit enumeration under the supplied generators."""
    seen = {seed}
    frontier = deque([seed])
    processed = 0
    started = time.monotonic()
    while frontier:
        key = frontier.popleft()
        processed += 1
        for sign_mask, generator_tables in zip(sign_masks, tables):
            image = transform_key(key, sign_mask, generator_tables)
            if image not in seen:
                seen.add(image)
                frontier.append(image)
                if limit is not None and len(seen) >= limit:
                    return seen
        if progress and processed % progress == 0:
            elapsed = time.monotonic() - started
            print(
                f"processed={processed:,} orbit={len(seen):,} "
                f"frontier={len(frontier):,} elapsed={elapsed:.1f}s",
                flush=True,
            )
    return seen


def save_packed(path: Path, orbit: set[int], n: int) -> None:
    """Save little-endian uint64 words, one row per orbit vector."""
    words = (n + 63) // 64
    packed = np.empty((len(orbit), words), dtype=np.uint64)
    mask = (1 << 64) - 1
    for row, key in enumerate(orbit):
        for word in range(words):
            packed[row, word] = (key >> (64 * word)) & mask
    np.save(path, packed)


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--witness",
        type=Path,
        default=ROOT / "evidence" / "k7_p13_cpsat_witness.json",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--packed", type=Path)
    parser.add_argument("--progress", type=int, default=100_000)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    witness = json.loads(args.witness.read_text())
    p = int(witness["p"])
    q = p * p
    n = q + 1
    y = np.ones(n, dtype=np.int8)
    witness_record = witness.get("best_Zpsi_seen") or witness.get(
        "outside_orbit_Zpsi"
    )
    if witness_record is None:
        raise RuntimeError("witness has no supported negative-index record")
    finite_negative = np.asarray(
        witness_record["negative_indices"], dtype=np.int64
    )
    y[1 + finite_negative] = -1
    seed = negative_key(y)

    conference = paley_conference_prime_power(p).astype(np.int64)
    if not np.array_equal(conference @ y.astype(np.int64), p * y):
        raise RuntimeError("stored witness is not in Max+")
    generators = aut_generators(p)
    for perm, signs in generators:
        matrix = signed_matrix(perm, signs, n).astype(np.int64)
        if not np.array_equal(matrix @ conference @ matrix.T, conference):
            raise RuntimeError("signed generator does not preserve C")
        direct = transform_key_direct(seed, perm, signs)
        image = matrix @ y.astype(np.int64)
        if direct != negative_key(image):
            raise RuntimeError("signed-bit action convention mismatch")
        if not np.array_equal(conference @ image, p * image):
            raise RuntimeError("signed generator does not preserve Max+")

    sign_masks, tables = transform_tables(generators, n)
    for index, (perm, signs) in enumerate(generators):
        if transform_key(seed, sign_masks[index], tables[index]) != (
            transform_key_direct(seed, perm, signs)
        ):
            raise RuntimeError("byte lookup action mismatch")

    started = time.monotonic()
    orbit = enumerate_orbit(
        seed, sign_masks, tables, args.progress, args.limit
    )
    elapsed = time.monotonic() - started
    projective_psl_order = q * (q * q - 1) // 2
    signed_lift_order = 2 * projective_psl_order
    complete = args.limit is None or len(orbit) < args.limit
    if complete and signed_lift_order % len(orbit):
        raise RuntimeError("orbit size does not divide the signed PSL lift")

    epsilon_plus = sum(1 for key in orbit if not (key & 1))
    report = {
        "p": p,
        "n": n,
        "algorithm": "signed-PSL BFS with packed negative-bit keys",
        "generator_count": len(generators),
        "generator_conference_audit": True,
        "orbit_size": len(orbit),
        "projective_psl_order": projective_psl_order,
        "signed_psl_lift_order": signed_lift_order,
        "signed_stabilizer_order": (
            signed_lift_order // len(orbit) if complete else None
        ),
        "epsilon_plus_count": epsilon_plus,
        "epsilon_minus_count": len(orbit) - epsilon_plus,
        "complete_orbit": complete,
        "elapsed_seconds": elapsed,
        "source_witness": str(args.witness),
    }
    if args.packed:
        args.packed.parent.mkdir(parents=True, exist_ok=True)
        save_packed(args.packed, orbit, n)
        report["packed_orbit"] = str(args.packed)
        report["packed_sha256"] = __import__("hashlib").sha256(
            args.packed.read_bytes()
        ).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    main()
