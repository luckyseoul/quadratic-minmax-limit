#!/usr/bin/env python3
"""Classify p=11 quartic profiles for square-circle channel moments.

The ordinary profile counter only needs the value histogram of a quartic
polynomial on F_11.  The channel trace additionally uses

    U4(a) = sum_c (sum_s eta(s-c) a_s)^4,

which remembers the cyclic placement of the polynomial values.  U4 is
unchanged by affine permutations of the input coordinate.  The 11^4
zero-constant quartics therefore collapse to 1,007 domain-affine profile
types.  (The zero-constant family is not itself invariant under input
translation, which is why the count is larger than the 147 orbits obtained
after also quotienting output translation.)  Allowing an affine
transformation of the output collapses those 1,007 types to 20
dynamic-programming tables; the output affine descriptor
reconstructs the character phase exactly.

This module builds both reductions without heuristic hashing.  It is used
by the rich tuple and marked-profile generators for the channel-resolved
p=11 R1 calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r1_p11_profile_dual_orbits import P


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_sequence(coefficients: tuple[int, int, int, int]) -> tuple[int, ...]:
    """Values of c1*x+...+c4*x^4 on F_11, in coordinate order."""
    return tuple(
        sum(coefficient * pow(value, degree, P) for degree, coefficient in enumerate(coefficients, 1))
        % P
        for value in range(P)
    )


def domain_images(sequence: tuple[int, ...]):
    """Yield all input-affine permutations x -> alpha*x+beta."""
    for alpha in range(1, P):
        for beta in range(P):
            yield tuple(sequence[(alpha * value + beta) % P] for value in range(P))


def domain_canonical(sequence: tuple[int, ...]) -> tuple[int, ...]:
    return min(domain_images(sequence))


def output_affine_canonical(
    sequence: tuple[int, ...],
) -> tuple[tuple[int, ...], int, int]:
    """Return canonical sequence and ``canonical=alpha*domain(seq)+beta``.

    The selected domain permutation need not be recorded: profile sums,
    energies, and U4 are invariant under it.  The output parameters are
    required to transport the additive-character phase.
    """
    best: tuple[int, ...] | None = None
    best_alpha = 0
    best_beta = 0
    for domain_image in domain_images(sequence):
        for alpha in range(1, P):
            scaled = tuple(alpha * value % P for value in domain_image)
            for beta in range(P):
                candidate = tuple((value + beta) % P for value in scaled)
                if best is None or candidate < best:
                    best = candidate
                    best_alpha = alpha
                    best_beta = beta
    if best is None:
        raise ArithmeticError("empty affine orbit")
    return best, best_alpha, best_beta


def rich_profile_data() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return packed lookup, 1,007 representatives, 20 tables, descriptors.

    ``lookup[c1+11*c2+11^2*c3+11^3*c4]`` is the rich domain-affine ID.
    A descriptor row is ``(table_id, alpha, beta)`` and means that, after
    an irrelevant domain permutation,

        table_sequence = alpha * rich_sequence + beta.
    """
    canonical_by_packed: list[tuple[int, ...]] = []
    for c4 in range(P):
        for c3 in range(P):
            for c2 in range(P):
                for c1 in range(P):
                    canonical_by_packed.append(
                        domain_canonical(polynomial_sequence((c1, c2, c3, c4)))
                    )
    representatives = sorted(set(canonical_by_packed))
    if len(representatives) != 1007:
        raise ArithmeticError(
            f"quartic domain-affine orbit count changed to {len(representatives)}"
        )
    rich_ids = {sequence: index for index, sequence in enumerate(representatives)}
    lookup = np.asarray(
        [rich_ids[sequence] for sequence in canonical_by_packed], dtype=np.uint16
    )

    table_representatives: list[tuple[int, ...]] = []
    table_ids: dict[tuple[int, ...], int] = {}
    descriptors: list[tuple[int, int, int]] = []
    for sequence in representatives:
        canonical, alpha, beta = output_affine_canonical(sequence)
        table_id = table_ids.get(canonical)
        if table_id is None:
            table_id = len(table_representatives)
            table_ids[canonical] = table_id
            table_representatives.append(canonical)
        descriptors.append((table_id, alpha, beta))
    if len(table_representatives) != 20:
        raise ArithmeticError(
            "quartic domain/output-affine table count changed to "
            f"{len(table_representatives)}"
        )
    return (
        lookup,
        np.asarray(representatives, dtype=np.uint8),
        np.asarray(table_representatives, dtype=np.uint8),
        np.asarray(descriptors, dtype=np.uint8),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    lookup, rich, tables, descriptors = rich_profile_data()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        lookup=lookup,
        rich_sequences=rich,
        canonical_table_sequences=tables,
        affine_descriptors=descriptors,
    )
    report = {
        "experiment": "r1_p11_channel_profile_types",
        "status": "complete_exact_affine_profile_classification",
        "p": P,
        "quartic_coefficient_tuples": P**4,
        "domain_affine_profile_types": int(len(rich)),
        "domain_output_affine_table_types": int(len(tables)),
        "lookup_dtype": str(lookup.dtype),
        "output": str(args.output),
        "output_sha256": sha256(args.output),
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
