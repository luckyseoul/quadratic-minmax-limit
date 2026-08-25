#!/usr/bin/env python3
"""Export complete exact shells of the Paley dual lattice.

Let ``B`` be a saturated integer basis of ``ker(C-pI)`` and
``Q=2p(B^T B)^-1``.  PARI's ``qfminim`` enumerates integral dual
coordinates ``k``.  The corresponding ambient vector is

    u = B Q k / (2p).

The archive stores one member of every ``+/-`` pair, both as ``k`` and as
the integral numerator ``B Q k``.  This is deliberately separate from the
sparse CUDA scanner: the result is a complete shell enumeration and can be
used to identify exactly which vectors the sparse parametrization misses.
"""
from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def paley_conference(p: int) -> np.ndarray:
    return np.rint(paley_conference_prime_power(p)).astype(np.int64)


def gp_matrix(matrix: np.ndarray) -> str:
    rows = [",".join(str(int(value)) for value in row) for row in matrix]
    return "[" + ";".join(rows) + "]"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--scaled-bound", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gp", default="gp")
    parser.add_argument("--pari-stack", default="2G")
    parser.add_argument("--max-half-vectors", type=int, default=1_000_000)
    args = parser.parse_args()

    p = args.p
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    program = f"""
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
Q={2 * p}*G^-1;
if(denominator(Q)!=1,error("nonintegral scaled dual Gram"));
DM=qfminim(Q,{args.scaled_bound},{args.max_half_vectors});
V=DM[3];
print("COUNT=",DM[1]);
print("STORED=",matsize(V)[2]);
print("MAXNORM=",DM[2]);
for(j=1,matsize(V)[2],k=V[,j];s=k~*Q*k;unum=B*Q*k;print("VECTOR=",s,"|",Vec(k),"|",Vec(unum)));
quit;
"""
    started = time.monotonic()
    proc = subprocess.run(
        [args.gp, "-fq", "-s", args.pari_stack],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    elapsed = time.monotonic() - started
    if "***" in proc.stderr:
        raise RuntimeError(f"PARI/GP failed: {proc.stderr.strip()}")

    signed_count = stored = maximum_norm = None
    norms: list[int] = []
    coordinates: list[list[int]] = []
    numerators: list[list[int]] = []
    for line in proc.stdout.splitlines():
        if line.startswith("COUNT="):
            signed_count = int(line.split("=", 1)[1])
        elif line.startswith("STORED="):
            stored = int(line.split("=", 1)[1])
        elif line.startswith("MAXNORM="):
            maximum_norm = int(line.split("=", 1)[1])
        elif line.startswith("VECTOR="):
            label, raw_coordinate, raw_numerator = line.split("|", 2)
            norms.append(int(label.split("=", 1)[1]))
            coordinates.append(ast.literal_eval(raw_coordinate))
            numerators.append(ast.literal_eval(raw_numerator))

    if signed_count is None or stored is None or maximum_norm is None:
        raise RuntimeError(f"could not parse PARI summary: {proc.stdout[-2000:]}")
    if signed_count != 2 * stored or stored != len(norms):
        raise RuntimeError(
            "qfminim output was truncated or malformed: "
            f"signed={signed_count}, stored={stored}, parsed={len(norms)}"
        )

    norm_array = np.asarray(norms, dtype=np.int64)
    coordinate_array = np.asarray(coordinates, dtype=np.int64)
    numerator_array = np.asarray(numerators, dtype=np.int64)
    exact_squared_numerators = np.einsum(
        "vi,vi->v", numerator_array, numerator_array, dtype=np.int64
    )
    expected_squared_numerators = 2 * p * norm_array
    if not np.array_equal(exact_squared_numerators, expected_squared_numerators):
        bad = int(np.flatnonzero(
            exact_squared_numerators != expected_squared_numerators
        )[0])
        raise ArithmeticError(
            "ambient exact-norm check failed at row "
            f"{bad}: {exact_squared_numerators[bad]} != "
            f"{expected_squared_numerators[bad]}"
        )

    metadata = {
        "experiment": "r1_dual_shell_export",
        "status": "complete_exact_qfminim_enumeration",
        "p": p,
        "dimension": len(C),
        "rank": coordinate_array.shape[1] if stored else (p * p + 1) // 2,
        "scaled_bound": args.scaled_bound,
        "maximum_scaled_norm": maximum_norm,
        "signed_count": signed_count,
        "stored_half_count": stored,
        "elapsed_seconds": elapsed,
        "ambient_denominator": 2 * p,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
        scaled_norm=norm_array,
        dual_coordinate=coordinate_array,
        ambient_numerator=numerator_array,
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
