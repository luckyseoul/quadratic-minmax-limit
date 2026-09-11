"""Cone-ranked integral repair on the finite canonical-majorizer diagnostic.

This consumes the complete *numerical* order-four cross-block table already
recorded by ``diagonal_majorizer_gap_scan.py``.  A state is a cross-block mask.
Its residual is the pair (relative SDP gap, diagonal dispersion).  At each
stage, all one-bit changes supply finite-difference residual directions; the
positive-cone primitive ranks integral one- and two-bit proposals, and the
correction-cycle driver recomputes the complete residual after each accepted
move.  No fractional cross block is ever accepted.

The table is a floating SDP diagnostic, so this program is a method test only:
it cannot certify an SDP optimum or an all-orders majorizer statement.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np

from correction_cycle.correction_cycle import CorrectionCycle
from stress_cone.stress_cone import NotInConeError, general_cone_representation


TABLE: np.ndarray
TOL = 1.0e-11


def load_table(path: str) -> None:
    global TABLE
    with gzip.open(path, "rt") as handle:
        rows = json.load(handle)["rows"]
    table = np.empty((1 << 16, 2), dtype=float)
    for row in rows:
        diagonal = np.asarray(row["diagonal"], dtype=float)
        total = diagonal.sum()
        dispersion = total * (1.0 / diagonal).sum() / 64.0 - 1.0
        table[int(row["mask"])] = (float(row["relative_gap"]), dispersion)
    if not np.isfinite(table).all():
        raise ValueError("incomplete or non-finite gap table")
    TABLE = table


def better(candidate: np.ndarray, current: np.ndarray) -> bool:
    """Strict lexicographic descent: gap first, then dispersion."""
    if candidate[0] < current[0] - TOL:
        return True
    return abs(candidate[0] - current[0]) <= TOL and candidate[1] < current[1] - TOL


def run_task(task: tuple[int, str]) -> dict:
    initial, mode = task
    last = [initial]
    cone_stages = 0
    evaluated = 0
    accepted_stages = 0
    stalled = [False]

    def residual(mask: int) -> np.ndarray:
        return TABLE[mask]

    def solve(current: np.ndarray, mask: int) -> tuple[int, ...]:
        nonlocal cone_stages, evaluated, accepted_stages
        one = np.array([mask ^ (1 << bit) for bit in range(16)], dtype=np.int64)
        directions = current - TABLE[one]
        candidates: list[tuple[int, ...]] = [(bit,) for bit in range(16)]
        ranking: np.ndarray
        if mode == "cone_pair":
            try:
                # The target is the current two-component defect.  A feasible
                # cone representation means the *linearized* flips can cancel
                # it; only their integer ranking is used below.
                rep = general_cone_representation(directions, current)
                ranking = np.argsort(-rep.coefficients, kind="stable")[:12]
                cone_stages += 1
            except NotInConeError:
                ranking = np.argsort(-directions[:, 0], kind="stable")[:12]
        elif mode == "index_pair":
            ranking = np.arange(12)
        elif mode == "random_pair":
            ranking = np.random.default_rng(9109 + 97 * mask).permutation(16)[:12]
        else:
            ranking = np.empty(0, dtype=int)
        candidates.extend(
            (int(bit_a), int(bit_b))
            for offset, bit_a in enumerate(ranking)
            for bit_b in ranking[offset + 1 :]
        )
        best_mask = mask
        best = current
        for bits in candidates:
            trial = mask
            for bit in bits:
                trial ^= 1 << bit
            candidate = TABLE[trial]
            evaluated += 1
            if better(candidate, best):
                best_mask, best = trial, candidate
        if best_mask == mask:
            stalled[0] = True
            return ()
        accepted_stages += 1
        return tuple(bit for bit in range(16) if (mask ^ best_mask) & (1 << bit))

    def apply(mask: int, bits: tuple[int, ...]) -> int:
        for bit in bits:
            mask ^= 1 << bit
        last[0] = mask
        return mask

    cycle = CorrectionCycle(
        residual_fn=residual,
        linear_solve=solve,
        apply_correction=apply,
        residual_norm=lambda value: -2.0 if stalled[0] else float(value[0]),
        tol=-1.0,
        max_stages=8,
    )
    result = cycle.run(initial)
    final = last[0]
    return {
        "mode": mode,
        "initial": initial,
        "final": final,
        "initial_residual": TABLE[initial].tolist(),
        "final_residual": TABLE[final].tolist(),
        "stages": accepted_stages,
        "cone_stages": cone_stages,
        "evaluated": evaluated,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--starts", type=int, default=65536)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("refusing to overwrite evidence")
    if not 1 <= args.workers <= len(os.sched_getaffinity(0)):
        raise ValueError("invalid worker budget")
    if not 1 <= args.starts <= 65536:
        raise ValueError("invalid start count")
    raw = args.input.read_bytes()
    starts = np.linspace(0, 65535, args.starts, dtype=np.int64).tolist()
    tasks = [(start, mode) for start in starts for mode in ("single_edge", "cone_pair", "index_pair", "random_pair")]
    began = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers, initializer=load_table, initargs=(str(args.input),)) as pool:
        rows = list(pool.map(run_task, tasks, chunksize=128))
    summary = {}
    for mode in ("single_edge", "cone_pair", "index_pair", "random_pair"):
        values = [row for row in rows if row["mode"] == mode]
        initial_gap = np.array([row["initial_residual"][0] for row in values])
        final_gap = np.array([row["final_residual"][0] for row in values])
        summary[mode] = {
            "strict_gap_improvements": int(np.sum(final_gap < initial_gap - TOL)),
            "mean_gap_reduction": float(np.mean(initial_gap - final_gap)),
            "best_final_relative_gap": float(np.min(final_gap)),
            "mean_stages": float(np.mean([row["stages"] for row in values])),
            "cone_linearizations": int(sum(row["cone_stages"] for row in values)),
        }
    args.output.write_text(json.dumps({
        "classification": "finite numerical majorizer-residual repair diagnostic; no SDP or all-orders claim",
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "starts": args.starts,
        "workers": args.workers,
        "seconds": time.monotonic() - began,
        "summary": summary,
        "rows": rows,
    }, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
