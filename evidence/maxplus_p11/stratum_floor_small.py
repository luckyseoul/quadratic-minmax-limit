#!/usr/bin/env python3
"""Check the Phi floor on each exhaustive profile stratum at p=5,7.

This diagnostic tests whether the floor can be proved separately on every
Aut-invariant profile stratum.  A stratum below 6 kills that route; a high-k
stratum above 6 motivates a weighted-mixture argument instead.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/nick/quadratic-minmax-limit/src")
sys.path.insert(0, "/mnt/storage/e1work/scripts")

from e1_gmin_m4_prop15588 import z_basis  # type: ignore
from stratum_delta_geometry_small import active_counts  # type: ignore


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=int, choices=(5, 7))
    args = parser.parse_args()
    p = args.p
    y = np.load(f"/tmp/maxplus_p{p}.npy")
    y = np.rint(y[y[:, 0] == 1]).astype(np.float64)
    kvals = active_counts(y.astype(np.int8), p)
    bases = z_basis(p)
    report: dict[str, object] = {"p": p, "dimZ": len(bases), "strata": {}}
    for k in sorted(set(map(int, kvals))):
        ys = y[kvals == k]
        qb = np.einsum("ai,tij,aj->at", ys, bases, ys, optimize=True)
        phi = qb.T @ qb / len(ys)
        eig = np.linalg.eigvalsh((phi + phi.T) / 2)
        report["strata"][f"k{k}"] = {
            "count_eps_plus": len(ys),
            "rank": int(np.linalg.matrix_rank(phi, tol=1e-8)),
            "lambda_min": float(eig[0]),
            "lambda_max": float(eig[-1]),
            "lambda_min_ge_6": bool(eig[0] >= 6 - 1e-8),
        }
    output = Path(f"/mnt/storage/e1work/maxplus_p11/stratum_floor_p{p}.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
