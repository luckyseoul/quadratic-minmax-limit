#!/usr/bin/env python3
"""MILP probes for residual-(ii) at k=4p.

Modes
-----
slice-max-far
    Impose S_G=-2 on the complete Max-minus slice f_01=-1, |G|=4p,
    and e=(0,1) not in G.  Maximise the number of edges disjoint from e.
    An optimum of zero would support sparse double-star rigidity.

official
    Impose the full linear leftover+splus conditions:
      S_G >= 2 on Max-plus,
      S_G <= -2 on Max-minus with f_01=-1,
      S_G <= -4 on Max-minus with f_01=+1,
      |G|=4p and e not in G.

This is a finite diagnostic, not a general-prime proof.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, vstack

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


DEFAULT_PLUS = {
    5: Path("/mnt/storage/e1work/qvar_nuka_caches/maxplus_p5.npy"),
    7: Path("/mnt/storage/e1work/qvar_nuka_caches/maxplus_p7.npy"),
}
DEFAULT_MINUS = {5: Path("/tmp/maxminus_p5.npy"), 7: Path("/tmp/maxminus_p7.npy")}


def feature_rows(Y: np.ndarray, C: np.ndarray) -> tuple[list[tuple[int, int]], np.ndarray]:
    n = C.shape[0]
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    left = np.fromiter((a for a, _b in edges), dtype=np.int32)
    right = np.fromiter((b for _a, b in edges), dtype=np.int32)
    signs = C[left, right].astype(np.int8)
    F = (Y[:, left] * Y[:, right] * signs).astype(np.int8)
    return edges, F


def unique_rows(F: np.ndarray) -> np.ndarray:
    return np.unique(np.ascontiguousarray(F), axis=0)


def solve(
    p: int,
    mode: str,
    plus_path: Path,
    minus_path: Path,
    time_limit: float,
    relax: bool = False,
) -> dict:
    t0 = time.time()
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    Ym = np.sign(np.load(minus_path, mmap_mode="r")).astype(np.int8)
    edges, Fm = feature_rows(Ym, C)
    ei = edges.index((0, 1))
    fe = Fm[:, ei]
    nvar = len(edges)

    rows = []
    lower = []
    upper = []

    # Cardinality and distinguished-edge exclusion.
    rows.append(csr_matrix(np.ones((1, nvar), dtype=np.float64)))
    lower.append(float(4 * p))
    upper.append(float(4 * p))
    erow = np.zeros((1, nvar), dtype=np.float64)
    erow[0, ei] = 1.0
    rows.append(csr_matrix(erow))
    lower.append(0.0)
    upper.append(0.0)

    if mode == "slice-max-far":
        FU = unique_rows(Fm[fe == -1])
        rows.append(csr_matrix(FU.astype(np.float64)))
        lower.extend([-2.0] * len(FU))
        upper.extend([-2.0] * len(FU))
        c = np.zeros(nvar, dtype=np.float64)
        for j, (a, b) in enumerate(edges):
            if a not in (0, 1) and b not in (0, 1):
                c[j] = -1.0
    elif mode == "official":
        Yp = np.sign(np.load(plus_path, mmap_mode="r")).astype(np.int8)
        edges_p, Fp = feature_rows(Yp, C)
        if edges_p != edges:
            raise RuntimeError("plus/minus edge orders differ")
        Fp = unique_rows(Fp)
        FU = unique_rows(Fm[fe == -1])
        FC = unique_rows(Fm[fe == 1])
        rows.extend(
            [
                csr_matrix(Fp.astype(np.float64)),
                csr_matrix(FU.astype(np.float64)),
                csr_matrix(FC.astype(np.float64)),
            ]
        )
        lower.extend([2.0] * len(Fp))
        upper.extend([np.inf] * len(Fp))
        lower.extend([-np.inf] * len(FU))
        upper.extend([-2.0] * len(FU))
        lower.extend([-np.inf] * len(FC))
        upper.extend([-4.0] * len(FC))
        c = np.zeros(nvar, dtype=np.float64)
    else:
        raise ValueError(mode)

    A = vstack(rows, format="csr")
    constraint = LinearConstraint(A, np.asarray(lower), np.asarray(upper))
    result = milp(
        c,
        integrality=(
            np.zeros(nvar, dtype=np.uint8)
            if relax
            else np.ones(nvar, dtype=np.uint8)
        ),
        bounds=Bounds(np.zeros(nvar), np.ones(nvar)),
        constraints=constraint,
        options={
            "time_limit": float(time_limit),
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )

    out = {
        "p": p,
        "mode": mode,
        "relaxation": bool(relax),
        "n_variables": nvar,
        "n_constraints": int(A.shape[0]),
        "status": int(result.status),
        "success": bool(result.success),
        "message": str(result.message),
        "objective": None if result.fun is None else float(result.fun),
        "mip_node_count": getattr(result, "mip_node_count", None),
        "mip_gap": getattr(result, "mip_gap", None),
        "seconds": round(time.time() - t0, 3),
    }
    if result.x is not None and relax:
        x = np.asarray(result.x)
        out.update(
            {
                "sum_x": float(x.sum()),
                "min_x": float(x.min()),
                "max_x": float(x.max()),
                "n_fractional": int(
                    np.count_nonzero((x > 1e-8) & (x < 1.0 - 1e-8))
                ),
            }
        )
    elif result.x is not None:
        x = np.rint(result.x).astype(np.int8)
        chosen = [list(edge) for edge, bit in zip(edges, x) if bit]
        far = sum(a not in (0, 1) and b not in (0, 1) for a, b in map(tuple, chosen))
        out.update(
            {
                "integral_max_error": float(np.max(np.abs(result.x - x))),
                "chosen_edges": chosen,
                "n_chosen": int(x.sum()),
                "n_far": int(far),
                "slice_scores": sorted({int(v) for v in (Fm[fe == -1] @ x).tolist()}),
            }
        )
        if mode == "official":
            out["minus_complement_scores"] = sorted(
                {int(v) for v in (Fm[fe == 1] @ x).tolist()}
            )
            out["plus_scores"] = sorted({int(v) for v in (Fp @ x).tolist()})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=7)
    ap.add_argument("--mode", choices=("slice-max-far", "official"), default="slice-max-far")
    ap.add_argument("--plus", type=Path)
    ap.add_argument("--minus", type=Path)
    ap.add_argument("--time-limit", type=float, default=600.0)
    ap.add_argument("--relax", action="store_true", help="solve the [0,1] LP relaxation")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    plus_path = args.plus if args.plus is not None else DEFAULT_PLUS[args.p]
    minus_path = args.minus if args.minus is not None else DEFAULT_MINUS[args.p]
    out = solve(
        args.p,
        args.mode,
        plus_path,
        minus_path,
        args.time_limit,
        relax=args.relax,
    )
    print(json.dumps(out, indent=2), flush=True)
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
