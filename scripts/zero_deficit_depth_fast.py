#!/usr/bin/env python3
"""Fast parallel enumeration backend for zero_deficit_transform_depth.py.

Reuses the reference Engine and its exact recursion semantics; only the
maximizer enumeration is replaced by a bit-packed popcount kernel and
multiprocessing across chunk ranges. Validated to reproduce the reference
depths on the pinned witnesses (see --selftest).

Usage:
    python3 zero_deficit_depth_fast.py MATRIX.json [--workers 84] [--chunk 16777216]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zero_deficit_transform_depth import Engine as RefEngine, load_matrix  # noqa: E402

_FULL64 = np.uint64(0xFFFFFFFFFFFFFFFF)


def _pack_rows(B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pack +/-1 rows and diagonal masks as uint64 (bit j = 1 iff entry is -1)."""
    m = B.shape[0]
    rows = np.zeros(m, dtype=np.uint64)
    for i in range(m):
        mask = np.uint64(0)
        for j in range(m):
            if j != i and B[i, j] == -1:
                mask |= np.uint64(1) << np.uint64(j)
        rows[i] = mask
    nd = np.array([_FULL64 ^ (np.uint64(1) << np.uint64(i)) for i in range(m)],
                  dtype=np.uint64)
    return rows, nd


def _q_range(z: np.ndarray, rows: np.ndarray, nd: np.ndarray, m: int) -> np.ndarray:
    """Exact Q_B(x) for packed candidate masks z (bit i = 1 iff spin i = -1)."""
    total = np.zeros(z.shape[0], dtype=np.int64)
    for i in range(m):
        cnt = np.bitwise_count((z ^ rows[i]) & nd[i]).astype(np.int64)
        deg = np.int64(m - 1) - 2 * cnt
        s_i = 1 - 2 * ((z >> np.uint64(i)) & np.uint64(1)).astype(np.int64)
        total += s_i * deg
    return total // 2


def _chunk_scan(payload) -> tuple[int, list, int, list]:
    rows, nd, m, lo, hi = payload
    z = np.arange(lo, hi, dtype=np.uint64) << np.uint64(1)
    q = _q_range(z, rows, nd, m)
    qmax = int(q.max())
    qmin = int(q.min())
    pos = z[q == qmax]
    neg = z[q == qmin]
    return qmax, [int(v) for v in pos.tolist()], qmin, [int(v) for v in neg.tolist()]


class FastEngine(RefEngine):
    def __init__(self, A: np.ndarray, chunk: int, gpu: bool = False,
                 workers: int = 1):
        super().__init__(A, chunk, gpu=False)
        self.workers = max(1, int(workers))
        self.pool = None

    def _get_pool(self):
        if self.pool is None and self.workers > 1:
            import multiprocessing as mp
            self.pool = mp.get_context("fork").Pool(self.workers)
        return self.pool

    def close(self):
        if self.pool is not None:
            self.pool.terminate()
            self.pool.join()
            self.pool = None

    def _enumerate_extrema(self, subset: int):
        if subset in self.max_cache:
            return self.max_cache[subset]
        verts = [i for i in range(self.n) if (subset >> i) & 1]
        m = len(verts)
        if m <= 1:
            out = (0, 0, [0], [0])
            self.max_cache[subset] = out
            return out

        B = self.A[np.ix_(verts, verts)].astype(np.int64, copy=False)
        rows, nd = _pack_rows(B)
        total = 1 << (m - 1)

        step = max(int(self.chunk), 1 << 20)
        ranges = [(lo, min(lo + step, total)) for lo in range(0, total, step)]

        pool = self._get_pool() if (self.workers > 1 and m >= 24 and len(ranges) > 1) else None
        if pool is not None:
            results = pool.map(_chunk_scan,
                               [(rows, nd, m, lo, hi) for lo, hi in ranges])
        else:
            results = [_chunk_scan((rows, nd, m, lo, hi)) for lo, hi in ranges]

        P = -10**18
        Qmin = 10**18
        pos_masks: list[int] = []
        neg_masks: list[int] = []
        for qmax, pos, qmin, neg in results:
            if qmax > P:
                P = qmax
                pos_masks = []
            if qmax == P:
                pos_masks.extend(pos)
            if qmin < Qmin:
                Qmin = qmin
                neg_masks = []
            if qmin == Qmin:
                neg_masks.extend(neg)

        out = (int(P), int(-Qmin), pos_masks, neg_masks)
        self.max_cache[subset] = out
        return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", nargs="?", default=None)
    ap.add_argument("--workers", type=int, default=84)
    ap.add_argument("--chunk", type=int, default=1 << 22)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--decisions-only", action="store_true")
    args = ap.parse_args()

    cases = [
        ("evidence/K15_exact_minimizer_20260911.json", 2, None),
        ("evidence/m16_phi30_witness_20260911.json", 2, 2),
        ("evidence/m19_phi39_witness_20260911.json", 2, None),
        ("evidence/coherent_BI_lift_20260919/best_n25.json", 3, 3),
        ("/tmp/paley_n14.json", 2, 2),
        ("/tmp/paley_n18.json", 2, 2),
        ("/tmp/paley_n26.json", 3, 3),
    ]
    if args.selftest:
        ok = True
        for path, want_pos, want_neg in cases:
            A = load_matrix(path)
            E = FastEngine(A, args.chunk, workers=args.workers)
            root = (1 << E.n) - 1
            F, P, N = E.phi(root)
            depths = {}
            for sigma, edge in ((1, P), (-1, N)):
                if edge == F:
                    d, _, _ = E.depth(root, sigma)
                    depths[sigma] = d
                else:
                    depths[sigma] = None
            E.close()
            got_pos, got_neg = depths[1], depths[-1]
            ok_case = ((want_pos is None or got_pos == want_pos)
                       and (want_neg is None or got_neg == want_neg))
            ok = ok and ok_case
            print(f"{path}: n={E.n} F={F} depth(+)= {got_pos} depth(-)= {got_neg} "
                  f"[want {want_pos}/{want_neg}] {'OK' if ok_case else 'MISMATCH'}",
                  flush=True)
        sys.exit(0 if ok else 1)

    if args.matrix is None:
        ap.error("matrix required unless --selftest")

    A = load_matrix(args.matrix)
    n = A.shape[0]
    E = FastEngine(A, args.chunk, workers=args.workers)
    t0 = time.time()
    root = (1 << n) - 1
    F, P, N = E.phi(root)
    r = F * ((1.0 + 1.0 / n) ** 1.5 - 1.0)

    phases = []
    for sigma, name, edge in ((1, "+", P), (-1, "-", N)):
        if edge != F:
            phases.append({
                "phase": name, "root_edge": edge, "root_deficit": F - edge,
                "zero_deficit_present": False, "depth": 0,
                "margin_r_minus_depth": r,
                "safe_for_pure_zero_deficit_branch": None, "chain": [],
            })
            continue
        d, p, ch = E.depth(root, sigma)
        phases.append({
            "phase": name, "root_edge": edge, "root_deficit": 0,
            "zero_deficit_present": True, "depth": d,
            "margin_r_minus_depth": r - d,
            "safe_for_pure_zero_deficit_branch": bool(d < r),
            "best_first_pivot": p,
            "chain": E.witness_chain(root, sigma),
        })
    out = {
        "matrix": args.matrix, "n": n, "Phi": F, "P": P, "N": N,
        "neutral_r": r,
        "criterion": "pure zero-deficit same-phase branch is safe iff minimax depth < r",
        "phases": phases,
        "principal_blocks_evaluated": len(E.max_cache),
        "depth_states_evaluated": len(E.depth_cache),
        "engine": "packed-popcount", "workers": E.workers,
        "seconds": round(time.time() - t0, 2),
    }
    print(json.dumps(out, indent=2))
    E.close()


if __name__ == "__main__":
    main()
