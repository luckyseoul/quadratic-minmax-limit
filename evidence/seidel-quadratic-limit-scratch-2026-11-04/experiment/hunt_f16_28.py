#!/usr/bin/env python3
"""Aggressive hunt for an order-16 signing with Phi <= 28."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

N = 16
PROJ = 1 << (N - 1)
OUT = Path(__file__).resolve().parents[1] / "results"


def edge_list():
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


EDGES = edge_list()
E = len(EDGES)


def spins_table() -> np.ndarray:
    x = np.ones((PROJ, N), dtype=np.int8)
    for i in range(1, N):
        bit = 1 << (i - 1)
        x[:, i] = np.where((np.arange(PROJ) & bit) == 0, 1, -1).astype(np.int8)
    return x


X = spins_table()
U = np.array([e[0] for e in EDGES], dtype=np.int16)
V = np.array([e[1] for e in EDGES], dtype=np.int16)
# products X[:,u]*X[:,v] for each edge: shape (E, PROJ)
XV = np.ascontiguousarray((X[:, U] * X[:, V]).T)  # (E, PROJ) wait -> make (PROJ, E)
XP = (X[:, U] * X[:, V]).astype(np.int8)  # (PROJ, E)


def scores_from_edges(edges: np.ndarray) -> np.ndarray:
    # edges shape (E,), values ±1
    return XP @ edges.astype(np.int32)


def phi_from_scores(scores: np.ndarray) -> int:
    return int(np.max(np.abs(scores)))


def flip(edges: np.ndarray, scores: np.ndarray, e: int) -> None:
    old = int(edges[e])
    edges[e] = -old
    scores -= (2 * old) * XP[:, e].astype(np.int32)


def greedy(edges: np.ndarray, scores: np.ndarray) -> int:
    cur = phi_from_scores(scores)
    improved = True
    while improved:
        improved = False
        order = np.random.permutation(E)
        for e in order:
            flip(edges, scores, int(e))
            v = phi_from_scores(scores)
            if v < cur:
                cur = v
                improved = True
                if cur <= 28:
                    return cur
            else:
                flip(edges, scores, int(e))
    return cur


def load_witness(path: Path) -> np.ndarray:
    data = json.loads(path.read_text(encoding="utf-8"))
    a = np.array(data["A"], dtype=np.int8)
    edges = np.array([a[i, j] for i, j in EDGES], dtype=np.int8)
    return edges


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wit_path = Path(
        r"C:\Users\nicho\seidel-quadratic-limit\_ref\quadratic-minmax-limit\evidence\m16_phi30_witness_20260911.json"
    )
    wit = load_witness(wit_path)
    scores = scores_from_edges(wit)
    assert phi_from_scores(scores) == 30

    best_edges = wit.copy()
    best = 30
    t0 = time.time()
    budget = 600.0
    evals = 0
    rng = np.random.default_rng(7)

    print("numpy hunt start", flush=True)
    # Phase 1: many greeds
    for trial in range(200):
        edges = best_edges.copy() if trial % 5 == 0 else rng.choice(np.array([-1, 1], dtype=np.int8), size=E)
        if trial % 5 == 0:
            # scramble k edges from best
            k = int(rng.integers(2, 15))
            idx = rng.choice(E, size=k, replace=False)
            edges[idx] *= -1
        sc = scores_from_edges(edges)
        g = greedy(edges, sc)
        evals += 1
        if g < best:
            best = g
            best_edges = edges.copy()
            print(f"[GREEDY] Phi={best} trial={trial} t={time.time()-t0:.1f}s", flush=True)
            if best <= 28:
                break

    # Phase 2: SA
    edges = best_edges.copy()
    scores = scores_from_edges(edges)
    cur = phi_from_scores(scores)
    T = 5.0
    while time.time() - t0 < budget and best > 28:
        e = int(rng.integers(0, E))
        flip(edges, scores, e)
        v = phi_from_scores(scores)
        evals += 1
        d = v - cur
        if d <= 0 or rng.random() < np.exp(-d / max(T, 1e-9)):
            cur = v
            if v < best:
                best = v
                best_edges = edges.copy()
                print(f"[SA] Phi={best} evals={evals} t={time.time()-t0:.1f}s", flush=True)
                if best <= 28:
                    break
        else:
            flip(edges, scores, e)
        T *= 0.99999
        if T < 0.03:
            T = 4.0
            edges = best_edges.copy()
            scores = scores_from_edges(edges)
            cur = best

    # Assemble matrix
    a = np.zeros((N, N), dtype=np.int8)
    for s, (i, j) in zip(best_edges, EDGES):
        a[i, j] = a[j, i] = s

    payload = {
        "best_Phi": best,
        "target_reached": best <= 28,
        "evals": evals,
        "seconds": round(time.time() - t0, 2),
        "A": a.tolist(),
        "claim": "F(16)<=28 hence F(16)=28 by monotonicity+parity" if best <= 28 else None,
    }
    out = OUT / "experiment5_numpy_hunt.json"
    out.write_text(json.dumps(payload), encoding="utf-8")
    print(f"FINAL Phi={best} target={best<=28} wrote {out}", flush=True)


if __name__ == "__main__":
    main()
