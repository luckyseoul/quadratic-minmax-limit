#!/usr/bin/env python3
import json
from pathlib import Path

import numpy as np

n = 16
proj = 1 << (n - 1)


def sylvester(k: int) -> np.ndarray:
    H = np.array([[1]], dtype=np.int8)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]]).astype(np.int8)
    return H


def phi_matrix(A: np.ndarray) -> int:
    best = 0
    for m in range(proj):
        x = np.ones(n, dtype=np.int8)
        for i in range(1, n):
            x[i] = 1 if ((m >> (i - 1)) & 1) == 0 else -1
        q = int(x @ A @ x) // 2
        best = max(best, abs(q))
    return best


H = sylvester(4)
A = H.copy()
np.fill_diagonal(A, 0)
print("Sylvester-Seidel Phi", phi_matrix(A), flush=True)

EDGES = [(i, j) for i in range(n) for j in range(i + 1, n)]
E = len(EDGES)
X = np.ones((proj, n), dtype=np.int8)
for i in range(1, n):
    bit = 1 << (i - 1)
    X[:, i] = np.where((np.arange(proj) & bit) == 0, 1, -1).astype(np.int8)
U = np.array([e[0] for e in EDGES])
V = np.array([e[1] for e in EDGES])
XP = (X[:, U] * X[:, V]).astype(np.int8)


def phi_edges(edges: np.ndarray) -> int:
    return int(np.max(np.abs(XP @ edges.astype(np.int32))))


rng = np.random.default_rng(42)
best = 10**9
best_edges = None
for trial in range(400):
    edges = rng.choice(np.array([-1, 1], dtype=np.int8), size=E)
    cur = phi_edges(edges)
    for _pass in range(6):
        improved = False
        for e in rng.permutation(E):
            edges[e] *= -1
            v = phi_edges(edges)
            if v < cur:
                cur = v
                improved = True
            else:
                edges[e] *= -1
        if not improved:
            break
    if cur < best:
        best = cur
        best_edges = edges.copy()
        print("best", best, "trial", trial, flush=True)
        if best <= 28:
            break

print("final random-greedy best", best, flush=True)
out = {
    "sylvester_phi": phi_matrix(A),
    "random_greedy_best": int(best),
    "target_28": best <= 28,
}
Path(__file__).resolve().parents[1].joinpath("results", "experiment5_sylvester_hunt.json").write_text(
    json.dumps(out, indent=2), encoding="utf-8"
)
