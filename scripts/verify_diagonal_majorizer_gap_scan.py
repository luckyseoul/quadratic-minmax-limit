"""Independent structural replay for the finite canonical-gap scan.

This verifies exact source and winning-matrix identities plus numerical
primal residual/counter invariants. It does not certify SDP optimality.
"""
import gzip
import itertools
import json
from pathlib import Path
import sys


path = Path(sys.argv[1])
reader = gzip.open if path.suffix == ".gz" else open
with reader(path, "rt") as f:
    data = json.load(f)
assert data["n"] == 4 and data["cross_blocks"] == 65536
n = 4
source = data["source"]
edges = [(i, j) for i in range(n) for j in range(i + 1, n)]

def phi(matrix):
    return max(abs(sum(matrix[i][j] * x[i] * x[j] for i, j in edges))
               for x in [(1,) + tail for tail in itertools.product((-1, 1), repeat=n - 1)])

assert phi(source) == 4
for bits in range(64):
    matrix = [[0] * n for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        matrix[i][j] = matrix[j][i] = -1 if (bits >> e) & 1 else 1
    assert phi(matrix) >= 4

rows = data["rows"]
best = rows[0]
mask = best["mask"]
B = [[-1 if (mask >> (i * n + j)) & 1 else 1 for j in range(n)] for i in range(n)]
K = [source[i] + B[i] for i in range(n)] + [
    [B[j][i] for j in range(n)] + [-source[i][j] for j in range(n)] for i in range(n)
]
K2 = [[sum(K[i][h] * K[h][j] for h in range(8)) for j in range(8)] for i in range(8)]
I = [[int(i == j) for j in range(8)] for i in range(8)]
def mul(X, Y):
    return [[sum(X[i][h] * Y[h][j] for h in range(8)) for j in range(8)] for i in range(8)]
P = mul([[K2[i][j] - I[i][j] for j in range(8)] for i in range(8)],
        [[K2[i][j] - 9 * I[i][j] for j in range(8)] for i in range(8)])
assert not any(map(any, P))
assert sum(K2[i][i] for i in range(8)) == 56
assert min(min(x["min_eig_D_minus_K"], x["min_eig_D_plus_K"]) for x in rows) >= -5e-7
assert sum(abs(x["relative_gap"]) < 1e-6 for x in rows) == 0
assert sum(abs(x["relative_gap"] - best["relative_gap"]) < 1e-7 for x in rows) == 104
print("PASS exact m4; winning K polynomial; numerical scan residual and count invariants")
