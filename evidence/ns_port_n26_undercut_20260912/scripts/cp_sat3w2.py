import sys, time, json
import numpy as np
from ortools.sat.python import cp_model
N = 26
d = np.load("live43.npz")
Q = d["Q"].astype(np.int64); X = d["X"].astype(np.int64)
A = np.array(json.load(open("/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json"))["A"], dtype=np.int64)
edges = [(i, j) for i in range(N) for j in range(i + 1, N)]
M = len(edges)
av = np.array([A[i, j] for (i, j) in edges], dtype=np.int64)
W = np.empty((X.shape[0], M), dtype=np.int64)
for k, (i, j) in enumerate(edges):
    W[:, k] = X[:, i] * X[:, j]
del X
K = int(sys.argv[1]); cap = int(sys.argv[2])
live = np.abs(Q) >= 60 - 2 * K
Ql = Q[live]; sig = np.where(Ql > 0, 1, -1)
pos = (sig[:, None] * W[live] * av[None, :]) > 0
tau = (np.abs(Ql) == 61).astype(np.int64)
need = (K + tau + 1) // 2
order = np.argsort(-tau)   # maximizer rows first
model = cp_model.CpModel()
y = [model.NewBoolVar(f"y{e}") for e in range(M)]
model.Add(sum(y) == K)
t0 = time.time()
for r in order:
    idx = np.flatnonzero(pos[r])
    model.Add(sum(y[e] for e in idx) >= int(need[r]))
print(f"K={K}: live={len(Ql)} rows; model built [{time.time()-t0:.0f}s]", flush=True)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = cap
solver.parameters.num_search_workers = 2
solver.parameters.linearization_level = 2
t0 = time.time()
status = solver.Solve(model)
st = solver.StatusName(status)
print(f"status: {st} [{time.time()-t0:.0f}s] conflicts={solver.NumConflicts()}", flush=True)
if st in ("OPTIMAL", "FEASIBLE"):
    chosen = [e for e in range(M) if solver.Value(y[e]) == 1]
    print("EDGES:", [edges[e] for e in chosen], flush=True)
    np.save(f"undercut_K{K}.npy", np.array([edges[e] for e in chosen]))
