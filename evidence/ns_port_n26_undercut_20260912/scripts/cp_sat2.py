import json, sys, time
import numpy as np
from ortools.sat.python import cp_model

N = 26
d = np.load("live43.npz")
Q = d["Q"].astype(np.int64); X = d["X"].astype(np.int64)
A = np.array(json.load(open(
    "/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/"
    "inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json"))["A"], dtype=np.int64)
edges = [(i, j) for i in range(N) for j in range(i + 1, N)]
M = len(edges)
av = np.array([A[i, j] for (i, j) in edges], dtype=np.int64)
W = np.empty((X.shape[0], M), dtype=np.int64)
for k, (i, j) in enumerate(edges):
    W[:, k] = X[:, i] * X[:, j]
del X

def build_rows(K):
    live = np.abs(Q) >= 60 - 2 * K
    Ql = Q[live]
    sig = np.where(Ql > 0, 1, -1)
    pos = (sig[:, None] * W[live] * av[None, :]) > 0   # aligned edges bool
    tau = (np.abs(Ql) == 61).astype(np.int64)
    # dedupe rows by (pos-bits, tau): keep max tau
    pack = np.packbits(pos, axis=1)                    # (L, 41) uint8
    seen = {}
    for r in range(pack.shape[0]):
        key = (pack[r].tobytes(), int(tau[r]))
        if key not in seen:
            seen[key] = None
    # rows as frozenset of aligned indices
    rows = []
    for (key, tr) in seen:
        bits = np.unpackbits(np.frombuffer(key, dtype=np.uint8))[:M].astype(bool)
        idx = np.flatnonzero(bits)
        need = (K + tr + 1) // 2
        rows.append((idx, need))
    return rows

def run_k(K, time_limit, workers=64):
    t0 = time.time()
    rows = build_rows(K)
    print(f"K={K}: unique rows = {len(rows)}  [{time.time()-t0:.0f}s build]", flush=True)
    model = cp_model.CpModel()
    y = [model.NewBoolVar(f"y{e}") for e in range(M)]
    model.Add(sum(y) == K)
    t0 = time.time()
    for idx, need in rows:
        model.Add(sum(y[e] for e in idx) >= need)
    print(f"  model built [{time.time()-t0:.0f}s]", flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    t0 = time.time()
    status = solver.Solve(model)
    st = solver.StatusName(status)
    print(f"  status: {st}  [{time.time()-t0:.0f}s]  conflicts={solver.NumConflicts()}", flush=True)
    if st in ("OPTIMAL", "FEASIBLE"):
        chosen = [e for e in range(M) if solver.Value(y[e]) == 1]
        print("  EDGES:", [edges[e] for e in chosen], flush=True)
        np.save(f"undercut_K{K}.npy", np.array([edges[e] for e in chosen]))
    return st

if __name__ == "__main__":
    for K in [int(x) for x in sys.argv[1:]]:
        tl = {4: 600, 5: 900, 6: 1500, 7: 2400, 8: 3600}.get(K, 900)
        run_k(K, tl)
