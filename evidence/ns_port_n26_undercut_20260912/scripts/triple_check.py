import numpy as np, time, json
d = np.load("stage01.npz")
Q = d["liveQ"].astype(np.int32); X = d["liveX"].astype(np.int32)
N=26
edges = [(i,j) for i in range(N) for j in range(i+1,N)]
A = np.array(json.load(open("/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json"))["A"], dtype=np.int32)
av = np.array([A[i,j] for (i,j) in edges], dtype=np.int32)
M = len(edges)
W = np.empty((M, Q.size), dtype=np.int32)
for k,(i,j) in enumerate(edges): W[k] = X[:,i]*X[:,j]
best = 10**6; arg = None; t0=time.time(); count=0
best_triple_min = {}
for k in range(M):
    Vk = Q - 2*av[k]*W[k]
    for f in range(k+1, M):
        V = Vk - 2*av[f]*W[f]
        if f+1 <= M-1:
            Gm = np.abs(V[None,:] - 2*av[f+1:,None]*W[f+1:]).max(axis=1)
            count += Gm.size
            m = int(Gm.min())
            if m < best:
                best = m; arg = (k, f, f+1+int(Gm.argmin()))
                print(f"new best {best} at {edges[k]} {edges[f]} {edges[arg[2]]}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE best =", best, "at", [edges[arg[0]], edges[arg[1]], edges[arg[2]]] if arg else None, f"| {count} triples | {time.time()-t0:.0f}s", flush=True)
