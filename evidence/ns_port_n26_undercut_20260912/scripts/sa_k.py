import numpy as np, json, time
d = np.load("stage01.npz")
Q = d["liveQ"].astype(np.int32); X = d["liveX"].astype(np.int32)
N=26
edges = [(i,j) for i in range(N) for j in range(i+1,N)]
A = np.array(json.load(open("/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json"))["A"], dtype=np.int32)
av = np.array([A[i,j] for (i,j) in edges], dtype=np.int32)
M = len(edges)
W = np.empty((M, Q.size), dtype=np.int32)
for k,(i,j) in enumerate(edges): W[k] = X[:,i]*X[:,j]
def score(S):
    return int(np.abs(Q - 2*np.sum(av[S,None]*W[S], axis=0)).max())
rng = np.random.default_rng(12345)
for K in [4,5,6,7,8]:
    best = 10**9; t0=time.time()
    for trial in range(200):
        S = rng.choice(M, size=K, replace=False)
        cur = score(S); it=0; stall=0
        while it < 3000 and stall < 1500:
            i = rng.integers(K)
            new = rng.integers(M)
            if new in S: continue
            T = S.copy(); T[i] = new
            v = score(T)
            if v <= cur: cur = v; S = T; stall = 0
            else: stall += 1
            it += 1
        best = min(best, cur)
    print(f"K={K}: SA best live-max = {best}  [{time.time()-t0:.0f}s]", flush=True)
