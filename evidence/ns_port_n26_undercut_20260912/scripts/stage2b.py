import numpy as np, time
N=26
d = np.load("stage01.npz")
A = d["A"]; liveQ = d["liveQ"].astype(np.int64); liveX = d["liveX"].astype(np.int64)
edges = [(i,j) for i in range(N) for j in range(i+1,N)]
M = len(edges)
L = liveX; Q = liveQ
W = np.empty((M, L.shape[0]), dtype=np.int64)
for k,(i,j) in enumerate(edges): W[k] = L[:,i]*L[:,j]
av = np.array([A[i,j] for (i,j) in edges], dtype=np.int64)
S = np.empty((M,M), dtype=np.int64)
for k in range(M):
    V = Q - 2*av[k]*W[k]
    S[k] = np.abs(V[None,:] - 2*av[:,None]*W).max(axis=1)
# invalidate diagonal
S[np.arange(M), np.arange(M)] = 10**6
iu = np.triu_indices(M, k=0)
Sv = S[iu]
order = np.argsort(Sv)
print("distinct best pairs:")
for t in order[:15]:
    k1,k2 = iu[0][t], iu[1][t]
    print(f"  score={Sv[t]}  {edges[k1]} {edges[k2]}")
print("counts by pair score:", {int(s): int((Sv==s).sum()) for s in sorted(set(Sv.tolist()))[:6]})
np.save("pairS.npy", S)
