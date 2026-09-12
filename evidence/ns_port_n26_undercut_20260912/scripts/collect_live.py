import json, numpy as np, cupy as cp, time
N=26
A = np.array(json.load(open("/home/nick/scratch/campaign/nuka_campaign/mo-nuka-source-20260906-C3PMKY/inputs/mo-nuka-20260906-cpu-seeds-v3/v100_paired_n26_phi61.json"))["A"], dtype=np.int64)
Ac = cp.asarray(A, dtype=cp.float32)
Ql=[]; Xl=[]
t0=time.time()
full = 1 << (N-1); chunk = 1 << 21
for lo in range(0, full, chunk):
    q = cp.arange(lo, min(lo+chunk, full), dtype=cp.int64)
    X = cp.empty((q.size, N), dtype=cp.int8); X[:,0]=1
    for c in range(1,N): X[:,c] = (1 - 2*cp.bitwise_and(q >> (c-1), 1)).astype(cp.int8)
    Xf = X.astype(cp.float32)
    Q = cp.sum(cp.matmul(Xf, Ac)*Xf, axis=1)*cp.float32(0.5)
    if float(cp.abs(Q).max()) >= 43:
        Qn = cp.asnumpy(Q); Xn = cp.asnumpy(X)
        m = np.abs(Qn) >= 43
        Ql.append(Qn[m].astype(np.int16)); Xl.append(Xn[m])
Ql = np.concatenate(Ql); Xl = np.concatenate(Xl)
print("collected", len(Ql), "states with |Q|>=43  [", time.time()-t0, "s]")
np.savez_compressed("live43.npz", Q=Ql, X=Xl)
# counts
for v in sorted({abs(int(v)) for v in np.unique(Ql)}, reverse=True):
    c = int((np.abs(Ql)==v).sum())
    print(v, c)
