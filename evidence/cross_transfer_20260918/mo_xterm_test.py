import numpy as np
rng = np.random.default_rng(42)
def makeA(n):
    A = np.zeros((n,n),dtype=np.int8)
    A[np.triu_indices(n,1)] = rng.choice([-1,1],size=n*(n-1)//2)
    A = A + A.T; return A.astype(np.float32)
for n in [16,20]:
    N = 2**n
    idx = np.arange(N)[:,None]
    X = (((idx >> np.arange(n)[None,:]) & 1) * 2 - 1).astype(np.float32)
    A = makeA(n)
    Qs=[]
    niter = 40 if n==16 else 20
    for _ in range(niter):
        Z = rng.choice([-1,1],size=(n,n)).astype(np.float32)
        vals = np.abs(X @ Z).sum(axis=1)
        xstar = X[int(np.argmax(vals))]
        Qs.append(0.5*(xstar @ (A @ xstar)))
    Qs=np.array(Qs)
    B=[]
    for _ in range(3000):
        x = rng.choice([-1,1],size=n).astype(np.float32); B.append(0.5*(x@(A@x)))
    B=np.array(B)
    print(f'n={n}:  x*=argmax||xZ||_1:  meanQ/n={Qs.mean()/n:.3f} p95/n={np.percentile(Qs,95)/n:.3f} maxQ/n={Qs.max()/n:.3f}  (n^{3/2}/n=sqrt(n)={np.sqrt(n):.2f})')
    print(f'         random x:         mean|Q|/n={np.abs(B).mean()/n:.3f} stdQ/n={B.std()/n:.3f}')
