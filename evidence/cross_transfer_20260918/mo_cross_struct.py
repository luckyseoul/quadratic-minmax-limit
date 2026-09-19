import numpy as np
rng = np.random.default_rng(42)
kappa = 2.0/np.pi
def legendre(q):
    c = {0:0}
    for k in range(1,q):
        t = pow(k,(q-1)//2,q)
        c[k] = 1 if t==1 else -1
    return c
def confA(q):
    chi = legendre(q)
    A = np.zeros((q,q))
    for i in range(q):
        for j in range(q):
            if i!=j: A[i,j] = chi[(i-j)%q]
    return A
def makeSigma(A, alpha, rho):
    n = len(A)
    ev = np.linalg.eigvalsh(A)
    a = ev.max(); b = -ev.min()
    H = np.kron(A,A) - alpha*(np.kron(A,np.eye(n)) + np.kron(np.eye(n),A))
    mu = a*b + alpha*(a-b)
    Sig = np.eye(n*n) + rho*H/mu
    return Sig
def cmax(A, alpha, rho, niter):
    Sig = makeSigma(A,alpha,rho)
    S = Sig.shape[0]
    L = np.linalg.cholesky(Sig + 1e-9*np.eye(S))
    N = len(A)
    idx = np.arange(2**N)[:,None]
    X = (((idx >> np.arange(N)[None,:]) & 1)*2 - 1).astype(np.float64)
    Cs = []
    for _ in range(niter):
        G = (L.T @ rng.normal(size=(S,))).reshape(N,N)
        W = rng.normal(size=(N,N))
        Z = np.sqrt(kappa)*G + np.sqrt(1-kappa)*W
        M = X @ Z
        Cs.append(np.abs(M).sum(axis=1).max())
    Cs = np.array(Cs)
    return Cs.mean()/N**1.5, np.percentile(Cs,5)/N**1.5, Cs.max()/N**1.5
for q in [13,17]:
    A = confA(q)
    a = np.linalg.eigvalsh(A).max(); bm = -np.linalg.eigvalsh(A).min()
    for alpha in [0.0, 0.0, (a-bm)/2]:
        m,p5,mx = cmax(A,alpha,0.5,30)
        print(f"n={q} alpha={alpha:.2f} rho=0.5: C_cross mean={m:.4f} p5={p5:.4f} max={mx:.4f}   thr=0.920")
