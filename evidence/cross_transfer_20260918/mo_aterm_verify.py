import numpy as np
rng = np.random.default_rng(3)
kappa = 2.0/np.pi
def confA(q):
    c={0:0}
    for k in range(1,q):
        t=pow(k,(q-1)//2,q); c[k]=1 if t==1 else -1
    A=np.zeros((q,q))
    for i in range(q):
        for j in range(q):
            if i!=j: A[i,j]=c[(i-j)%q]
    return A
def buildL(A,alpha,rho,n):
    ev=np.linalg.eigvalsh(A); a=ev.max(); b=-ev.min()
    H=np.kron(A,A)-alpha*(np.kron(A,np.eye(n))+np.kron(np.eye(n),A))
    mu=a*b+alpha*(a-b)
    Sig=np.eye(n*n)+rho*H/mu
    return np.linalg.cholesky(Sig+1e-9*np.eye(Sig.shape[0]))
for n in [13,17]:
    A=confA(n); ev=np.linalg.eigvalsh(A); a=ev.max()
    L=buildL(A, a, 0.995, n)
    idx=np.arange(2**n)[:,None]
    X=((idx>>np.arange(n)[None,:])&1*2-1).astype(np.float64)
    Qxs=[]; Qys=[]
    niter=30 if n==13 else 15
    for _ in range(niter):
        G=(L.T@rng.normal(size=(n*n,1))).reshape(n,n)
        W=rng.normal(size=(n,n))
        Z=np.sqrt(kappa)*G+np.sqrt(1-kappa)*W
        XZ=np.abs(X@Z).sum(axis=1)
        iarg=np.argmax(XZ)
        xstar=X[iarg]
        Zrow=xstar@Z
        ystar=np.sign(Zrow)
        ystar[ystar==0]=1.0
        Qxs.append(0.5*(xstar@(A@xstar)))
        Qys.append(0.5*(ystar@(A@ystar)))
    Qxs=np.array(Qxs); Qys=np.array(Qys)
    print("n=%d alpha=a rho=0.995:" % n)
    print("  Qx*/n: mean=%.3f p95=%.3f max=%.3f" % (Qxs.mean()/n, np.percentile(np.abs(Qxs),95)/n, np.abs(Qxs).max()/n))
    print("  Qy*/n: mean=%.3f p95=%.3f max=%.3f" % (Qys.mean()/n, np.percentile(np.abs(Qys),95)/n, np.abs(Qys).max()/n))
    print("  |Qx*-Qy*|/n: p95=%.3f max=%.3f" % (np.percentile(np.abs(Qxs-Qys),95)/n, np.abs(Qxs-Qys).max()/n))
