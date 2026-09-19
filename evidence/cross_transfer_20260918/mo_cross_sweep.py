import numpy as np
rng = np.random.default_rng(1)
kappa = 2.0/np.pi
def legendre(q):
    c={0:0}
    for k in range(1,q):
        t=pow(k,(q-1)//2,q); c[k]=1 if t==1 else -1
    return c
def confA(q):
    chi=legendre(q); A=np.zeros((q,q))
    for i in range(q):
        for j in range(q):
            if i!=j: A[i,j]=chi[(i-j)%q]
    return A
def cmax(A,alpha,rho,niter,n):
    ev=np.linalg.eigvalsh(A); a=ev.max(); b=-ev.min()
    H=np.kron(A,A)-alpha*(np.kron(A,np.eye(n))+np.kron(np.eye(n),A))
    mu=a*b+alpha*(a-b)
    Sig=np.eye(n*n)+rho*H/mu
    L=np.linalg.cholesky(Sig+1e-9*np.eye(Sig.shape[0]))
    idx=np.arange(2**n)[:,None]
    X=((idx>>np.arange(n)[None,:])&1*2-1).astype(np.float64)
    Cs=[]
    for _ in range(niter):
        G=(L.T@rng.normal(size=(n*n,1))).reshape(n,n)
        W=rng.normal(size=(n,n))
        Z=np.sqrt(kappa)*G+np.sqrt(1-kappa)*W
        Cs.append(np.abs(X@Z).sum(axis=1).max())
    Cs=np.array(Cs); return Cs.mean()/n**1.5
n=13
A=confA(n)
ev=np.linalg.eigvalsh(A); a=ev.max(); b=-ev.min()
alphas=[0.0,(a-b)/2]
rhos=[0.0,0.1,0.25,0.5,0.75,0.9,0.99]
for alpha in alphas:
    for rho in rhos:
        m=cmax(A,alpha,rho,40,n)
        tag = " YES(<0.92)" if m<0.92 else ""
        print("n=13 alpha=%.2f rho=%.2f: C=%.3f%s" % (alpha,rho,m,tag))
print("=== done ===")
