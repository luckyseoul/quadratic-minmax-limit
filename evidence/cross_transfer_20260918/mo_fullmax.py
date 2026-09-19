import numpy as np
rng = np.random.default_rng(4)
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
n=13
A=confA(n); ev=np.linalg.eigvalsh(A); a=ev.max(); b=-ev.min()
def buildL(alpha,rho):
    H=np.kron(A,A)-alpha*(np.kron(A,np.eye(n))+np.kron(np.eye(n),A))
    mu=a*b+alpha*(a-b)
    Sig=np.eye(n*n)+rho*H/mu
    return np.linalg.cholesky(Sig+1e-9*np.eye(Sig.shape[0]))
thr=2*(2**0.5)*0.32584
print("threshold/n^1.5 = %.4f (n=%d)" % (thr, n))
idx=np.arange(2**n)[:,None]
X=((idx>>np.arange(n)[None,:])&1*2-1).astype(np.float64)
d=0.5*np.diag(X@(A@X.T))  # Qx for all x
niter=20
Cs=[]
for it in range(niter):
    L=buildL(a, 0.995)
    G=(L.T@rng.normal(size=(n*n,1))).reshape(n,n)
    W=rng.normal(size=(n,n))
    Z=np.sqrt(kappa)*G+np.sqrt(1-kappa)*W
    M=(X@Z)@X.T                     # (2^n,2^n) : M[x,y]=x^T Z y
    dxy=d[:,None]-d[None,:]
    V=np.abs(M + dxy)
    Cs.append(V.max())
Cs=np.array(Cs)
print("full max |xZy + Qx - Qy| / n^1.5:")
print("  mean=%.4f  min=%.4f  p5=%.4f  max=%.4f" % (Cs.mean()/n**1.5, Cs.min()/n**1.5, np.percentile(Cs,5)/n**1.5, Cs.max()/n**1.5))
print("  all below threshold: %s" % (np.all(Cs < thr*n**1.5)))
