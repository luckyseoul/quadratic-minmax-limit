#!/usr/bin/env python3
"""Data-free check: <W, Phi_part W> = lambda_bar ||W||^2 for random W in Z."""
import numpy as np, itertools, sys
p=int(sys.argv[1]); q=p*p; n=q+1
r=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
def mul(u,v):
    a1,b1=divmod(u,p); a2,b2=divmod(v,p)
    return p*((a1*a2+r*b1*b2)%p)+((a1*b2+a2*b1)%p)
def sub(u,v):
    a1,b1=divmod(u,p); a2,b2=divmod(v,p)
    return p*((a1-a2)%p)+((b1-b2)%p)
sq=set(mul(x,x) for x in range(1,q))
C=np.zeros((n,n),dtype=np.int64); C[0,1:]=1; C[1:,0]=1
for e1 in range(q):
    for e2 in range(q):
        if e1!=e2: C[1+e1,1+e2]= 1 if sub(e1,e2) in sq else -1
Cf=C.astype(np.float64); Pp=(np.eye(n)+Cf/p)/2
lbar=8*(n-2)/(n-6)
S4=np.array(list(itertools.combinations(range(n),4)),dtype=np.int32)
i_,j_,k_,l_=S4[:,0],S4[:,1],S4[:,2],S4[:,3]
kap=(C[i_,j_]*C[k_,l_]+C[i_,k_]*C[j_,l_]+C[i_,l_]*C[j_,k_]).astype(np.float64)
# phi and star vectorized
phi=np.empty(len(S4),dtype=np.float64)
CH=200000
for lo in range(0,len(S4),CH):
    hi=min(lo+CH,len(S4)); a,b,c,d=i_[lo:hi],j_[lo:hi],k_[lo:hi],l_[lo:hi]
    pr=(C[:,a]*C[:,b]*C[:,c]*C[:,d])
    idx=np.arange(hi-lo)
    phi[lo:hi]=pr.sum(axis=0)-(pr[a,idx]+pr[b,idx]+pr[c,idx]+pr[d,idx])
    del pr
star=(C[j_,i_]*C[k_,i_]*C[l_,i_]+C[i_,j_]*C[k_,j_]*C[l_,j_]
     +C[i_,k_]*C[j_,k_]*C[l_,k_]+C[i_,l_]*C[j_,l_]*C[k_,l_]).astype(np.float64)
D=p*p*(p*p-5)
m4p=((p*p-1)*kap-2*phi-2*p*star)/D
del kap,phi,star
rng=np.random.default_rng(0)
print(f"p={p} n={n} four-sets={len(S4)} lbar={lbar:.9f}",flush=True)
for trial in range(3):
    M=rng.standard_normal((n,n)); M=(M+M.T)/2
    W=Pp@M@Pp
    for _ in range(200):                      # alternating projection onto Z
        np.fill_diagonal(W,0.0); W=Pp@W@Pp
    np.fill_diagonal(W,0.0)
    res=max(np.abs(Cf@W-p*W).max(), np.abs(np.diag(W)).max())
    nrm=(W*W).sum()
    t=(W[i_,j_]*W[k_,l_]+W[i_,k_]*W[j_,l_]+W[i_,l_]*W[j_,k_])
    qf=8*np.sum(m4p*t)+2*nrm+4*np.sum((W@W)*(Cf/p))
    print(f"  trial{trial}: Z-residual={res:.2e}  <W,Phi_part W>/||W||^2 = {qf/nrm:.9f}"
          f"  dev={abs(qf/nrm-lbar):.2e}",flush=True)
