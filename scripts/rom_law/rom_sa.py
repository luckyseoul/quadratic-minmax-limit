# Heuristic (simulated annealing) estimates of E max_x x^T J x / n for Haar flat-spectrum J at larger n.
# Lower estimates of the true maxima; validated against exhaustive enumeration at n=30 (same seeds as rom.py).
import numpy as np, sys, json, struct, subprocess, os
from rom import haar
def sa_max(J,R,sweeps,rng):
    n=J.shape[0]; X=rng.choice([-1.0,1.0],size=(R,n)); H=X@J
    for T in np.geomspace(1.2,0.03,sweeps):
        for i in rng.permutation(n):
            d=-4*X[:,i]*(H[:,i]-J[i,i]*X[:,i])          # change of x^T J x when flipping i
            acc=(d>=0)|(rng.random(R)<np.exp(np.minimum(d,0)/T))
            idx=np.nonzero(acc)[0]
            if idx.size:
                H[idx]-=2*X[idx,i][:,None]*J[i][None,:]; X[idx,i]*=-1
    imp=True
    while imp:
        imp=False
        for i in range(n):
            d=-4*X[:,i]*(H[:,i]-J[i,i]*X[:,i]); idx=np.nonzero(d>1e-12)[0]
            if idx.size:
                H[idx]-=2*X[idx,i][:,None]*J[i][None,:]; X[idx,i]*=-1; imp=True
    return np.einsum('ri,ri->r',X,H).max()
out=[]
for n,k in ((30,24),(34,24),(38,24),(42,24),(50,24)):
    vals=[]
    for s in range(k):
        rng=np.random.default_rng(1000*n+s); O=haar(n,rng); d=np.array([1.0]*(n//2)+[-1.0]*(n-n//2)); J=(O*d)@O.T
        r2=np.random.default_rng(7+s)
        vals+= [sa_max(J,300,250,r2)/n, sa_max(-J,300,250,r2)/n]
    v=np.array(vals); r=dict(n=n,instances=k,sa_one_sided_mean=float(v.mean()),se=float(v.std()/np.sqrt(len(v))))
    print(json.dumps(r),flush=True); out.append(r)
json.dump(out,open('rom_sa_results.json','w'),indent=1)
