#!/usr/bin/env python3
"""p=7 GATE: is Gamma_delta = Gamma - lbar*psi_Z quantized?
U_g is a signed permutation, so (U Q U^T)_{ij} = d_{pi^-1 i} d_{pi^-1 j} Q_{pi^-1 i, pi^-1 j}
-- pure reindexing, no matmul."""
import numpy as np, sys, itertools
from collections import deque
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from e1_gmin_m4_prop15590 import MuLab, signed_generators
p=int(sys.argv[1]); NS=int(sys.argv[2]) if len(sys.argv)>2 else 3000
lab=MuLab(p,with_deg6=False); n=lab.n; C=lab.C; Y=lab.Yp.astype(np.int32)
lbar=8*(n-2)/(n-6)
gens=signed_generators(p,C); auto=[(pi,d) for pi,d,s in gens if s==1]
# random-walk sample of G+ (Gamma is a class function, so sampling suffices)
rng0=np.random.default_rng(7)
elems=[(tuple(range(n)),tuple([1]*n))]
pa=np.arange(n); da=np.ones(n,dtype=np.int64)
seen={elems[0]}
while len(elems)<NS:
    pi,d=auto[rng0.integers(len(auto))]
    pa=pi[pa]; da=da*d[pa] if False else da[np.argsort(np.argsort(pa))]*1
    pa2=pi[np.arange(n)]
    # proper composition: new = gen o cur
    cur_p,cur_d=elems[-1]
    cp=np.array(cur_p); cd=np.array(cur_d)
    np_=pi[cp]; nd=cd*d[cp]
    e=(tuple(np_),tuple(nd))
    elems.append(e); seen.add(e)
    pa=np_; da=nd
print(f"p={p} n={n} |Max+|={len(Y)}",flush=True)
# Gamma for ALL elements
elems=list(dict.fromkeys(elems))
print(f"sampled {len(elems)} distinct group elements",flush=True)
gam=np.empty(len(elems))
for i,(pi,d) in enumerate(elems):
    pa=np.array(pi); da=np.array(d,dtype=np.int32)
    s=(Y*(Y[:,pa]*da[None,:])).sum(axis=1).astype(np.float64)
    gam[i]=(s*s).mean()-2*n
print(f"distinct Gamma values: {len(set(np.round(gam,6)))}   Gamma(e)={gam[0]:.1f} (n(n-2)={n*(n-2)})",flush=True)
# orthonormal Z basis
Cf=C.astype(np.float64); pairs=list(itertools.combinations(range(n),2)); A=[]
for (i,j) in pairs:
    E=np.zeros((n,n)); E[i,j]=E[j,i]=1
    A.append((Cf@E-p*E).reshape(-1))
Am=np.array(A).T; sv=np.linalg.svd(Am,compute_uv=False); _,_,Vt=np.linalg.svd(Am)
Q=[]
for v in Vt[np.sum(sv>1e-8):]:
    B=np.zeros((n,n))
    for k,(i,j) in enumerate(pairs): B[i,j]=B[j,i]=v[k]
    for W in Q: B=B-(B*W).sum()*W
    nr=np.linalg.norm(B)
    if nr>1e-8: Q.append(B/nr)
Qa=np.array(Q); dZ=len(Q); print(f"dim Z = {dZ}",flush=True)
# sample elements covering all distinct Gamma values, then psi_Z by reindexing
rng=np.random.default_rng(0)
byval={}
for i,v in enumerate(np.round(gam,6)): byval.setdefault(v,[]).append(i)
samp=[]
for v,ix in byval.items(): samp+= list(rng.choice(ix,size=min(3,len(ix)),replace=False))
samp=sorted(set(samp)); print(f"sampling {len(samp)} elements (>=3 per Gamma value)",flush=True)
psi=np.empty(len(samp))
for t,i in enumerate(samp):
    pi,d=elems[i]; pa=np.array(pi); da=np.array(d,dtype=np.float64)
    ip=np.argsort(pa); s2=da[ip]
    Qg=Qa[:,ip,:][:,:,ip]*(s2[None,:,None]*s2[None,None,:])
    psi[t]=float((Qa*Qg).sum())
gd=gam[samp]-lbar*psi
print(f"psi_Z(e)={psi[list(samp).index(0)] if 0 in samp else 'n/a'}  (dim Z={dZ})",flush=True)
uniq=sorted(set(np.round(gd,5)))
print(f"distinct Gamma_delta values: {len(uniq)}")
print("values:",uniq[:20])
# quantization test against ||delta||^2/24
d2=19180800/1840091; u=d2/24
print(f"\n||delta||^2(p=7) = {d2:.6f}   u = ||delta||^2/24 = {u:.6f}")
r=[v/u for v in uniq]
print("Gamma_delta / u :",[round(x,4) for x in r])
print("ALL INTEGERS:",all(abs(round(x)-x)<1e-3 for x in r))
