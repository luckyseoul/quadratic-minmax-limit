#!/usr/bin/env python3
"""Step 1: Gamma(g) = E_y[s(y,g.y)^2] - 2n over the auto subgroup G+.
Enumerate G+ by BFS on (perm, signs); test class-function property and
whether Gamma depends only on the fixed-point count."""
import numpy as np, sys
from collections import deque, Counter
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from e1_gmin_m4_prop15590 import MuLab, paley_conference, signed_generators
p=int(sys.argv[1])
lab=MuLab(p,with_deg6=False); n=lab.n; C=lab.C; Y=lab.Yp.astype(np.int64)
gens=signed_generators(p,C)
auto=[(pi,d) for pi,d,s in gens if s==1]
print(f"p={p} n={n} |Max+|={len(Y)}  auto generators: {len(auto)}",flush=True)
# BFS over G+ : element = (perm tuple, sign tuple)
idp=(tuple(range(n)),tuple([1]*n))
seen={idp}; dq=deque([idp]); elems=[idp]
while dq:
    pi0,d0=dq.popleft()
    pa=np.array(pi0); da=np.array(d0)
    for pi,d in auto:
        # compose: (g after h)(x) = g(h(x)); signs multiply appropriately
        np2=tuple(pi[pa]); nd2=tuple(da*d[pa])
        e=(np2,nd2)
        if e not in seen:
            seen.add(e); dq.append(e); elems.append(e)
    if len(elems)>400000: break
print(f"|G+| = {len(elems)}",flush=True)
# Gamma per element + fixed-point count
Yl=Y.astype(np.int32)
gam={}; fixc={}
for idx,(pi,d) in enumerate(elems):
    pa=np.array(pi); da=np.array(d,dtype=np.int32)
    Z=Yl[:,pa]*da[None,:]
    s=(Yl*Z).sum(axis=1)
    gam[idx]=float((s.astype(np.float64)**2).mean())-2*n
    fixc[idx]=int((pa==np.arange(n)).sum())
vals=np.array([gam[i] for i in range(len(elems))])
fx=np.array([fixc[i] for i in range(len(elems))])
print(f"distinct Gamma values: {len(set(np.round(vals,6)))}",flush=True)
# does Gamma depend only on fixed-point count?
byfix={}
for v,f in zip(np.round(vals,6),fx): byfix.setdefault(f,set()).add(v)
print("fixed-count -> #distinct Gamma:", {k:len(v) for k,v in sorted(byfix.items())})
print("fixed-count -> Gamma values (first 4):",
      {k:sorted(v)[:4] for k,v in sorted(byfix.items())})
# class-function check: Gamma(h g h^-1) == Gamma(g) on random samples
rng=np.random.default_rng(0); bad=0
eidx={e:i for i,e in enumerate(elems)}
for _ in range(300):
    i,j=rng.integers(0,len(elems),2)
    (pg,dg)=elems[i]; (ph,dh)=elems[j]
    pg,dg,ph,dh=map(np.array,(pg,dg,ph,dh))
    phi_inv=np.argsort(ph); dh_inv=dh[phi_inv]
    a_p=pg[ph]; a_d=dh*dg[ph]                    # g o h
    c_p=phi_inv[a_p]; c_d=a_d*dh_inv[a_p]        # h^-1 o g o h
    k=eidx.get((tuple(c_p),tuple(c_d)))
    if k is None: continue
    if abs(gam[k]-gam[int(i)])>1e-6: bad+=1
print(f"class-function violations: {bad}/300",flush=True)

# ---- Gamma_delta = Gamma - lbar * psi_Z ; how structured is it?
import itertools
lbar=8*(n-2)/(n-6)
Cf=C.astype(np.float64)
pairs=list(itertools.combinations(range(n),2)); A=[]
for (i,j) in pairs:
    E=np.zeros((n,n)); E[i,j]=E[j,i]=1
    A.append((Cf@E-p*E).reshape(-1))
Am=np.array(A).T; sv=np.linalg.svd(Am,compute_uv=False); _,_,Vt=np.linalg.svd(Am)
Q=[]
for v in Vt[np.sum(sv>1e-8):]:
    B=np.zeros((n,n))
    for idx2,(i,j) in enumerate(pairs): B[i,j]=B[j,i]=v[idx2]
    for W in Q: B=B-(B*W).sum()*W
    nr=np.linalg.norm(B)
    if nr>1e-8: Q.append(B/nr)
Qa=np.array(Q); dZ=len(Q)
psi=np.zeros(len(elems))
for idx,(pi,d) in enumerate(elems):
    pa=np.array(pi); da=np.array(d,dtype=np.float64)
    Ug=np.zeros((n,n)); Ug[pa,np.arange(n)]=da      # (Ug W Ug^T)
    tr=0.0
    for a in range(dZ):
        Wg=Ug@Qa[a]@Ug.T
        tr+=float((Qa[a]*Wg).sum())
    psi[idx]=tr
gd=vals-lbar*psi
print(f"\npsi_Z(e) = {psi[0]:.1f} (dim Z = {dZ})")
print(f"distinct Gamma_delta values: {len(set(np.round(gd,6)))}")
nz=np.abs(gd)>1e-6
print(f"Gamma_delta nonzero on {nz.sum()}/{len(elems)} elements ({100*nz.mean():.1f}%)")
print(f"Gamma_delta(e) = {gd[0]:.4f}   (= 24||delta||^2 ? -> {24*23.630769:.4f})")
print("distinct Gamma_delta values:", sorted(set(np.round(gd,4))))
