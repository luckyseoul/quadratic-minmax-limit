#!/usr/bin/env python3
"""Which conjugacy classes share a Gamma value? Full enumeration at p=5.
Classify each g by: order, cycle type on P^1, fixed points, sign pattern.
Test whether Gamma is constant on the classical element families."""
import numpy as np, sys
from collections import deque, Counter, defaultdict
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from e1_gmin_m4_prop15590 import MuLab, signed_generators
p=int(sys.argv[1])
lab=MuLab(p,with_deg6=False); n=lab.n; C=lab.C; Y=lab.Yp.astype(np.int32)
gens=signed_generators(p,C); auto=[(pi,d) for pi,d,s in gens if s==1]
idp=(tuple(range(n)),tuple([1]*n)); seen={idp}; dq=deque([idp]); elems=[idp]
while dq:
    pi0,d0=dq.popleft(); pa=np.array(pi0); da=np.array(d0)
    for pi,d in auto:
        e=(tuple(pi[pa]),tuple(da*d[pa]))
        if e not in seen: seen.add(e); dq.append(e); elems.append(e)
print(f"p={p} |G+|={len(elems)}",flush=True)
eidx={e:i for i,e in enumerate(elems)}
# Gamma
gam=np.empty(len(elems))
for i,(pi,d) in enumerate(elems):
    pa=np.array(pi); da=np.array(d,dtype=np.int32)
    s=(Y*(Y[:,pa]*da[None,:])).sum(axis=1).astype(np.float64)
    gam[i]=(s*s).mean()-2*n
# true conjugacy classes via orbit under conjugation
lab_cls=np.full(len(elems),-1,dtype=np.int64); ncls=0
for i in range(len(elems)):
    if lab_cls[i]>=0: continue
    cid=ncls; ncls+=1; stack=[i]; lab_cls[i]=cid
    while stack:
        j=stack.pop()
        pg,dg=map(np.array,elems[j])
        for pi,d in auto:
            ipi=np.argsort(pi); dinv=d[ipi]
            a_p=pg[pi]; a_d=d*dg[pi]
            c_p=ipi[a_p]; c_d=a_d*dinv[a_p]
            k=eidx.get((tuple(c_p),tuple(c_d)))
            if k is not None and lab_cls[k]<0:
                lab_cls[k]=cid; stack.append(k)
print(f"conjugacy classes: {ncls}",flush=True)
# per-class: Gamma, size, order, cycle type, fixed pts
rows=[]
for c in range(ncls):
    ix=np.flatnonzero(lab_cls==c); i0=int(ix[0])
    pa=np.array(elems[i0][0]); da=np.array(elems[i0][1])
    # order of the signed element
    o=1; cp,cd=pa.copy(),da.copy()
    while not (np.array_equal(cp,np.arange(n)) and np.array_equal(cd,np.ones(n,dtype=int))):
        cp2=pa[cp]; cd2=cd*da[cp]; cp,cd=cp2,cd2; o+=1
        if o>200: break
    # cycle type of the permutation
    ct=[]; vis=np.zeros(n,bool)
    for st in range(n):
        if vis[st]: continue
        L=0; x=st
        while not vis[x]: vis[x]=True; x=pa[x]; L+=1
        ct.append(L)
    rows.append((round(gam[i0],6), len(ix), o, tuple(sorted(Counter(ct).items())),
                 int((pa==np.arange(n)).sum()), int((da==-1).sum())))
byg=defaultdict(list)
for r in rows: byg[r[0]].append(r)
print(f"\ndistinct Gamma values: {len(byg)}  (over {ncls} classes)")
print(f"{'Gamma':>12} {'#classes':>9} {'tot size':>9}   (order, cycletype, fix, #neg) per class")
for g in sorted(byg):
    rs=byg[g]
    print(f"{g:12.4f} {len(rs):9d} {sum(r[1] for r in rs):9d}   "
          f"{[(r[2],r[4],r[5]) for r in rs][:6]}")
