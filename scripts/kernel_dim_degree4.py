#!/usr/bin/env python3
"""K_4(p): equivariant kernel dim of the complete degree-4 mu/delta contraction
system (Prop 15.590).  Data-free: orbits + sign cocycles only, no Max+/-.
Usage: KP=<prime> python3 scripts/kernel_dim_degree4.py"""
import numpy as np, itertools, time, sys
from collections import deque
from fractions import Fraction
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')

import os
p=int(os.environ["KP"]); q=p*p; n=q+1
t0=time.time()
# F_121 = F_11[t]/(t^2-r), r nonresidue mod 11 -> 2
r=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
a=np.arange(q)//p; b=np.arange(q)%p     # e = p*a+b <-> a+bt
# difference table  e1-e2  encoded, then chi
A1=a[:,None]; B1=b[:,None]; A2=a[None,:]; B2=b[None,:]
DA=(A1-A2)%p; DB=(B1-B2)%p
DEnc=DA*p+DB
# squares of F_121: x^2 = (a^2+r b^2) + (2ab) t
SqEnc=np.unique(((a*a+r*b*b)%p)*p+(2*a*b)%p)
SqEnc=SqEnc[SqEnc!=0] if (0 in SqEnc) else SqEnc
isq=np.zeros(q,dtype=bool); isq[SqEnc]=True; isq[0]=False
C=np.zeros((n,n),dtype=np.int64)
C[0,1:]=1; C[1:,0]=1
blk=np.where(DEnc==0,0,np.where(isq[DEnc],1,-1))
C[1:,1:]=blk
assert (C==C.T).all()
assert (C.astype(np.int64)@C.astype(np.int64)==q*np.eye(n,dtype=np.int64)).all()
print(f"C ok {time.time()-t0:.0f}s")

# field ops for generators
def fmul(e1,e2):
    a1,b1=divmod(e1,p); a2,b2=divmod(e2,p)
    return p*((a1*a2+r*b1*b2)%p)+((a1*b2+a2*b1)%p)
def fadd(e1,e2):
    a1,b1=divmod(e1,p); a2,b2=divmod(e2,p)
    return p*((a1+a2)%p)+((b1+b2)%p)
one=p
finv=[0]*q
for e in range(1,q): finv[e]=next(x for x in range(1,q) if fmul(e,x)==one)
def order_of(e):
    x,o=e,1
    while x!=one: x=fmul(x,e); o+=1
    return o
gen=next(e for e in range(2,q) if order_of(e)==q-1); g2=fmul(gen,gen)
def mkperm(fn,inf_to=0):
    pi=np.zeros(n,dtype=np.int64); pi[0]=inf_to
    for e in range(q): pi[1+e]=fn(e)
    return pi
def frob(e):
    x=e
    for _ in range(p-1): x=fmul(x,e)
    return x
perms=[mkperm(lambda e:1+fadd(e,one)),mkperm(lambda e:1+fadd(e,1)),
 mkperm(lambda e:1+fmul(e,g2)),mkperm(lambda e:1+fmul(e,gen)),
 mkperm(lambda e:(1+finv[e]) if e else 0,inf_to=1),
 mkperm(lambda e:1+frob(e))]
def signed_lift(pi):
    for s in (1,-1):
        d=np.zeros(n,dtype=np.int64); d[0]=1
        d[1:]=s*C[pi[0],pi[1:]]*C[0,1:]
        ok=(s*d[:,None]*d[None,:]*C==C[np.ix_(pi,pi)]); np.fill_diagonal(ok,True)
        if ok.all(): return d,s
    return None,0
gens=[]
for pi in perms:
    d,s=signed_lift(pi); assert s!=0
    gens.append((pi,d,s)); ipi=np.argsort(pi); gens.append((ipi,d[ipi],s))
print(f"gens ok {time.time()-t0:.0f}s")

S4=np.array(list(itertools.combinations(range(n),4)),dtype=np.int16)
mult=np.array([n**3,n*n,n,1],dtype=np.int64)
enc=(S4.astype(np.int64)*mult[None,:]).sum(axis=1)
print(f"{len(S4)} four-sets")

def orbits(twist):
    lab=np.arange(len(S4),dtype=np.int64); sg=np.ones(len(S4),dtype=np.int8)
    dead=np.zeros(len(S4),dtype=bool)
    tgts=[];epss=[];sss=[]
    for pi,d,s in gens:
        img=np.sort(pi[S4.astype(np.int64)],axis=1)
        tgts.append(np.searchsorted(enc,(img*mult[None,:]).sum(axis=1)).astype(np.int32))
        epss.append(d[S4.astype(np.int64)].prod(axis=1).astype(np.int8)); sss.append(s)
    for it in range(300):
        ch=False
        for tgt,eps,s in zip(tgts,epss,sss):
            e=(eps*(s if twist else 1)).astype(np.int8)
            lt=lab[tgt]; st=(sg[tgt]*e).astype(np.int8)
            m=lt<lab
            if m.any(): lab=np.where(m,lt,lab); sg=np.where(m,st,sg); ch=True
            m2=(lt==lab)&(~m)&(st!=sg)
            if m2.any(): dead|=m2
        l2=lab[lab]; s2=(sg[lab]*sg).astype(np.int8)
        if (l2!=lab).any(): lab,sg=l2,s2; ch=True
        if not ch: break
    return lab,sg,set(int(x) for x in np.unique(lab[dead]))
t0=time.time()
labM,sgM,deadM=orbits(False)
labD,sgD,deadD=orbits(True)
uM=np.unique(labM); uD=np.unique(labD)
print(f"orbits {time.time()-t0:.0f}s: mu {len(uM)} (dead {len(deadM)}), delta {len(uD)} (dead {len(deadD)})")

liveM=sorted(set(uM.tolist())-deadM); liveD=sorted(set(uD.tolist())-deadD)
iM={l:i for i,l in enumerate(liveM)}; iD={l:len(liveM)+i for i,l in enumerate(liveD)}
NV=len(liveM)+len(liveD)
print(f"unknowns {NV} = {len(liveM)} mu + {len(liveD)} delta")

def i4(*args):
    x=sorted(args)
    return int(np.searchsorted(enc,x[0]*int(mult[0])+x[1]*int(mult[1])+x[2]*int(mult[2])+x[3]))

# marked-triple reps for star rows
S3=np.array(list(itertools.combinations(range(n),3)),dtype=np.int16)
m3=np.array([n*n,n,1],dtype=np.int64)
e3=(S3.astype(np.int64)*m3[None,:]).sum(axis=1)
lab3=np.arange(len(S3),dtype=np.int64)
tg3=[]
for pi,d,s in gens:
    img=np.sort(pi[S3.astype(np.int64)],axis=1)
    tg3.append(np.searchsorted(e3,(img*m3[None,:]).sum(axis=1)).astype(np.int32))
for it in range(200):
    ch=False
    for tgt in tg3:
        lt=lab3[tgt]; m=lt<lab3
        if m.any(): lab3=np.where(m,lt,lab3); ch=True
    l2=lab3[lab3]
    if (l2!=lab3).any(): lab3=l2; ch=True
    if not ch: break
reps3=[S3[int(np.where(lab3==l)[0][0])] for l in np.unique(lab3)]
print(f"3-set orbits: {len(reps3)}")

rows=set()
# RHS carried exactly so we also get consistency; scale mu by N later irrelevant for kernel
# unknowns are mu (units of mu itself), RHS: (★mu)=C_ij ; (★δ)=-(2/p)C_ik C_jk ; (out) as in prop
def addrow(c,rhs): rows.add((tuple(c),rhs))
FR=Fraction
for T in reps3:
    for k in map(int,T):
        i,j=[int(x) for x in T if x!=k]
        cM=[FR(0)]*NV; cD=[FR(0)]*NV
        for l in range(n):
            if l in (i,j,k) or C[k,l]==0: continue
            si=i4(i,j,k,l)
            lm=int(labM[si]); ld=int(labD[si])
            if lm not in deadM: cM[iM[lm]]+=int(C[k,l])*int(sgM[si])
            if ld not in deadD: cD[iD[ld]]+=int(C[k,l])*int(sgD[si])
        addrow(cM,FR(int(C[i,j])))
        addrow(cD,FR(-2*int(C[i,k])*int(C[j,k]),p))
for l4 in np.unique(labM):
    T=[int(x) for x in S4[int(np.where(labM==l4)[0][0])]]
    for l_ in T:
        tri=[x for x in T if x!=l_]
        st=i4(*T)
        cM=[FR(0)]*NV
        for x in range(n):
            if x in T or C[l_,x]==0: continue
            si=i4(*tri,x); lm=int(labM[si])
            if lm not in deadM: cM[iM[lm]]+=int(C[l_,x])*int(sgM[si])
        ld=int(labD[st])
        if ld not in deadD: cD_=iD[ld]; cM[cD_]+=FR(-p)*int(sgD[st])
        addrow(cM,FR(0))
        cD=[FR(0)]*NV
        for x in range(n):
            if x in T or C[l_,x]==0: continue
            si=i4(*tri,x); ld2=int(labD[si])
            if ld2 not in deadD: cD[iD[ld2]]+=int(C[l_,x])*int(sgD[si])
        lm=int(labM[st])
        if lm not in deadM: cD[iM[lm]]+=FR(-p)*int(sgM[st])
        a_,b_,c_=tri
        corr=int(C[l_,a_])*int(C[b_,c_])+int(C[l_,b_])*int(C[a_,c_])+int(C[l_,c_])*int(C[a_,b_])
        addrow(cD,FR(-corr,p))
print(f"rows: {len(rows)}")
M=[list(rw)+[rhs] for rw,rhs in rows]
rk=0; piv=[]
for c in range(NV):
    pv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
    if pv is None: continue
    M[rk],M[pv]=M[pv],M[rk]; M[rk]=[x/M[rk][c] for x in M[rk]]
    for i2 in range(len(M)):
        if i2!=rk and M[i2][c]!=0: M[i2]=[u-M[i2][c]*v for u,v in zip(M[i2],M[rk])]
    piv.append(c); rk+=1
incons=any(all(x==0 for x in row[:-1]) and row[-1]!=0 for row in M)
print(f"\nK_4({p}): rank={rk}/{NV}, KERNEL DIM={NV-rk}, inconsistent={incons}")
# kappa of each mu-orbit rep for reporting
i_,j_,k_,l_=S4[:,0].astype(np.int64),S4[:,1].astype(np.int64),S4[:,2].astype(np.int64),S4[:,3].astype(np.int64)
kap=C[i_,j_]*C[k_,l_]+C[i_,k_]*C[j_,l_]+C[i_,l_]*C[j_,k_]
free=[c for c in range(NV) if c not in piv]
names=[f"mu(orb{l},|k|={abs(int(kap[int(np.where(labM==l)[0][0])]))})" for l in liveM]+[f"delta(orb{l})" for l in liveD]
for fc in free:
    kv=[FR(0)]*NV; kv[fc]=FR(1)
    for ri,c in enumerate(piv): kv[c]=-M[ri][fc]
    touch=[names[i5] for i5,v in enumerate(kv) if v!=0 and i5<len(liveM)]
    print(f"  free {names[fc]}: mu-coords touched: {len(touch)} -> {touch[:6]}{'...' if len(touch)>6 else ''}")
print("orbit inventory: mu orbits by |kappa|:", 
      {kk:sum(1 for l in liveM if abs(int(kap[int(np.where(labM==l)[0][0])]))==kk) for kk in (1,3)})
