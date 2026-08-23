#!/usr/bin/env python3
"""Per-factor bad-set sizes: for each in-U class element, which irreducible
f | g divides the content? Compare |B_f|/|inU| with the model 2^{-deg f}."""
import sys,json,time; sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from concurrent.futures import ProcessPoolExecutor
import numpy as np
p=int(sys.argv[1])
def reps():
    out=[]
    for a in range(p):
        for b in range(p):
            for c in range(p):
                if (a*a+b*c)%p==1%p:
                    t=(a,b,c); tn=((-a)%p,(-b)%p,(-c)%p)
                    if t<=tn: out.append((a,b,c,(-a)%p))
    return out
def work(t):
    from e1_gmin_m4_prop15626 import (_switched, _g_factors, _poly_gcd,
        named_z, _mobius_perm, named_gamma, krylov_g)
    from e1_gmin_m4_prop15626 import paley_conference_prime_power
    A,B,C,D=t
    z,bits,eigen,inU,q,mul,add,chi,sig=named_z(p)
    pi=_mobius_perm(p,A,B,C,D); inv=np.empty_like(pi); inv[pi]=np.arange(len(pi))
    y=np.zeros_like(z)
    for j in range(q+1):
        src=int(inv[j])
        if j==0:
            sw=chi(C) if C else 1
            if sw==0: sw=1
        else:
            lin=add(mul(C,j-1),D); sw=chi(lin)
            if sw==0: sw=1
        y[j]=np.int8(int(sw)*int(z[src]))
    Cmat=paley_conference_prime_power(p)
    em=bool(np.max(np.abs(Cmat@y.astype(np.float64)+p*y.astype(np.float64)))<1e-6)
    yb=((1-y)//2).astype(np.uint8)
    if not (em and int(yb[0])==1 and int(yb[1])==0): return None
    omega=None
    gamma,_,_,_=named_gamma(p)
    _,facs=_g_factors(p)
    from e1_gmin_m4_prop15626 import _primitive
    om=_primitive(mul,q); gen=mul(om,om)
    d=(bits^yb)&1
    wfn=d[1:1+q].copy()
    if d[0]: wfn^=1
    c=krylov_g(wfn,gamma,mul,gen,q,(q-1)//2)
    if c is None: return ("NOCONTENT",)
    cl=list(map(int,c))
    return tuple(_poly_gcd(cl,f)!=[1] for f in facs)
if __name__=="__main__":
    from e1_gmin_m4_prop15626 import _g_factors
    _,facs=_g_factors(p)
    degs=[len(f)-1 for f in facs]
    R=reps()
    res=[]
    t0=time.time()
    with ProcessPoolExecutor(max_workers=60) as ex:
        for r in ex.map(work,R,chunksize=4):
            if r is not None: res.append(r)
    noc=sum(1 for r in res if r[0]=="NOCONTENT")
    good=[r for r in res if r[0]!="NOCONTENT"]
    nU=len(good)
    print(f"p={p}: inU={nU} nocontent={noc}  [{time.time()-t0:.0f}s]")
    print(f"{'factor deg':>10} {'model 2^-d':>11} {'|B_f|/inU':>10} {'|B_f|':>6}  ratio")
    out={}
    for i,dg in enumerate(degs):
        bf=sum(1 for r in good if r[i])
        mdl=2.0**-dg
        print(f"{dg:10d} {mdl:11.4f} {bf/nU:10.4f} {bf:6d}  {bf/nU/mdl:6.2f}")
        out[f"f{i}_deg{dg}"]={"B_f":bf,"model":mdl,"rate":bf/nU}
    miss=sum(1 for r in good if any(r))
    print(f"any-factor miss rate: {miss/nU:.4f}  (1-model = {1-np.prod([1-2.0**-d for d in degs]):.4f})")
    json.dump({"p":p,"inU":nU,"factors":out},open(f"/tmp/w2_pf_{p}.json","w"))
