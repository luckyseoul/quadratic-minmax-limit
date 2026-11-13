"""Certified windows  L(q) <= Phi(C_{q+1}) <= U(q)  for Paley conference matrices, q prime, q = 1 mod 4.

L: the best explicit state among
     - the torus state (Theorem 4 of the parity-ceiling note):      (q+1) b / 2,
     - the quartic square wave (quartic-resonance note):             (q-1) a / 2 + 1,
     - every state fixed, up to sign, by a subgroup of PSL(2,q) with at most ORB_MAX orbits: cyclic subgroups
       of both tori, A4, S4 (q = +-1 mod 8) and A5 (q = +-1 mod 5); each orbit quotient is maximised exactly
       and the maximiser is rebuilt and re-evaluated on the full matrix.
U: the parity ceiling (integral; with the row identity when the window is at most 12).
Both ends are rigorous; L is a lower bound only, so a wide window says nothing about Phi.

    gcc -O2 -o /tmp/qopt scripts/parity_ceiling/qopt.c -lm
    QOPT=/tmp/qopt python3 scripts/parity_ceiling/windows.py 5 1000 [--json out.json]
"""
from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time

import numpy as np

from ceiling import closed_form, levels_by_class
from paley_lib import chi_table, is_prime, paley, prim_root, two_squares
from sym import _mul, _pow, signed_perm

ORB_MAX = 26
QOPT = os.environ.get("QOPT", "./qopt")
I2=((1,0),(0,1))
def neg(A,q): return tuple(tuple((-v)%q for v in r) for r in A)
def is_pm_identity(A,q): return A==I2 or A==neg(I2,q)
def psl_order(A,q,cap=200):
    P=A
    for k in range(1,cap+1):
        if is_pm_identity(P,q): return k
        P=_mul(P,A,q)
    return None
def inv2(A,q):
    (a,b),(c,d)=A; return ((d%q,(-b)%q),((-c)%q,a%q))
def cyclic_generators(q,chi):
    """yield (label, matrix) for cyclic subgroups of PSL(2,q) of both tori with at most ORB_MAX orbits."""
    g0=prim_root(q)
    for L in range(1,(q-1)//2+1):
        if ((q-1)//2)%L==0 and 2+(q-1)//L<=ORB_MAX:
            t=pow(g0,(q-1)//(2*L),q); yield (f"split C{L}", ((t,0),(0,pow(t,q-2,q))))
    primes=[p for p in range(2,q+2) if (q+1)%p==0 and is_prime(p)]
    M=None
    for tau in range(q):
        if chi[(tau*tau-4)%q]==-1:
            cand=((0,q-1),(1,tau))
            if all(_pow(cand,(q+1)//p,q)!=I2 for p in primes) and _pow(cand,q+1,q)==I2: M=cand; break
    for L in range(1,(q+1)//2+1):
        if ((q+1)//2)%L==0 and (q+1)//L<=ORB_MAX:
            yield (f"nonsplit C{L}", _pow(M,(q+1)//(2*L),q))
def triangle_group(q,k,rng,tries=20000):
    A=((0,q-1),(1,0)); B0=((0,q-1),(1,q-1))
    for _ in range(tries):
        a,b,c=rng.randrange(q),rng.randrange(q),rng.randrange(q)
        if a==0: continue
        d=((1+b*c)*pow(a,q-2,q))%q; g=((a,b),(c,d))
        B=_mul(_mul(g,B0,q),inv2(g,q),q)
        if psl_order(_mul(A,B,q),q,cap=k)==k: return A,B
    return None
def group_quotients(C,gens,q,chi):
    """all Boolean semi-invariant families: for each sign character (eps per generator) with consistent orbits."""
    n=q+1; sp=[signed_perm(g,q,chi) for g in gens]; out=[]
    for mask in range(1<<len(gens)):
        eps=[-1 if (mask>>i)&1 else 1 for i in range(len(gens))]
        val=np.zeros(n,dtype=np.int64); comp=-np.ones(n,dtype=np.int64); ok=True; U=[]
        for P in range(n):
            if comp[P]>=0: continue
            comp[P]=len(U); val[P]=1; stack=[P]
            while stack and ok:
                cur=stack.pop()
                for (pi,s),e in zip(sp,eps):
                    nxt=pi[cur]; v=val[cur]*e*s[cur]
                    if comp[nxt]<0: comp[nxt]=comp[P]; val[nxt]=v; stack.append(nxt)
                    elif val[nxt]!=v: ok=False; break
            if not ok: break
            U.append(None)
        if not ok: continue
        k=len(U); Um=np.zeros((k,n),dtype=np.int64); Um[comp,np.arange(n)]=val
        out.append((eps,Um))
    return out
def exact_max(C,Um):
    R=Um@C@Um.T; m=len(Um); const=int(np.trace(R))//2
    with tempfile.NamedTemporaryFile("w",suffix=".txt",delete=False) as f:
        f.write(f"{m}\n")
        for row in R: f.write(" ".join(map(str,row.tolist()))+"\n")
        fn=f.name
    out=subprocess.run([QOPT,"exact",fn],capture_output=True,text=True).stdout.split(); os.unlink(fn)
    d=dict(t.split("=") for t in out if "=" in t)
    res=[]
    for key,gk in (("max","gmax"),("min","gmin")):
        g=int(d[gk]); gray=g^(g>>1)
        e=np.array([-1 if (gray>>k)&1 else 1 for k in range(m)],dtype=np.int64)
        x=e@Um; val=int(x@C@x)//2
        assert val==int(d[key])+const, (val,d[key],const)
        res.append((val,x))
    return res
def window(q,rng):
    C=paley(q); chi=chi_table(q); n=q+1; a,b=two_squares(q)
    cands=[("torus (Thm 4)",n*b//2,None),("quartic square wave",(q-1)*a//2+1,None)]
    groups=[(lab,[g]) for lab,g in cyclic_generators(q,chi)]
    for k,name,cond in ((3,"A4",True),(4,"S4",q%8 in (1,7)),(5,"A5",q%5 in (1,4))):
        if cond:
            tg=triangle_group(q,k,rng)
            if tg: groups.append((name,list(tg)))
    for lab,gens in groups:
        for eps,Um in group_quotients(C,gens,q,chi):
            if 2<=len(Um)<=ORB_MAX:
                for val,x in exact_max(C,Um):
                    cands.append((f"{lab} ({len(Um)} orbits)",abs(val),x if val>0 else None))
    best=max(cands,key=lambda z:z[1])
    wit=next((x for lab,v,x in cands if v==best[1] and x is not None),None)
    Uint=levels_by_class(q,1,False)[0]["Q"]
    U=Uint
    if Uint-best[1]<=12 and Uint!=best[1]:
        U=levels_by_class(q,1,True)[0]["Q"]
    return {"q":q,"n":n,"a":a,"b":b,"spectral":0.5*n*math.sqrt(q),"L":best[1],"L_source":best[0],"U":U,"U_integral":Uint,
            "closed_form":closed_form(q),"exact":best[1]==U,"witness_x":None if wit is None else "".join('+' if v>0 else '-' for v in wit.tolist())}
if __name__=="__main__":
    lo,hi=int(sys.argv[1]),int(sys.argv[2]); rng=random.Random(20261112); rows=[]
    for q in range(lo,hi,4):
        if not is_prime(q): continue
        t=time.time(); r=window(q,rng); rows.append(r)
        S=r["spectral"]
        print(f"q={q:4d} (a,b)=({r['a']:2d},{r['b']:2d}) S={S:8.1f} L={r['L']:6d} ({r['L']/S:.4f}) U={r['U']:6d} ({r['U']/S:.4f}) {'EXACT' if r['exact'] else 'window '+str(r['U']-r['L'])}  via {r['L_source']}  [{time.time()-t:.0f}s]",flush=True)
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump(rows, f, indent=1)
