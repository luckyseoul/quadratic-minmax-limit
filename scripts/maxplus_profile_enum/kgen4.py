"""Staged general-k (flip solver replaced by flipsolve.flip_solve) enumerator:
Stage A (per outer): row filters, u-join, v-join -> candidate sigma tables.
Stage B (batched numpy per chunk): point tests for flip-free candidates.
Flip candidates: line-sum propagation solver + tiny DFS (no towers).
"""
import numpy as np,itertools
from math import comb
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen3 import prep_subset  # reuse
from flipsolve import flip_solve as _flip_solve_new

def _flip_solve(p,k,sigstd,idxs,fvec,Tm,thi,tlo,out):
    """sigstd: k x p (std lifts, +p at class p-1). Choose fvec[j] flips (p -> -p)
    among idxs[j] per profile so all point sums land in {tlo,thi}."""
    q=Tm.shape[1]
    A=np.zeros(q,dtype=np.int64)
    for j in range(k): A=A+sigstd[j][Tm[j]]
    diff=A-thi
    if (diff%(2*p)).any(): return
    g=diff//(2*p)   # required Xi in {g, g+1}
    if (g>k).any() or (g+1<0).any(): return
    # unknown positions
    unk=[(j,int(s)) for j in range(k) for s in idxs[j]]
    F=sum(fvec)
    # bounds arrays for Xi contributions: assigned values
    fixed={}   # (j,s)->0/1
    # propagation loop
    changed=True
    while changed:
        changed=False
        for j in range(k):
            Fj=F-fvec[j]
            # refined: use current fixed info for other profiles? keep simple line rule:
            for s in idxs[j]:
                if (j,s) in fixed: continue
                line=np.where(Tm[j]==s)[0]
                Gl=int(g[line].sum())
                lo=Gl-Fj; hi=Gl+p-Fj
                # p*xi in [lo,hi]
                can0=(lo<=0<=hi); can1=(lo<=p<=hi)
                if can0 and not can1: fixed[(j,s)]=0; changed=True
                elif can1 and not can0: fixed[(j,s)]=1; changed=True
                elif not can0 and not can1: return
        # weight constraints
        for j in range(k):
            ones=sum(1 for s in idxs[j] if fixed.get((j,int(s)))==1)
            zer=sum(1 for s in idxs[j] if fixed.get((j,int(s)))==0)
            free=[int(s) for s in idxs[j] if (j,int(s)) not in fixed]
            if ones>fvec[j] or ones+len(free)<fvec[j]: return
            if free and ones==fvec[j]:
                for s in free: fixed[(j,s)]=0
                changed=True
            elif free and ones+len(free)==fvec[j]:
                for s in free: fixed[(j,s)]=1
                changed=True
    free=[(j,int(s)) for j in range(k) for s in idxs[j] if (j,int(s)) not in fixed]
    if len(free)>22: raise RuntimeError(f"flip solver: {len(free)} free bits")
    # brute force remaining (small)
    for bits in range(1<<len(free)):
        asg=dict(fixed)
        okw=True
        for bi,(j,s) in enumerate(free): asg[(j,s)]=(bits>>bi)&1
        for j in range(k):
            if sum(asg.get((j,int(s)),0) for s in idxs[j])!=fvec[j]: okw=False;break
        if not okw: continue
        sg=sigstd.copy()
        for (j,s),val in asg.items():
            if val: sg[j,s]-=2*p
        PS=np.zeros(q,dtype=np.int64)
        for j in range(k): PS=PS+sg[j][Tm[j]]
        if ((PS==thi)|(PS==tlo)).all():
            out.append(np.where(PS==thi,1,-1).astype(np.int8))

def enum_chunk(args):
    ctx,lo,hi,eps=args
    p=ctx['p'];q=ctx['q'];k=ctx['k'];deg=ctx['deg'];Tm=ctx['Tm']
    kern=ctx['kern'];UU=ctx['UU'];c0=ctx['c0'];outer_dims=ctx['outer_dims']
    s_ar=np.arange(p,dtype=np.int64)
    spow={d:(s_ar**d)%p for d in range(0,deg+1)}
    two_p=2*p
    thi=(k-1)*eps+p; tlo=(k-1)*eps-p
    u_ar=np.arange(p,dtype=np.int64)
    f0_batch=[]   # sigma tables k x p, no flips
    sols=[]
    for oi in range(lo,hi):
        cc=oi; coeffs={}
        for d,dim in outer_dims:
            code=cc%(p**dim); cc//=(p**dim)
            vec=np.zeros(k,dtype=np.int64)
            for b in range(dim):
                vec=(vec+(code%p)*kern[d][b])%p
                code//=p
            coeffs[d]=vec
        upper=np.zeros((k,p),dtype=np.int64)
        for d in range(2,deg+1):
            upper=(upper+np.outer(coeffs[d],spow[d]))%p
        allow_u=np.zeros((k,p),dtype=bool)
        vlists=[[None]*p for _ in range(k)]
        dead=False
        for j in range(k):
            W=(upper[j][None,:]+np.outer(u_ar,s_ar))%p
            sW=W.sum(1)
            any_u=False
            for u in range(p):
                hist=np.bincount(W[u],minlength=p)
                cs=np.cumsum(hist[::-1])
                lst=[]
                for v in range(p):
                    nwrap=int(cs[v-1]) if v>0 else 0
                    sR=int(sW[u])+v*p-p*nwrap
                    fl=2*sR-p*p+2*p-p*eps
                    if fl<0 or fl%two_p: continue
                    f=fl//two_p
                    cnt=int(hist[(p-1-v)%p])
                    if f>cnt: continue
                    lst.append((v,f))
                if lst:
                    allow_u[j,u]=True; vlists[j][u]=lst; any_u=True
            if not any_u: dead=True; break
        if dead: continue
        mask=allow_u[0,UU[:,0]]
        for j in range(1,k):
            mask&=allow_u[j,UU[:,j]]
        cand=np.where(mask)[0]
        for ci in cand:
            uu=UU[ci]
            vl=[vlists[j][int(uu[j])] for j in range(k)]
            base=[(upper[j]+uu[j]*s_ar)%p for j in range(k)]
            for vcombo in itertools.product(*vl[:-1]):
                vs=sum(v for v,f in vcombo)%p
                vlast=(c0-vs)%p
                hit=[vf for vf in vl[-1] if vf[0]==vlast]
                if not hit: continue
                full=list(vcombo)+[hit[0]]
                R=np.stack([(base[j]+full[j][0])%p for j in range(k)],0)
                sigstd=2*R-p+2
                fvec=[full[j][1] for j in range(k)]
                if sum(fvec)==0:
                    f0_batch.append(sigstd.astype(np.int16))
                else:
                    idxs=[np.where(R[j]==p-1)[0] for j in range(k)]
                    _flip_solve(p,k,sigstd,idxs,fvec,Tm,thi,tlo,sols)
    # stage B: batched point test for flip-free candidates
    if f0_batch:
        S=np.stack(f0_batch)               # M x k x p
        PS=np.zeros((len(S),q),dtype=np.int32)
        for j in range(k):
            PS+=S[:,j,:][:,Tm[j]]
        ok=((PS==thi)|(PS==tlo)).all(axis=1)
        for r in np.where(ok)[0]:
            sols.append(np.where(PS[r]==thi,1,-1).astype(np.int8))
    # activity filter
    out=[]
    for y in sols:
        act=True
        for j in range(k):
            pr=np.array([int(y[Tm[j]==s].sum()) for s in range(p)])
            if (pr==eps).all(): act=False;break
        if act: out.append(y)
    return ctx['subset'],lo,hi,out

_flip_solve_old=_flip_solve
_flip_solve=_flip_solve_new
