"""Parallel general-k enumerator; DFS flip-branching with per-point interval pruning."""
import numpy as np,itertools
from math import comb
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen import field_ctx,square_coords,kernel_modp

def prep_subset(p,subset,forms,coords):
    q=p*p;k=len(subset);deg=k-2
    F=[forms[j] for j in subset]; Tm=np.stack([coords[j] for j in subset],0)
    kern={}
    for d in range(2,deg+1):
        A=[[comb(d,i)*pow(F[j][0],i)*pow(F[j][1],d-i) for j in range(k)] for i in range(d+1)]
        B,_=kernel_modp(A,p)
        kern[d]=np.array(B,dtype=np.int64) if B else np.zeros((0,k),dtype=np.int64)
    A1=np.array([[F[j][0] for j in range(k)],[F[j][1] for j in range(k)]],dtype=np.int64)
    U,_=kernel_modp(A1,p)
    U=np.array(U,dtype=np.int64); du=len(U)
    alphas=np.array(list(itertools.product(range(p),repeat=du)),dtype=np.int64)
    UU=(alphas@U)%p if du>0 else np.zeros((1,k),dtype=np.int64)
    c0=((k*p-k-1+p)//2)%p
    outer_dims=[(d,int(kern[d].shape[0])) for d in range(deg,1,-1)]
    outer_total=1
    for d,dim in outer_dims: outer_total*=p**dim
    return dict(p=p,q=q,k=k,deg=deg,Tm=Tm,kern=kern,UU=UU,c0=c0,
                outer_dims=outer_dims,outer_total=outer_total,subset=tuple(subset))

def _final_tests(p,k,upper,uu,vf,Tm,eps,sols):
    """Build sigma tables; DFS over flip choices with per-point interval pruning."""
    s_ar=np.arange(p,dtype=np.int64)
    two_p=2*p
    thi=(k-1)*eps+p; tlo=(k-1)*eps-p
    R=np.stack([(upper[j]+uu[j]*s_ar+vf[j][0])%p for j in range(k)],0)
    sig=2*R-p+2
    fvec=[vf[j][1] for j in range(k)]
    idxs=[np.where(R[j]==p-1)[0] for j in range(k)]
    order=sorted(range(k),key=lambda j:comb(len(idxs[j]),fvec[j]))
    sig_min=sig.copy(); sig_max=sig.copy()
    for j in range(k):
        if fvec[j]>0:
            sig_min[j,idxs[j]]=-p; sig_max[j,idxs[j]]=p
    q=Tm.shape[1]
    remin=np.zeros((k+1,q),dtype=np.int64); remax=np.zeros((k+1,q),dtype=np.int64)
    for lvl in range(k-1,-1,-1):
        j=order[lvl]
        remin[lvl]=remin[lvl+1]+sig_min[j][Tm[j]]
        remax[lvl]=remax[lvl+1]+sig_max[j][Tm[j]]
    def dfs(lvl,PS):
        if lvl==k:
            if ((PS==thi)|(PS==tlo)).all():
                sols.append(np.where(PS==thi,1,-1).astype(np.int8))
            return
        j=order[lvl]
        for pick in itertools.combinations(idxs[j],fvec[j]):
            sg=sig[j].copy()
            for s_ in pick: sg[s_]-=two_p
            PS2=PS+sg[Tm[j]]
            lo=PS2+remin[lvl+1]; hivec=PS2+remax[lvl+1]
            ok=((lo<=thi)&(hivec>=thi))|((lo<=tlo)&(hivec>=tlo))
            if ok.all(): dfs(lvl+1,PS2)
    dfs(0,np.zeros(q,dtype=np.int64))

def enum_chunk(args):
    ctx,lo,hi,eps=args
    p=ctx['p'];q=ctx['q'];k=ctx['k'];deg=ctx['deg'];Tm=ctx['Tm']
    kern=ctx['kern'];UU=ctx['UU'];c0=ctx['c0'];outer_dims=ctx['outer_dims']
    s_ar=np.arange(p,dtype=np.int64)
    spow={d:(s_ar**d)%p for d in range(0,deg+1)}
    two_p=2*p
    u_ar=np.arange(p,dtype=np.int64)
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
            for vcombo in itertools.product(*vl[:-1]):
                vs=sum(v for v,f in vcombo)%p
                vlast=(c0-vs)%p
                hit=[vf for vf in vl[-1] if vf[0]==vlast]
                if not hit: continue
                full=list(vcombo)+[hit[0]]
                _final_tests(p,k,upper,uu,full,Tm,eps,sols)
    out=[]
    for y in sols:
        act=True
        for j in range(k):
            pr=np.array([y[Tm[j]==s].sum() for s in range(p)])
            if (pr==eps).all(): act=False;break
        if act: out.append(y)
    return ctx['subset'],lo,hi,out
