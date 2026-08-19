"""Numba-JIT inner join for the staged general-k enumerator (v2: split prep)."""
import numpy as np,itertools
from numba import njit
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen3 import prep_subset
from kgen4 import _flip_solve

@njit(cache=True)
def _prep_tables(p,k,upper,eps):
    two_p=2*p
    bases=np.zeros((k,p,p),np.int64)
    for j in range(k):
        for u in range(p):
            for s in range(p):
                bases[j,u,s]=(upper[j,s]+u*s)%p
    av=np.full((k,p,p),-1,np.int64)
    af=np.zeros((k,p,p),np.int64)
    an=np.zeros((k,p),np.int64)
    aull=np.zeros((k,p),np.bool_)
    for j in range(k):
        for u in range(p):
            hist=np.zeros(p,np.int64)
            sW=0
            for s in range(p):
                w=bases[j,u,s]; hist[w]+=1; sW+=w
            for v in range(p):
                nwrap=0
                for w in range(p-v,p):
                    nwrap+=hist[w]
                sR=sW+v*p-p*nwrap
                fl=2*sR-p*p+2*p-p*eps
                if fl<0 or fl%two_p: continue
                f=fl//two_p
                c_=hist[(p-1-v)%p]
                if f>c_: continue
                m=an[j,u]
                av[j,u,m]=v; af[j,u,m]=f; an[j,u]=m+1
                aull[j,u]=True
    return bases,av,af,an,aull

@njit(cache=True)
def _join_outer(p,k,q,bases,av,af,an,aull,UU,Tm,c0,eps,
                sol_buf,flip_sig,flip_f,counts,u_lo,u_hi):
    thi=(k-1)*eps+p; tlo=(k-1)*eps-p
    two_p=2*p
    nsol=counts[0]; nflip=counts[1]
    idx=np.zeros(k,np.int64)
    lens=np.zeros(k,np.int64)
    gbuf=np.zeros(q,np.int64)
    fv=np.zeros(k,np.int64)
    for ci in range(u_lo,u_hi):
        ok=True
        for j in range(k):
            if not aull[j,UU[ci,j]]: ok=False;break
        if not ok: continue
        for j in range(k):
            idx[j]=0
            lens[j]=an[j,UU[ci,j]]
        while True:
            vs=0; fs=0
            for j in range(k-1):
                vs+=av[j,UU[ci,j],idx[j]]
                fs+=af[j,UU[ci,j],idx[j]]
            vlast=(c0-vs)%p
            jl=k-1; ul=UU[ci,jl]
            pos=-1
            for t in range(lens[jl]):
                if av[jl,ul,t]==vlast: pos=t;break
            if pos>=0:
                ftot=fs+af[jl,ul,pos]
                if ftot==0:
                    good=True
                    for x in range(q):
                        ps=0
                        for j in range(k):
                            u=UU[ci,j]
                            v=av[j,u,idx[j]] if j<k-1 else vlast
                            b=bases[j,u,Tm[j,x]]+v
                            if b>=p: b-=p
                            ps+=2*b-p+2
                        if ps!=thi and ps!=tlo:
                            good=False;break
                    if good:
                        if nsol<sol_buf.shape[0]:
                            for x in range(q):
                                ps=0
                                for j in range(k):
                                    u=UU[ci,j]
                                    v=av[j,u,idx[j]] if j<k-1 else vlast
                                    b=bases[j,u,Tm[j,x]]+v
                                    if b>=p: b-=p
                                    ps+=2*b-p+2
                                sol_buf[nsol,x]=1 if ps==thi else -1
                        nsol+=1
                else:
                    # inline flip pre-validation: class + g bounds + line integrality
                    valid=True
                    for x in range(q):
                        ps=0
                        for j in range(k):
                            u=UU[ci,j]
                            v=av[j,u,idx[j]] if j<k-1 else vlast
                            b=bases[j,u,Tm[j,x]]+v
                            if b>=p: b-=p
                            ps+=2*b-p+2
                        d=ps-thi
                        if d%two_p!=0: valid=False;break
                        gx=d//two_p
                        if gx>k or gx+1<0: valid=False;break
                        gbuf[x]=gx
                    if valid:
                        F=0
                        for j in range(k):
                            fv[j]=af[j,UU[ci,j],idx[j]] if j<k-1 else af[jl,ul,pos]
                            F+=fv[j]
                        Gl=np.zeros((k,p),np.int64)
                        for x in range(q):
                            gx=gbuf[x]
                            for j in range(k):
                                Gl[j,Tm[j,x]]+=gx
                        for j in range(k):
                            if not valid: break
                            Fj=F-fv[j]
                            n1=0; n0=0
                            for s in range(p):
                                lo2=Gl[j,s]-Fj; hi2=Gl[j,s]+p-Fj
                                can0=(lo2<=0) and (0<=hi2)
                                can1=(lo2<=p) and (p<=hi2)
                                if not (can0 or can1):
                                    valid=False;break
                                if can1 and not can0: n1+=1
                                if can0 and not can1: n0+=1
                            if valid and (n1>fv[j] or fv[j]>p-n0):
                                valid=False
                    if valid:
                        if nflip<flip_sig.shape[0]:
                            for j in range(k):
                                u=UU[ci,j]
                                v=av[j,u,idx[j]] if j<k-1 else vlast
                                for s in range(p):
                                    b=bases[j,u,s]+v
                                    if b>=p: b-=p
                                    flip_sig[nflip,j,s]=2*b-p+2
                                flip_f[nflip,j]=fv[j]
                        nflip+=1
            c2=0
            while c2<k-1:
                idx[c2]+=1
                if idx[c2]<lens[c2]: break
                idx[c2]=0; c2+=1
            if c2>=k-1: break
    counts[0]=nsol; counts[1]=nflip

def process_outer(p,k,q,upper,UU,Tm,c0,eps,sol_buf,flip_sig,flip_f,sols,thi,tlo):
    """Run one outer completely, chunking uu so buffers can't overflow."""
    bases,av,af,an,aull=_prep_tables(p,k,upper,eps)
    worst=1
    for j in range(k-1): worst*=p
    UCH=max(1,sol_buf.shape[0]//max(1,worst))
    nu=UU.shape[0]
    for ulo in range(0,nu,UCH):
        counts=np.zeros(2,np.int64)
        _join_outer(p,k,q,bases,av,af,an,aull,UU,Tm,c0,eps,
                    sol_buf,flip_sig,flip_f,counts,ulo,min(ulo+UCH,nu))
        ns,nf=int(counts[0]),int(counts[1])
        if ns>sol_buf.shape[0] or nf>flip_sig.shape[0]:
            raise RuntimeError(f"overflow ns={ns} nf={nf}")
        for r in range(ns):
            sols.append(sol_buf[r].copy())
        if nf:
            S=flip_sig[:nf]
            PS=np.zeros((nf,q),dtype=np.int64)
            for j in range(k):
                PS+=S[:,j,:][:,Tm[j]]
            diff=PS-thi
            ok=((diff%(2*p))==0).all(axis=1)
            g=diff//(2*p)
            ok&=~((g>k).any(axis=1)|((g+1)<0).any(axis=1))
            for r in np.where(ok)[0]:
                sigstd=flip_sig[r].copy()
                fvec=[int(x) for x in flip_f[r]]
                idxs=[np.where(sigstd[j]==p)[0] for j in range(k)]
                _flip_solve(p,k,sigstd,idxs,fvec,Tm,thi,tlo,sols)

def _activity_filter(sols,Tm,p,k,eps):
    if not sols: return []
    Y=np.stack(sols).astype(np.int32)          # M x q
    keep=np.ones(len(Y),dtype=bool)
    for j in range(k):
        LS=np.zeros((len(Y),p),dtype=np.int32)
        for s in range(p):
            LS[:,s]=Y[:,Tm[j]==s].sum(axis=1)
        keep&=~(LS==eps).all(axis=1)
    A=np.stack(sols).astype(np.int8)
    return [A[i] for i in np.where(keep)[0]]

def enum_chunk(args):
    ctx,lo,hi,eps=args
    p=ctx['p'];q=ctx['q'];k=ctx['k'];deg=ctx['deg'];Tm=ctx['Tm']
    kern=ctx['kern'];UU=ctx['UU'];c0=ctx['c0'];outer_dims=ctx['outer_dims']
    s_ar=np.arange(p,dtype=np.int64)
    spow={d:(s_ar**d)%p for d in range(0,deg+1)}
    thi=(k-1)*eps+p; tlo=(k-1)*eps-p
    sols=[]
    SOLCAP=400000; FLIPCAP=400000
    sol_buf=np.zeros((SOLCAP,q),np.int8)
    flip_sig=np.zeros((FLIPCAP,k,p),np.int64)
    flip_f=np.zeros((FLIPCAP,k),np.int64)
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
        process_outer(p,k,q,upper,UU,Tm,c0,eps,sol_buf,flip_sig,flip_f,sols,thi,tlo)
    out=_activity_filter(sols,Tm,p,k,eps)
    return ctx['subset'],lo,hi,out
