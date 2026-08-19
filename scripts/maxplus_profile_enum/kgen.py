"""General-k Max+ enumerator over polynomial profile coefficients.
A k-active (eps=+1) solution on direction subset J: profiles sigma_j, j in J,
sigma_j = 2*rhohat_j - p + 2 with rhohat_j an integer lift of a degree<=k-2
polynomial mod p, level-d coefficient vectors in ker(sum_j c_j t_j^d),
sum_s sigma_j(s) = p, and point sums in {(k-1) +- p}."""
import numpy as np,itertools
from math import comb

def field_ctx(p):
    q=p*p
    def is_irr(a,b): return all((x*x-a*x-b)%p!=0 for x in range(p))
    ia=ib=None
    for a in range(p):
        for b in range(p):
            if is_irr(a,b): ia,ib=a,b; break
        if ia is not None: break
    def mul(u,v):
        c0,c1=u%p,u//p; d0,d1=v%p,v//p
        return (c0*d0+c1*d1*ib)%p + ((c0*d1+c1*d0+c1*d1*ia)%p)*p
    def powm(u,e):
        r,base=1,u
        while e:
            if e&1: r=mul(r,base)
            base=mul(base,base); e>>=1
        return r
    def chi(x):
        return 0 if x==0 else (1 if powm(x,(q-1)//2)==1 else -1)
    def tr(x): return (2*(x%p)+ia*(x//p))%p
    return q,mul,chi,tr

def square_coords(p):
    q,mul,chi,tr=field_ctx(p)
    dirs=[];seen=set()
    for g in range(1,q):
        if g in seen: continue
        line=[mul(tt,g) for tt in range(1,p)]
        seen.update(line)
        if chi(g)==1: dirs.append(g)
    coords=[];forms=[]
    for g in dirs:
        cj=next(c for c in range(1,q) if tr(mul(c,g))==0)
        coords.append(np.array([tr(mul(c_,x)) for c_,x in ((cj,x) for x in range(q))],dtype=np.int64))
        forms.append((tr(mul(cj,1)),tr(mul(cj,p))))
    return dirs,forms,coords

def kernel_modp(A,p):
    A=np.array(A,dtype=np.int64)%p
    rows,cols=A.shape
    M=A.copy(); piv=[]; r=0
    for c in range(cols):
        pr=None
        for rr in range(r,rows):
            if M[rr,c]%p: pr=rr;break
        if pr is None: continue
        M[[r,pr]]=M[[pr,r]]
        M[r]=(M[r]*pow(int(M[r,c]),p-2,p))%p
        for rr in range(rows):
            if rr!=r and M[rr,c]:
                M[rr]=(M[rr]-M[rr,c]*M[r])%p
        piv.append(c); r+=1
    free=[c for c in range(cols) if c not in piv]
    basis=[]
    for fc in free:
        v=np.zeros(cols,dtype=np.int64); v[fc]=1
        for ri,c in enumerate(piv):
            v[c]=(-M[ri,fc])%p
        basis.append(v%p)
    return basis,piv

def enum_subset(p,subset,forms,coords,eps=1):
    """Enumerate all k-active eps=+1 solutions with support exactly `subset`.
    Returns list of int8 arrays y over the q finite points."""
    q=p*p;k=len(subset);deg=k-2
    F=[forms[j] for j in subset]; Tm=np.stack([coords[j] for j in subset],0)
    # kernels per degree
    kern={}
    for d in range(2,deg+1):
        A=[[comb(d,i)*pow(F[j][0],i)*pow(F[j][1],d-i) for j in range(k)] for i in range(d+1)]
        B,_=kernel_modp(A,p)
        kern[d]=B
    A1=[[F[j][0] for j in range(k)],[F[j][1] for j in range(k)]]
    U,_=kernel_modp(A1,p)
    assert len(U)==k-2
    # determined profiles for u-join: last two with invertible 2x2 of A1
    Adet=None
    detpair=None
    for pair in itertools.combinations(range(k),2):
        M=np.array([[A1[0][pair[0]],A1[0][pair[1]]],[A1[1][pair[0]],A1[1][pair[1]]]],dtype=np.int64)%p
        det=(M[0,0]*M[1,1]-M[0,1]*M[1,0])%p
        if det: detpair=pair; Minv=(pow(int(det),p-2,p)*np.array([[M[1,1],-M[0,1]],[-M[1,0],M[0,0]]],dtype=np.int64))%p; break
    freeidx=[j for j in range(k) if j not in detpair]
    c0=((k*p-k-1+p)//2)%p
    s_ar=np.arange(p,dtype=np.int64)
    spow={d:(s_ar**d)%p for d in range(0,deg+1)}
    # outer loop over kernel coords for degrees deg..2 (includes all-zero)
    outer_spaces=[]
    for d in range(deg,1,-1):
        dim=len(kern[d])
        outer_spaces.append((d,dim))
    def outer_iter():
        ranges=[range(p**dim) for d,dim in outer_spaces]
        for combo in itertools.product(*ranges):
            coeffs={}
            for (d,dim),code in zip(outer_spaces,combo):
                vec=np.zeros(k,dtype=np.int64)
                cc=code
                for b in range(dim):
                    vec=(vec+ (cc%p)*np.array(kern[d][b]))%p
                    cc//=p
                coeffs[d]=vec
            yield coeffs
    sols=[]
    two_p=2*p
    target_hi=(k-1)*eps+p; target_lo=(k-1)*eps-p
    for coeffs in outer_iter():
        # upper table per profile
        upper=np.zeros((k,p),dtype=np.int64)
        for d in range(2,deg+1):
            upper=(upper+np.outer(coeffs[d],spow[d]))%p
        # per profile allowed (u,v): pass row filter
        allowed=[dict() for _ in range(k)]  # u -> list of (v, f, cnt_idx info)
        for j in range(k):
            for u in range(p):
                W=(upper[j]+u*s_ar)%p
                hist=np.bincount(W,minlength=p)
                sW=int(W.sum())
                # cumulative: #\{s: W[s] >= p-v\}
                cs=np.cumsum(hist[::-1])  # cs[t-1]=#\{W>=p-t\}
                lst=[]
                for v in range(p):
                    nwrap=int(cs[v-1]) if v>0 else 0
                    sR=sW+v*p-p*nwrap
                    rows=2*sR-p*p+2*p
                    fl=rows-p*eps
                    if fl<0 or fl%two_p: continue
                    f=fl//two_p
                    cnt=int(hist[(p-1-v)%p])
                    if f>cnt: continue
                    lst.append((v,f))
                if lst: allowed[j][u]=lst
            if not allowed[j]: break
        else:
            # u-join: iterate free profiles' u from allowed keys
            freelists=[sorted(allowed[j].keys()) for j in freeidx]
            for ucombo in itertools.product(*freelists):
                # determined u's: solve A1 . u = 0 for detpair
                rhs0=rhs1=0
                for j,uv in zip(freeidx,ucombo):
                    rhs0=(rhs0+A1[0][j]*uv)%p; rhs1=(rhs1+A1[1][j]*uv)%p
                ud0=int((-(Minv[0,0]*rhs0+Minv[0,1]*rhs1))%p)
                ud1=int((-(Minv[1,0]*rhs0+Minv[1,1]*rhs1))%p)
                if ud0 not in allowed[detpair[0]] or ud1 not in allowed[detpair[1]]: continue
                uu=np.zeros(k,dtype=np.int64)
                for j,uv in zip(freeidx,ucombo): uu[j]=uv
                uu[detpair[0]]=ud0; uu[detpair[1]]=ud1
                # v-join: product of v-lists with sum == c0 mod p
                vlists=[allowed[j][int(uu[j])] for j in range(k)]
                for vcombo in itertools.product(*vlists[:-1]):
                    vs=sum(v for v,f in vcombo)%p
                    vlast=(c0-vs)%p
                    hit=[vf for vf in vlists[-1] if vf[0]==vlast]
                    if not hit: continue
                    full=list(vcombo)+[hit[0]]
                    # build sigma tables + flips branching
                    R=np.stack([(upper[j]+uu[j]*s_ar+full[j][0])%p for j in range(k)],0)
                    sig=2*R-p+2
                    fvec=[full[j][1] for j in range(k)]
                    idxs=[np.where(R[j]==p-1)[0] for j in range(k)]
                    choices=[list(itertools.combinations(idxs[j],fvec[j])) for j in range(k)]
                    for pick in itertools.product(*choices):
                        sg=sig.copy()
                        for j in range(k):
                            for s_ in pick[j]: sg[j,s_]-=two_p
                        # activity: sigma_j not identically eps
                        if any((sg[j]==eps).all() for j in range(k)): continue
                        PS=sg[0][Tm[0]]
                        for j in range(1,k): PS=PS+sg[j][Tm[j]]
                        mn=PS.min(); mx=PS.max()
                        if (mn==target_lo or mn==target_hi) and (mx==target_lo or mx==target_hi) and \
                           np.isin(PS,(target_lo,target_hi)).all():
                            y=np.where(PS==target_hi,1,-1).astype(np.int8)
                            sols.append(y)
    return sols
