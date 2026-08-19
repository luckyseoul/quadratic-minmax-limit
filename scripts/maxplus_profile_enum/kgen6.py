"""Gauged general-k enumerator (top stratum: next-to-top level = 0, expand by
translations). Phases:
  T: top coeff = lam*w0 (lam != 0), level deg-1 = 0, lower levels full.
  L: top = 0, recurse structure (here: full enumeration of remaining levels).
Counts are exact after q-translation expansion of phase-T reps (free action,
one rep per orbit; asserted via divisibility cross-checks by callers)."""
import numpy as np,itertools
from math import comb
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen import field_ctx
from kgen3 import prep_subset
from kgen5 import process_outer
from kgen4 import _flip_solve

def _run_outers(ctx,outer_iter,eps=1):
    """outer_iter yields dicts {level d: coeff vec}. Returns solutions."""
    p=ctx['p'];q=ctx['q'];k=ctx['k'];deg=ctx['deg'];Tm=ctx['Tm']
    UU=ctx['UU'];c0=ctx['c0']
    s_ar=np.arange(p,dtype=np.int64)
    spow={d:(s_ar**d)%p for d in range(0,deg+1)}
    thi=(k-1)*eps+p; tlo=(k-1)*eps-p
    SOLCAP=400000; FLIPCAP=400000
    sol_buf=np.zeros((SOLCAP,q),np.int8)
    flip_sig=np.zeros((FLIPCAP,k,p),np.int64)
    flip_f=np.zeros((FLIPCAP,k),np.int64)
    sols=[]
    for coeffs in outer_iter:
        upper=np.zeros((k,p),dtype=np.int64)
        for d,vec in coeffs.items():
            upper=(upper+np.outer(vec,spow[d]))%p
        process_outer(p,k,q,upper,UU,Tm,c0,eps,sol_buf,flip_sig,flip_f,sols,thi,tlo)
    from kgen5 import _activity_filter
    return _activity_filter(sols,Tm,p,k,eps)

def _lattice(basis,p):
    if len(basis)==0: return [np.zeros(basis.shape[1] if hasattr(basis,'shape') else 0,dtype=np.int64)]
    basis=np.array(basis)
    dim=len(basis)
    out=[]
    for combo in itertools.product(range(p),repeat=dim):
        v=np.zeros(basis.shape[1],dtype=np.int64)
        for c,b in zip(combo,basis): v=(v+c*b)%p
        out.append(v)
    return out

def enum_gauged_task(args):
    """One task: (ctx, phase, chunkspec). Phase 'T': lam-chunk with level deg-1 = 0.
    Phase 'L': top=0 full enumeration chunk over levels deg-1..2."""
    ctx,phase,lo,hi,eps,translations=args
    p=ctx['p'];k=ctx['k'];deg=ctx['deg'];kern=ctx['kern'];q=ctx['q']
    if phase=='T':
        assert len(kern[deg])==1
        w0=kern[deg][0]
        assert (w0%p!=0).all(), "top kernel vector must have full support"
        # levels deg-1 = 0; levels deg-2..2 full lattice
        low_levels=[d for d in range(deg-2,1,-1)]
        lows=[_lattice(kern[d],p) for d in low_levels]
        def it():
            cnt=0
            for lam in range(1,p):
                for combo in itertools.product(*lows):
                    if not (lo<=cnt<hi): cnt+=1; continue
                    cnt+=1
                    coeffs={deg:(lam*w0)%p}
                    if deg-1>=2: coeffs[deg-1]=np.zeros(k,dtype=np.int64)
                    for d,vec in zip(low_levels,combo): coeffs[d]=vec
                    yield coeffs
        reps=_run_outers(ctx,it(),eps)
        out=[]
        for y in reps:
            out.extend(np.asarray(y)[translations].astype(np.int8))
        return ctx['subset'],phase,lo,hi,out,len(reps)
    else:
        # top = 0: full lattices for levels deg-1 .. 2
        low_levels=[d for d in range(deg-1,1,-1)]
        lows=[_lattice(kern[d],p) for d in low_levels]
        def it():
            cnt=0
            for combo in itertools.product(*lows):
                if not (lo<=cnt<hi): cnt+=1; continue
                cnt+=1
                coeffs={deg:np.zeros(k,dtype=np.int64)}
                for d,vec in zip(low_levels,combo): coeffs[d]=vec
                yield coeffs
        sols=_run_outers(ctx,it(),eps)
        return ctx['subset'],phase,lo,hi,sols,len(sols)

def translation_tables(p):
    q=p*p
    _,mul,chi,tr=field_ctx(p)
    def add(u,v): return (u%p+v%p)%p + (((u//p)+(v//p))%p)*p
    T=np.zeros((q,q),dtype=np.int64)
    for s in range(q):
        for x in range(q):
            T[s,x]=add(x,s)
    return T
