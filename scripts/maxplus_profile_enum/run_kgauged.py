"""Generic dilation+translation-gauged GPU enumeration of the k-stratum.

Usage: python3 run_kgauged.py <k> <workers> [validate]
Phase T: top = lam*K[deg][0] (lam != 0), level deg-1 = 0 (translation gauge),
         lower levels full lattices; expand reps by transversal x translations.
Phase L: top = 0, levels deg-1..2 full lattices; expand by transversal only.
Orbits under the 120-element dilation/Frobenius group across subsets.
"""
import numpy as np, sys, time, itertools, os, pickle
sys.path.insert(0,'/tmp/e1work')
os.environ.setdefault("OMP_NUM_THREADS","1")
from kgen import square_coords
from kgen3 import prep_subset
from kgen5 import _activity_filter
from kgen6 import translation_tables
from dilation import build_group, orbits
from multiprocessing import Pool

p=11; q=p*p

def lattice(basis):
    if len(basis)==0:
        return [np.zeros(0,dtype=np.int64)]
    basis=np.array(basis); dim=len(basis)
    out=[]
    for combo in itertools.product(range(p),repeat=dim):
        v=np.zeros(basis.shape[1],dtype=np.int64)
        for c,b in zip(combo,basis): v=(v+c*b)%p
        out.append(v)
    return out

def build_states(k,ctxs,subsets):
    deg=k-2
    T=[];L=[]
    for sub in subsets:
        kern=ctxs[sub]['kern']
        low_levels=[d for d in range(deg-2,1,-1)]
        lows=[lattice(kern[d]) for d in low_levels]
        for lam in range(1,p):
            for combo in itertools.product(*lows):
                cf={deg:(lam*kern[deg][0])%p}
                if deg-1>=2: cf[deg-1]=np.zeros(k,dtype=np.int64)
                for d,vec in zip(low_levels,combo): cf[d]=vec
                T.append((sub,cf))
        levs=[d for d in range(deg-1,1,-1)]
        lats=[lattice(kern[d]) for d in levs]
        for combo in itertools.product(*lats):
            cf={deg:np.zeros(k,dtype=np.int64)}
            for d,vec in zip(levs,combo): cf[d]=vec
            L.append((sub,cf))
    return T,L

def worker(args):
    phase,sub,cf,tvidx=args
    global _ctxs,_testers,_k
    from gpu_inner import GpuTester,process_outer_gpu
    ctx=_ctxs[sub]
    if sub not in _testers:
        _testers[sub]=GpuTester(p,_k,ctx['Tm'],ctx['UU'])
    tester=_testers[sub]
    s_ar=np.arange(p,dtype=np.int64)
    upper=np.zeros((_k,p),dtype=np.int64)
    for d,vec in cf.items():
        upper=(upper+np.outer(vec,(s_ar**d)%p))%p
    sols=[]
    t0=time.time()
    nc=process_outer_gpu(p,_k,q,upper,ctx['UU'],ctx['Tm'],ctx['c0'],1,tester,sols)
    reps=_activity_filter(sols,ctx['Tm'],p,_k,1)
    R=np.stack(reps).astype(np.int8) if reps else np.zeros((0,q),np.int8)
    return phase,sub,tvidx,R,nc,time.time()-t0

def initw(k,ctxs):
    global _ctxs,_testers,_k
    _ctxs=ctxs;_testers={};_k=k

if __name__=='__main__':
    k=int(sys.argv[1]); nw=int(sys.argv[2])
    validate=len(sys.argv)>3 and sys.argv[3]=='validate'
    dirs,forms,coords=square_coords(p)
    m=len(dirs)
    subsets=list(itertools.combinations(range(m),k))
    ctxs={sub:prep_subset(p,list(sub),forms,coords) for sub in subsets}
    group,gcoords=build_group(p)
    TT=translation_tables(p)
    t0=time.time()
    T,L=build_states(k,ctxs,subsets)
    print(f"states: T={len(T)} L={len(L)}",flush=True)
    orbT=orbits(T,group,p); orbL=orbits(L,group,p)
    print(f"orbits: T={len(orbT)} L={len(orbL)}  ({time.time()-t0:.0f}s)",flush=True)
    outdir=f'/tmp/e1work/k{k}_gpu_out'
    os.makedirs(outdir,exist_ok=True)
    tasks_all=[]
    tvstore=[]
    for phase,orb in (('T',orbT),('L',orbL)):
        for (sub,cf),tv in orb:
            tvstore.append(tv)
            tasks_all.append((phase,sub,cf,len(tvstore)-1))
    # RESUME: tvidx assignment is deterministic (orbits() iterates a
    # dict built by enumerate() over a deterministically-ordered states
    # list), so orb{tvidx}.npy from a prior interrupted run lines up with
    # this run's task list exactly. Skip anything already on disk.
    total=0
    done_prior=0
    tasks=[]
    for t in tasks_all:
        tvidx=t[3]
        fp=f'{outdir}/orb{tvidx}.npy'
        if os.path.exists(fp):
            total+=len(np.load(fp))
            done_prior+=1
        else:
            tasks.append(t)
    if done_prior:
        print(f"RESUMED: {done_prior}/{len(tasks_all)} outers already on disk "
              f"({total} solutions); {len(tasks)} remaining",flush=True)
    perm_inv={id(g):np.argsort(g[0]) for g in group}
    done=done_prior
    with Pool(nw,initializer=initw,initargs=(k,ctxs)) as pool:
        for phase,sub,tvidx,R,nc,dt in pool.imap_unordered(worker,tasks):
            tv=tvstore[tvidx]
            cnt_here=0
            outs=[]
            if len(R):
                for g in tv:
                    pi=np.argsort(g[0])
                    Rg=R[:,pi]
                    if phase=='T':
                        for y in Rg:
                            outs.append(y[TT])
                    else:
                        outs.append(Rg)
            if outs:
                A=np.concatenate([o if o.ndim==2 else o[None,:] for o in outs])
                cnt_here=len(A)
                fp=f'{outdir}/orb{tvidx}.npy'
                np.save(fp,A.astype(np.int8))
                with open(fp,'rb') as fh:
                    os.fsync(fh.fileno())
            else:
                # zero-solution orbit: still mark it done so a restart skips it
                np.save(f'{outdir}/orb{tvidx}.npy',np.zeros((0,q),np.int8))
            total+=cnt_here
            done+=1
            print(f"[{done}/{len(tasks_all)}] {phase} reps={len(R)} |tv|={len(tv)} -> {cnt_here}  cand={nc} {dt:.0f}s  cum={total}",flush=True)
    print(f"k={k} gauged GPU TOTAL = {total}",flush=True)
    if validate and k==5:
        print("comparing against k5_p11_full.npy as a SET ...",flush=True)
        import glob
        parts=[np.load(f) for f in glob.glob(f'{outdir}/orb*.npy')]
        A=np.concatenate(parts)
        S1=set(map(tuple,A.tolist()))
        B=np.load('/tmp/e1work/k5_p11_full.npy')
        S2=set(map(tuple,B.tolist()))
        print(f"sets equal: {S1==S2}  |S1|={len(S1)} |S2|={len(S2)}",flush=True)
