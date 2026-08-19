import numpy as np,itertools,sys,time,os
sys.path.insert(0,'/tmp/e1work')
os.environ.setdefault("OMP_NUM_THREADS","1")
from kgen import square_coords
from kgen4 import enum_chunk
from kgen3 import prep_subset
from multiprocessing import Pool

def run(p,k,workers=84,chunk=None,label=""):
    dirs,forms,coords=square_coords(p)
    m=len(dirs)
    tasks=[]
    for sub in itertools.combinations(range(m),k):
        ctx=prep_subset(p,list(sub),forms,coords)
        tot=ctx['outer_total']
        step=chunk or max(1,tot//(workers*4))
        for lo in range(0,tot,step):
            tasks.append((ctx,lo,min(lo+step,tot),1))
    t0=time.time()
    per_sub={}
    with Pool(workers) as pool:
        for sub,lo,hi,sols in pool.imap_unordered(enum_chunk,tasks):
            per_sub.setdefault(sub,[]).extend(sols)
    allsols=[]
    for sub in sorted(per_sub):
        print(f"  p={p} k={k} subset {sub}: {len(per_sub[sub])}",flush=True)
        allsols.extend(per_sub[sub])
    print(f"p={p} k={k}{label}: TOTAL={len(allsols)}  time={time.time()-t0:.1f}s",flush=True)
    if allsols:
        A=np.stack(allsols)
        d=len(set(map(tuple,A.tolist())))
        print(f"   distinct={d}",flush=True)
        np.save(f'/tmp/e1work/k{k}_p{p}_full.npy',A)
    return allsols

if __name__=='__main__':
    which=sys.argv[1]
    if which=='validate':
        s=run(7,4)          # expect 4410
        s=run(11,4)         # expect >= 58080 (lam!=0 gave 58080)
    elif which=='k5':
        run(11,5)
    elif which=='k6':
        print('k6 skipped (use run_gauged.py)')
