import numpy as np,itertools,sys,time,os
sys.path.insert(0,'/tmp/e1work')
os.environ.setdefault("OMP_NUM_THREADS","1")
from kgen import square_coords
from kgen3 import prep_subset
from kgen6 import enum_gauged_task,translation_tables
from multiprocessing import Pool

def run(p,k,workers=84):
    dirs,forms,coords=square_coords(p)
    m=len(dirs)
    TT=translation_tables(p)
    tasks=[]
    for sub in itertools.combinations(range(m),k):
        ctx=prep_subset(p,list(sub),forms,coords)
        deg=k-2
        # phase T size: (p-1) * prod lattices deg-2..2
        tsize=(p-1)
        for d in range(deg-2,1,-1): tsize*=p**len(ctx['kern'][d])
        step=max(1,tsize//(workers))
        for lo in range(0,tsize,step):
            tasks.append((ctx,'T',lo,min(lo+step,tsize),1,TT))
        # phase L size
        lsize=1
        for d in range(deg-1,1,-1): lsize*=p**len(ctx['kern'][d])
        step=max(1,lsize//(workers))
        for lo in range(0,lsize,step):
            tasks.append((ctx,'L',lo,min(lo+step,lsize),1,TT))
    t0=time.time()
    per={}
    reps={}
    with Pool(workers) as pool:
        for sub,phase,lo,hi,sols,nre in pool.imap_unordered(enum_gauged_task,tasks):
            per.setdefault(sub,[]).extend(sols)
            reps[sub]=reps.get(sub,0)+(nre if phase=='T' else 0)
    allsols=[]
    for sub in sorted(per):
        print(f"  p={p} k={k} subset {sub}: {len(per[sub])} (T-reps {reps.get(sub,0)})",flush=True)
        allsols.extend(per[sub])
    print(f"p={p} k={k} gauged: TOTAL={len(allsols)} time={time.time()-t0:.1f}s",flush=True)
    if allsols:
        A=np.stack(allsols)
        ds=len(set(map(tuple,A.tolist())))
        print(f"   distinct={ds}",flush=True)
        np.save(f'/tmp/e1work/k{k}_p{p}_full.npy',A)
    return allsols

if __name__=='__main__':
    which=sys.argv[1]
    if which=='k4v':
        run(7,4); run(11,4)
    elif which=='k5':
        run(11,5)
    elif which=='k6':
        run(11,6)
