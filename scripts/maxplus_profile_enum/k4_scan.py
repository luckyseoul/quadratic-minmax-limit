import numpy as np,itertools,sys,time,collections
sys.path.insert(0,'/tmp/e1work')
from kgen import square_coords
from kgen4 import enum_chunk
from kgen3 import prep_subset
from multiprocessing import Pool

def scan(p,workers=8):
    dirs,forms,coords=square_coords(p)
    m=len(dirs)
    tasks=[]
    for sub in itertools.combinations(range(m),4):
        ctx=prep_subset(p,list(sub),forms,coords)
        tasks.append((ctx,0,ctx['outer_total'],1))
    t0=time.time()
    per={}
    with Pool(workers) as pool:
        for sub,lo,hi,sols in pool.imap_unordered(enum_chunk,tasks):
            per[sub]=len(sols)
    q=p*p
    dist=collections.Counter(v//q if v%q==0 else -v for v in per.values())
    tot=sum(per.values())
    print(f"p={p}: k4 total={tot} = {tot/q:.3f} q ; per-subset counts/q: {dict(dist)} ; time {time.time()-t0:.0f}s",flush=True)
    return per

if __name__=='__main__':
    for p in (17,19,23):
        scan(p)
