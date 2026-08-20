import numpy as np,itertools,sys,time,collections
sys.path.insert(0,'.')
from kgen import square_coords
from kgen3 import prep_subset
from kgen5 import _activity_filter
from gpu_inner import GpuTester, process_outer_gpu

def scan(p):
    dirs,forms,coords=square_coords(p)
    m=len(dirs); q=p*p
    s_ar=np.arange(p,dtype=np.int64)
    per={}
    t0=time.time()
    for sub in itertools.combinations(range(m),4):
        ctx=prep_subset(p,list(sub),forms,coords)
        tester=GpuTester(p,4,ctx['Tm'],ctx['UU'])
        sols=[]
        for lam in range(p):
            upper=(lam*np.outer(ctx['kern'][2][0],(s_ar**2)%p))%p
            process_outer_gpu(p,4,q,upper,ctx['UU'],ctx['Tm'],ctx['c0'],1,tester,sols)
        out=_activity_filter(sols,ctx['Tm'],p,4,1)
        per[sub]=len(out)
    dist=collections.Counter(v//q if v%q==0 else -v for v in per.values())
    tot=sum(per.values())
    print(f"p={p}: k4 total={tot} = {tot/q:.3f} q ; per-subset counts/q: {dict(dist)} ; time {time.time()-t0:.0f}s",flush=True)
    return per

if __name__=='__main__':
    primes=[int(x) for x in sys.argv[1:]] or [29,31,37,41]
    for p in primes:
        scan(p)
