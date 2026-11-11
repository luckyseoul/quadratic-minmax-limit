import numpy as np, struct, subprocess, sys, json, time, os
from multiprocessing import Pool
def haar(n,rng):
    G=rng.standard_normal((n,n)); Q,R=np.linalg.qr(G); return Q*np.sign(np.diag(R))
def one(args):
    n,seed=args
    rng=np.random.default_rng(seed)
    O=haar(n,rng); d=np.array([1.0]*(n//2)+[-1.0]*(n-n//2))
    J=(O*d)@O.T
    name=f'tmp_{n}_{seed}.bin'
    with open(name,'wb') as f:
        f.write(struct.pack('i',n)); f.write(np.ascontiguousarray(J).tobytes())
    out=subprocess.check_output(['./gs2',name]).decode().split(); os.remove(name)
    return float(out[0])/n, -float(out[1])/n
if __name__=='__main__':
    n=int(sys.argv[1]); k=int(sys.argv[2])
    t=time.time()
    with Pool(2) as p: res=p.map(one,[(n,1000*n+s) for s in range(k)])
    a=np.array(res); one_sided=a.ravel(); two=a.max(axis=1)
    r=dict(n=n,instances=k,one_sided_mean=float(one_sided.mean()),one_sided_sd=float(one_sided.std()),
           two_sided_mean=float(two.mean()),two_sided_min=float(two.min()),one_sided_min=float(one_sided.min()),secs=round(time.time()-t,1))
    print(json.dumps(r)); open('rom_results.jsonl','a').write(json.dumps(r)+'\n')
