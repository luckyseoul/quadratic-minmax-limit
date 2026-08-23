#!/usr/bin/env python3
import sys,json,time; sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from concurrent.futures import ProcessPoolExecutor
p=int(sys.argv[1])
def reps():
    out=[]
    for a in range(p):
        for b in range(p):
            for c in range(p):
                if (a*a+b*c)%p==1%p:
                    t=(a,b,c); tn=((-a)%p,(-b)%p,(-c)%p)
                    if t<=tn: out.append((a,b,c,(-a)%p))
    return out
def work(t):
    from e1_gmin_m4_prop15626 import _switched
    A,B,C,D=t
    r=_switched(p,A,B,C,D)
    return (r["eigen_minus"],r["inU_y"],r["W2"])
if __name__=="__main__":
    R=reps(); print(f"p={p} class size {len(R)} (expect {p*(p+1)//2})",flush=True)
    t0=time.time(); res=[]
    with ProcessPoolExecutor(max_workers=60) as ex:
        for i,r in enumerate(ex.map(work,R,chunksize=4)):
            res.append(r)
            if i%200==199: print(f"  {i+1}/{len(R)} {time.time()-t0:.0f}s",flush=True)
    ne=sum(1 for e,u,w in res if e); nu=sum(1 for e,u,w in res if e and u)
    nw=sum(1 for e,u,w in res if e and u and w)
    print(f"p={p}: class={len(R)} eigen={ne} inU={nu} W2={nw}  rate_inU={nw/max(nu,1):.4f}",flush=True)
    json.dump({"p":p,"n_class":len(R),"n_eigen":ne,"n_inU":nu,"n_W2":nw},
              open(f"/tmp/w2_p{p}.json","w"))
