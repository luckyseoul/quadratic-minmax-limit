import numpy as np,sys,time,itertools
sys.path.insert(0,'/tmp/e1work')
from kgen import square_coords
from kgen3 import prep_subset
from kgen5 import _activity_filter,process_outer as cpu_process_outer
from kgen6 import translation_tables
from dilation import build_group,orbits
from multiprocessing import Pool

p=11;k=5;q=p*p
dirs,forms,coords=square_coords(p)
m=len(dirs)
group,gcoords=build_group(p)
TT=translation_tables(p)
subsets=list(itertools.combinations(range(m),k))
ctxs={sub:prep_subset(p,list(sub),forms,coords) for sub in subsets}
s_ar=np.arange(p,dtype=np.int64)

def states_phaseT():
    out=[]
    for sub in subsets:
        K3=ctxs[sub]['kern'][3]
        assert len(K3)==1
        for lam in range(1,p):
            out.append((sub,{3:(lam*K3[0])%p,2:np.zeros(k,dtype=np.int64)}))
    return out

def states_phaseL():
    out=[]
    for sub in subsets:
        K2=ctxs[sub]['kern'][2]
        for e1 in range(p):
            for e2 in range(p):
                out.append((sub,{3:np.zeros(k,dtype=np.int64),
                                 2:(e1*K2[0]+e2*K2[1])%p}))
    return out

def enum_rep(args):
    sub,coeffs=args
    ctx=ctxs[sub]
    upper=np.zeros((k,p),dtype=np.int64)
    for d,vec in coeffs.items():
        upper=(upper+np.outer(vec,(s_ar**d)%p))%p
    sols=[]
    SOL=400000
    sol_buf=np.zeros((SOL,q),np.int8)
    fs=np.zeros((SOL,k,p),np.int64); ff=np.zeros((SOL,k),np.int64)
    thi=(k-1)+p; tlo=(k-1)-p
    cpu_process_outer(p,k,q,upper,ctx['UU'],ctx['Tm'],ctx['c0'],1,
                      sol_buf,fs,ff,sols,thi,tlo)
    return sub,coeffs,_activity_filter(sols,ctx['Tm'],p,k,1)

if __name__=='__main__':
    t0=time.time()
    total=0
    for phase,states in (('T',states_phaseT()),('L',states_phaseL())):
        orbs=orbits(states,group,p)
        print(f"phase {phase}: {len(states)} states -> {len(orbs)} orbits "
              f"(sizes {[len(tv) for _,tv in orbs]})",flush=True)
        with Pool(min(60,len(orbs))) as pool:
            results=pool.map(enum_rep,[st for st,tv in orbs])
        for (st,tv),(sub,cf,reps) in zip(orbs,results):
            mult=len(tv)*q if phase=='T' else len(tv)
            # phase L outers were never translation-gauged: solutions are full
            if phase=='T':
                total+=len(reps)*len(tv)*q
            else:
                total+=len(reps)*len(tv)
            print(f"   orbit |tv|={len(tv)}: reps={len(reps)} -> {len(reps)*mult if phase=='T' else len(reps)*len(tv)}",flush=True)
    print(f"k=5 gauged TOTAL = {total}  (expect 1306800)  time={time.time()-t0:.0f}s",flush=True)
