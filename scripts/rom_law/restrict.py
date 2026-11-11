import numpy as np, sys, subprocess, time
from paley import paley
from lifts import write
def phi(A,tag):
    n=A.shape[0]; write(A,tag+'.bin')
    if n<=30:
        o=subprocess.check_output(['./gs2',tag+'.bin']).decode().split(); return max(float(o[0]),-float(o[1]))/2
    ps=[subprocess.Popen(['./gs2',tag+'.bin',str(k)],stdout=subprocess.PIPE) for k in (0,1)]
    outs=[p.communicate()[0].decode().split() for p in ps]
    return max(max(float(o[0]) for o in outs),-min(float(o[1]) for o in outs))/2
if __name__=='__main__':
    q=int(sys.argv[1]); dels=[int(s) for s in sys.argv[2:]]
    C=paley(q); n=q+1
    for d in dels:
        keep=list(range(d,n))   # delete the first d points (includes the point at infinity)
        t=time.time(); A=C[np.ix_(keep,keep)]; v=phi(A,f'r{q}_{d}')
        print(f"C{n} minus {d} points -> order {n-d}: Phi = {v:.0f}   ({time.time()-t:.0f}s)",flush=True)
