import numpy as np, struct, subprocess, sys, time
from paley import paley
from lifts import write, lift
C=paley(17); d=int(sys.argv[1]); A=lift(C,d); n=A.shape[0]
assert n==36 and np.allclose(A,A.T) and np.all(np.diag(A)==0)
write(A,f'lift36_{d}.bin'); t=time.time()
ps=[subprocess.Popen(['./gs2',f'lift36_{d}.bin',str(k)],stdout=subprocess.PIPE) for k in (0,1)]
outs=[p.communicate()[0].decode().split() for p in ps]
mx=max(float(o[0]) for o in outs); mn=min(float(o[1]) for o in outs)
print(f"order 36 from C18, D={d:+d}I: max Q={mx/2:.0f} min Q={mn/2:.0f} Phi={max(mx,-mn)/2:.0f}  ({time.time()-t:.0f}s)  S=108, ROM*S~103.4")
