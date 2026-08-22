#!/usr/bin/env python3
"""Per-VECTOR constituent energies ||P_c B_y||^2 for y in Max+.
Fair share = n(n-2)/dimZ * mult_c = lbar*mult_c. Need >= 6*mult_c.
Diagnostic: is the spread pointwise-good, or only good on average?"""
import numpy as np, itertools, sys
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
from e1_gmin_m4_prop15590 import MuLab
p=int(sys.argv[1])
lab=MuLab(p,with_deg6=False); n=lab.n; C=lab.C.astype(np.float64)
Y=lab.Yp.astype(np.float64); Pp=(np.eye(n)+C/p)/2
lbar=8*(n-2)/(n-6)
# orthonormal basis of Z
pairs=list(itertools.combinations(range(n),2)); A=[]
for (i,j) in pairs:
    E=np.zeros((n,n)); E[i,j]=E[j,i]=1
    A.append((C@E-p*E).reshape(-1))
Am=np.array(A).T; s=np.linalg.svd(Am,compute_uv=False); _,_,Vt=np.linalg.svd(Am)
Q=[]
for v in Vt[np.sum(s>1e-8):]:
    B=np.zeros((n,n))
    for idx,(i,j) in enumerate(pairs): B[i,j]=B[j,i]=v[idx]
    for W in Q: B=B-(B*W).sum()*W
    nr=np.linalg.norm(B)
    if nr>1e-8: Q.append(B/nr)
d=len(Q); Qa=np.array(Q)
# Phi and its eigen-decomposition -> constituent projectors via eigenspaces
qf=np.empty((len(Y),d))
for c in range(d): qf[:,c]=((Y@Qa[c])*Y).sum(axis=1)
Phi=(qf.T@qf)/len(Y)
ev,U=np.linalg.eigh(Phi)
# cluster eigenvalues into constituents
cl=[]
for idx in np.argsort(ev):
    e=ev[idx]
    if cl and abs(e-cl[-1][0])<1e-6: cl[-1][1].append(idx)
    else: cl.append([e,[idx]])
print(f"p={p} n={n} dimZ={d} lbar={lbar:.4f}  |Max+|={len(Y)}")
print(f"constituents: {[(round(a,4),len(ix)) for a,ix in cl]}")
# per-vector energy in each constituent: ||P_c B_y||^2 = sum_{k in c} <B_y,u_k>^2
coef=qf@U                                       # y x d in eigenbasis
print(f"\nper-vector ||P_c B_y||^2 (fair share = lbar*mult; need >= 6*mult):")
print(f"{'lambda_c':>10} {'mult':>5} {'need6':>9} {'fair':>9} {'min_y':>10} {'mean_y':>10} {'max_y':>10}  pointwise>=6?")
allok=True
for a,ix in cl:
    m=len(ix)
    Ec=(coef[:,ix]**2).sum(axis=1)
    ok=Ec.min()>=6*m-1e-9
    allok&=ok
    print(f"{a:10.4f} {m:5d} {6*m:9.1f} {lbar*m:9.2f} {Ec.min():10.3f} {Ec.mean():10.3f} {Ec.max():10.3f}  {ok}")
print(f"\nPOINTWISE bound holds for every y and every constituent: {allok}")
tot=(coef**2).sum(axis=1)
print(f"total energy per y: min={tot.min():.3f} max={tot.max():.3f} (exact n(n-2)={n*(n-2)})")
