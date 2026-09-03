#!/usr/bin/env python3
import argparse, time
import numpy as np
from ortools.sat.python import cp_model

ap = argparse.ArgumentParser()
ap.add_argument('--directions', default='v,0,1,2,3,4')
ap.add_argument('--m-start', type=int, default=4)
ap.add_argument('--m-end', type=int, default=60)
ap.add_argument('--seconds', type=float, default=3.0)
ap.add_argument('--tag', default='oa')
args = ap.parse_args()
q = 11
D = {None if x == 'v' else int(x) for x in args.directions.split(',')}
assert len(D) == 6
pts = [(a,b) for a in range(q) for b in range(q)]
A = np.zeros((q*q,q*q), dtype=np.int16)
for i,(a,b) in enumerate(pts):
    for j,(c,d) in enumerate(pts):
        if i == j: continue
        da=(c-a)%q; db=(d-b)%q
        slope = None if da == 0 else (db*pow(da,-1,q))%q
        A[i,j] = int(slope in D)
assert np.all(A.sum(1)==60)
assert np.array_equal(A@A,30*np.ones((121,121),dtype=np.int16)+30*np.eye(121,dtype=np.int16)-A)

def eig(B, sign, secs=30.0):
    n=len(B); mdl=cp_model.CpModel(); z=[mdl.new_bool_var(f"z{i}") for i in range(n)]
    if sign==1:
        mdl.add(sum(z)==55)
        for i in range(n): mdl.add(sum(int(B[i,j])*z[j] for j in range(n))+6*z[i]==30)
    else:
        mdl.add(sum(z)==66)
        for i in range(n): mdl.add(sum(int(B[i,j])*z[j] for j in range(n))-5*z[i]==30)
    sol=cp_model.CpSolver(); sol.parameters.max_time_in_seconds=secs; sol.parameters.num_workers=16
    st=sol.solve(mdl)
    return sol.status_name(st), ([i for i,v in enumerate(z) if sol.value(v)] if st in (cp_model.FEASIBLE,cp_model.OPTIMAL) else None)

def switched(B, W):
    C=B.copy(); W=set(W); m=len(W)
    for v in range(len(B)):
        if v in W: continue
        if sum(int(B[v,w]) for w in W)*2==m:
            for w in W: C[v,w]=C[w,v]=1-C[v,w]
    return C

def find_sets(B,m,r,limit=20,secs=10.0):
    n=len(B); mdl=cp_model.CpModel(); w=[mdl.new_bool_var(f"w{i}") for i in range(n)]
    mdl.add(sum(w)==m); mdl.add(w[0]==1)
    ds=[]
    for i in range(n):
        di=mdl.new_int_var(0,m,f"d{i}"); ds.append(di)
        mdl.add(di==sum(int(B[i,j])*w[j] for j in range(n)))
        mdl.add_allowed_assignments([w[i],di],[(1,r),(0,0),(0,m//2),(0,m)])
    out=[]; deadline=time.time()+secs
    while len(out)<limit and time.time()<deadline:
        sol=cp_model.CpSolver(); sol.parameters.max_time_in_seconds=max(.1,deadline-time.time()); sol.parameters.num_workers=16
        st=sol.solve(mdl)
        if st not in (cp_model.FEASIBLE,cp_model.OPTIMAL): break
        W=[i for i,x in enumerate(w) if sol.value(x)]; out.append(W)
        mdl.add(sum(w[i] for i in W)<=m-1)
    return out

seen=set()
for m in range(args.m_start, args.m_end + 1, 2):
  for r in range(m):
    sets=find_sets(A,m,r,secs=args.seconds)
    if sets: print("SETS",m,r,len(sets),flush=True)
    for W in sets:
      key=tuple(W)
      if key in seen: continue
      seen.add(key); B=switched(A,W)
      if np.array_equal(B,A) or not np.all(B.sum(1)==60): continue
      if not np.array_equal(B@B,30*np.ones((121,121),dtype=np.int16)+30*np.eye(121,dtype=np.int16)-B): continue
      ep=eig(B,1); em=eig(B,-1)
      print("CAND",m,r,W,"eig+",ep[0],"eig-",em[0],flush=True)
      if ep[0]=="INFEASIBLE" and em[0]=="INFEASIBLE":
        np.savetxt(f'/tmp/{args.tag}_nonregular.csv',B,fmt='%d',delimiter=',')
        print("FOUND",m,r,W,flush=True); raise SystemExit(10)
print("DONE",len(seen))

