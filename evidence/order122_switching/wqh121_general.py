#!/usr/bin/env python3
"""Exact WQH partition finder for small equal cells in a 121 conference graph."""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0,str(Path(__file__).resolve().parent))
import gm121_scan as gs
import peisert121_exact as pe

def run(A,ell,kin,cross,seconds,workers,seed):
 n=len(A); m=cp_model.CpModel(); a=[m.NewBoolVar(f'a{i}') for i in range(n)]; b=[m.NewBoolVar(f'b{i}') for i in range(n)]
 m.Add(sum(a)==ell); m.Add(sum(b)==ell); m.Add(a[0]==1)
 for i in range(n):
  m.Add(a[i]+b[i]<=1)
  da=sum(int(A[i,j])*a[j] for j in range(n)); db=sum(int(A[i,j])*b[j] for j in range(n))
  m.Add(da==kin).OnlyEnforceIf(a[i]); m.Add(db==kin).OnlyEnforceIf(b[i])
  m.Add(db==cross).OnlyEnforceIf(a[i]); m.Add(da==cross).OnlyEnforceIf(b[i])
 pos=[]; neg=[]
 for i in range(n):
  eq=m.NewBoolVar(f'e{i}'); pp=m.NewBoolVar(f'p{i}'); nn=m.NewBoolVar(f'n{i}')
  m.Add(eq+pp+nn==1-a[i]-b[i])
  da=sum(int(A[i,j])*a[j] for j in range(n)); db=sum(int(A[i,j])*b[j] for j in range(n))
  m.Add(da==db).OnlyEnforceIf(eq)
  m.Add(da==ell).OnlyEnforceIf(pp); m.Add(db==0).OnlyEnforceIf(pp)
  m.Add(da==0).OnlyEnforceIf(nn); m.Add(db==ell).OnlyEnforceIf(nn)
  pos.append(pp); neg.append(nn)
 m.Add(sum(pos)+sum(neg)>=1)
 s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=seconds; s.parameters.num_search_workers=workers; s.parameters.random_seed=seed; s.parameters.randomize_search=True; s.parameters.symmetry_level=3
 t=time.time(); st=s.Solve(m); out={'status':s.StatusName(st),'seconds':time.time()-t,'ell':ell,'kin':kin,'cross':cross}
 if st in (cp_model.OPTIMAL,cp_model.FEASIBLE):
  C1=[i for i in range(n) if s.Value(a[i])]; C2=[i for i in range(n) if s.Value(b[i])]; U=set(C1+C2); B=A.copy(); changed=[]
  for v in range(n):
   if v in U: continue
   d1=int(A[v,C1].sum());d2=int(A[v,C2].sum())
   if (d1,d2) in ((ell,0),(0,ell)):
    B[v,C1+C2]=1-B[v,C1+C2];B[C1+C2,v]=B[v,C1+C2];changed.append(v)
  common=B.astype(np.int64)@B.astype(np.int64); exact=bool(np.all(B.sum(1)==60) and np.all(common[B==1]==29) and np.all(common[(B==0)&(~np.eye(n,dtype=bool))]==30))
  C=gs.conf_from_graph(B); conf=bool(np.array_equal(C@C,121*np.eye(122,dtype=np.int64)))
  out.update(C1=C1,C2=C2,changed=changed,srg=exact,conference=conf,sha256=hashlib.sha256(C.astype(np.int8).tobytes()).hexdigest())
  if conf:
   ep=pe.solve(C,'eigplus',seconds,workers);em=pe.solve(C,'eigminus',seconds,workers)
   out['plus_status']=ep['status'];out['minus_status']=em['status']
   if 'bits' in ep: out['plus_negative']=[i for i,x in enumerate(ep['bits']) if x]
   if 'bits' in em: out['minus_negative']=[i for i,x in enumerate(em['bits']) if x]
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('which',choices=['paley','peisert']);ap.add_argument('ell',type=int);ap.add_argument('kin',type=int);ap.add_argument('cross',type=int);ap.add_argument('--seconds',type=float,default=60);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seed',type=int,default=1)
 q=ap.parse_args();print('RESULT_JSON='+json.dumps({'which':q.which,**run(gs.cayley(q.which),q.ell,q.kin,q.cross,q.seconds,q.workers,q.seed)},sort_keys=True))
if __name__=='__main__':main()
