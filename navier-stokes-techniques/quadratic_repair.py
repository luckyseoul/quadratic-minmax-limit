"""Apply the imported cone and correction-cycle primitives to complete signings.

Finite local search only. Cone coefficients rank integral edge flips; they
never replace a signing or certify a Boolean norm. Every accepted move is
checked on all projective spin states, including previously inactive ones.
"""
from concurrent.futures import ProcessPoolExecutor
import argparse
import hashlib
import json
import os
from pathlib import Path
import tarfile
import time
import numpy as np
from correction_cycle.correction_cycle import CorrectionCycle
from stress_cone.stress_cone import general_cone_representation, NotInConeError

class NoImprovement(Exception):
    pass

def attempt(task):
    matrix, target, mode = task
    matrix=np.asarray(matrix,dtype=np.int64)
    n=len(matrix); i,j=np.triu_indices(n,1)
    masks=np.arange(1<<(n-1),dtype=np.int64)
    spins=np.ones((len(masks),n),dtype=np.int64)
    for c in range(1,n): spins[:,c]=1-2*((masks>>(c-1))&1)
    products=(spins[:,i]*spins[:,j]).T
    initial=matrix[i,j].copy()
    records=[]; evaluations=0
    def residual(a): return a@products
    def norm(q): return max(0,int(np.abs(q).max())-target)
    def solve(q,a):
        nonlocal evaluations
        before=int(np.abs(q).max())
        delta=-2*a[:,None]*products
        candidates=[(e,) for e in range(len(a))]
        cone_status='unused'
        if mode in ('index_pair','random_pair'):
            ranked=np.arange(len(a))[:12] if mode=='index_pair' else np.random.default_rng(731).permutation(len(a))[:12]
            candidates += [(int(e),int(f)) for pos,e in enumerate(ranked) for f in ranked[pos+1:]]
        if mode=='cone_pair':
            active=np.flatnonzero(np.abs(q)>target)
            signs=np.sign(q[active]); defect=np.abs(q[active])-target
            directions=np.concatenate([-delta[:,active]*signs,-np.eye(len(active))],axis=0)
            try:
                rep=general_cone_representation(directions,defect)
                assert np.max(np.abs(rep.residual()))<1e-6
                ranked=np.argsort(-rep.coefficients[:len(a)],kind='stable')[:12]
                candidates += [(int(e),int(f)) for pos,e in enumerate(ranked) for f in ranked[pos+1:]]
                cone_status='fractional ranking only'
            except NotInConeError:
                cone_status='no usable cone representation'
        best=before; chosen=None
        for edges in candidates:
            trial=q+delta[list(edges)].sum(axis=0)
            phi=int(np.abs(trial).max());evaluations+=1
            if phi<best: best=phi;chosen=edges
        records.append({'before':before,'after':best,'edges':None if chosen is None else [[int(i[e]),int(j[e])] for e in chosen],
                        'cone':cone_status,'evaluations':len(candidates)})
        if chosen is None:raise NoImprovement()
        return chosen
    last=[initial]
    def apply(a,edges):
        b=a.copy();b[list(edges)]*=-1;last[0]=b;return b
    cycle=CorrectionCycle(residual_fn=residual,linear_solve=solve,apply_correction=apply,
                          residual_norm=norm,tol=0.5,max_stages=8)
    start=time.monotonic()
    try: cycle.run(initial)
    except NoImprovement:pass
    final=last[0]
    # Independently recompute the accepted matrix using a direct pair loop.
    direct=np.zeros(len(masks),dtype=np.int64)
    for e in range(len(final)):direct+=final[e]*spins[:,i[e]]*spins[:,j[e]]
    assert np.array_equal(direct,residual(final))
    return {'mode':mode,'initial_phi':int(np.abs(residual(initial)).max()),'final_phi':int(np.abs(direct).max()),
            'target':target,'reached':norm(direct)==0,'records':records,'evaluations':evaluations,
            'seconds':time.monotonic()-start,'final_upper':final.tolist()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--workers',type=int,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--modes',nargs='+',choices=['single_edge','cone_pair','index_pair','random_pair'],default=['single_edge','cone_pair'])
    a=p.parse_args()
    if a.output.exists():raise ValueError('Refusing to overwrite evidence')
    if not 1<=a.workers<=len(os.sched_getaffinity(0)):raise ValueError('Invalid worker budget')
    with tarfile.open(a.archive) as tar:
        raw=tar.extractfile('original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json').read()
    d=json.loads(raw);B=np.array(d['aligned_full_target_matrix']);A=np.array(d['source_matrix']);k=len(A)
    tasks=[]
    for v in range(k,len(B)):
        ids=list(range(k))+[v];matrix=B[np.ix_(ids,ids)].copy();matrix[:k,:k]=A
        for mode in a.modes:tasks.append((matrix.tolist(),18,mode))
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=min(a.workers,len(tasks))) as pool:results=list(pool.map(attempt,tasks))
    result={'classification':'finite local correction experiment; no convergence assertion',
            'input_sha256':hashlib.sha256(raw).hexdigest(),'workers':min(a.workers,len(tasks)),
            'seconds':time.monotonic()-start,'results':results}
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'seconds':result['seconds'],'outcomes':[(r['mode'],r['initial_phi'],r['final_phi']) for r in results]}),flush=True)

if __name__=='__main__':main()
