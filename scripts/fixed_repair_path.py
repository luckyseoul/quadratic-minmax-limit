"""One archived source/target: exact minimax over all addition/repair orders.

This bounded diagnostic does not choose targets or establish convergence.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import time

MEMBER = 'original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--scorer', type=Path, required=True)
    p.add_argument('--workers', type=int, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('Refusing to overwrite evidence')
    with tarfile.open(a.archive) as tar:
        raw = tar.extractfile(MEMBER).read()
    d = json.loads(raw)
    B, A, repairs = d['aligned_full_target_matrix'],d['source_matrix'],d['mismatch_edges']
    k,N = len(A),len(B)
    m,r = N-k,len(repairs)
    assert [(i,j) for i in range(k) for j in range(i+1,k) if A[i][j]!=B[i][j]] == [tuple(e) for e in repairs]
    wire = f'{k} {N} {r}\n'+'\n'.join(' '.join(map(str,row)) for row in B+repairs)+'\n'
    start = time.monotonic()
    run = subprocess.run([str(a.scorer.resolve()),str(a.workers)],input=wire,text=True,capture_output=True,check=True)
    elapsed = time.monotonic()-start
    rows = [list(map(int,line.split())) for line in run.stdout.splitlines()]
    assert len(rows)==1<<(m+r)
    costs = []
    for s,n,lo,hi,w in rows:
        assert s==len(costs) and n==k+(s&((1<<m)-1)).bit_count()
        costs.append(Fraction(max(-lo,hi)**2,n**3))
    assert max(-rows[0][2],rows[0][3])==17
    assert max(-rows[-1][2],rows[-1][3])==44
    results={}
    for mode in ('interleaved','repairs_first','repairs_last'):
        best={0:costs[0]}; pred={}
        for s in range(1,len(rows)):
            candidates=[]
            for bit in range(m+r):
                if not s&(1<<bit): continue
                old=s^(1<<bit)
                if old not in best: continue
                if mode=='repairs_first' and bit<m and old>>m!=(1<<r)-1: continue
                if mode=='repairs_last' and bit>=m and old&((1<<m)-1)!=(1<<m)-1: continue
                candidates.append((max(best[old],costs[s]),old))
            if candidates: best[s],pred[s]=min(candidates)
        s=len(rows)-1; path=[s]
        while s: s=pred[s];path.append(s)
        path.reverse()
        results[mode]={'peak_squared':str(best[len(rows)-1]),'peak':float(best[len(rows)-1])**0.5,
                       'path':path,'norms':[max(-rows[s][2],rows[s][3]) for s in path],
                       'orders':[rows[s][1] for s in path]}
    a.output.write_text(json.dumps({'classification':'finite fixed-witness path optimum; convergence open',
        'input_sha256':hashlib.sha256(raw).hexdigest(),'scorer_sha256':hashlib.sha256(a.scorer.read_bytes()).hexdigest(),
        'workers':a.workers,'scoring_seconds':elapsed,'states':rows,'results':results},indent=2)+'\n')
    print(json.dumps({'scoring_seconds':elapsed,'results':results},indent=2),flush=True)

if __name__=='__main__': main()
