"""Independent direct-vector scoring of recovered paths and threshold reachability."""
import json
from fractions import Fraction
from pathlib import Path
import sys
import tarfile
import numpy as np

archive,result=map(Path,sys.argv[1:])
with tarfile.open(archive) as tar:
    d=json.load(tar.extractfile('original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json'))
z=json.loads(result.read_text()); rows=z['states']
B=np.array(d['aligned_full_target_matrix'],dtype=np.int64)
k=len(d['source_matrix']);m=len(B)-k;r=len(d['mismatch_edges'])
chosen=sorted({s for res in z['results'].values() for s in res['path']})
for s in chosen:
    ids=list(range(k))+[k+i for i in range(m) if s&(1<<i)]
    A=B[np.ix_(ids,ids)].copy();n=len(A)
    for b,(i,j) in enumerate(d['mismatch_edges']):
        if not s&(1<<(m+b)): A[i,j]*=-1;A[j,i]*=-1
    # Binary enumeration, direct pair sum: no Gray-code field updates.
    low,high=10**9,-10**9
    for start in range(0,1<<(n-1),16384):
        masks=np.arange(start,min(start+16384,1<<(n-1)),dtype=np.int64)
        x=np.ones((n,len(masks)),dtype=np.int64)
        for i in range(1,n):x[i]=1-2*((masks>>(i-1))&1)
        q=np.zeros(len(masks),dtype=np.int64)
        for i in range(n):
            for j in range(i+1,n):q+=A[i,j]*x[i]*x[j]
        low=min(low,int(q.min()));high=max(high,int(q.max()))
    assert [low,high]==rows[s][2:4],(s,low,high,rows[s])
cost=[Fraction(max(-row[2],row[3])**2,row[1]**3) for row in rows]
for mode,res in z['results'].items():
    threshold=Fraction(res['peak_squared'])
    path=res['path']
    assert path[0]==0 and path[-1]==len(rows)-1
    assert max(cost[s] for s in path)==threshold
    for old,new in zip(path,path[1:]):
        delta=new^old
        assert delta.bit_count()==1 and old&new==old
        if mode=='repairs_first' and delta<(1<<m): assert old>>m==(1<<r)-1
        if mode=='repairs_last' and delta>=(1<<m): assert old&((1<<m)-1)==(1<<m)-1
    # No path exists using only states STRICTLY below the claimed optimum.
    seen={0} if cost[0]<threshold else set()
    for s in range(len(rows)):
        if s not in seen:continue
        for b in range(m+r):
            if s&(1<<b):continue
            if mode=='repairs_first' and b<m and s>>m!=(1<<r)-1:continue
            if mode=='repairs_last' and b>=m and s&((1<<m)-1)!=(1<<m)-1:continue
            t=s|(1<<b)
            if cost[t]<threshold:seen.add(t)
    assert len(rows)-1 not in seen
print(json.dumps({'status':'PASS','independently_scored_path_states':len(chosen),'exact_threshold_obstructions':3}))
