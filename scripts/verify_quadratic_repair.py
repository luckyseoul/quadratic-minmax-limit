"""Replay all repair trajectories with independent direct integer spin sums."""
import itertools
import json
import sys
import tarfile
import numpy as np

with tarfile.open(sys.argv[1]) as tar:
    d=json.load(tar.extractfile('original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json'))
B=np.array(d['aligned_full_target_matrix']);A=np.array(d['source_matrix']);k=len(A)
x=np.array([(1,)+tail for tail in itertools.product((-1,1),repeat=k)],dtype=np.int64)
checks=0
def phi(M):
    q=np.zeros(len(x),dtype=np.int64)
    for i in range(len(M)):
        for j in range(i+1,len(M)):q+=M[i,j]*x[:,i]*x[:,j]
    return int(np.abs(q).max())
for filename in sys.argv[2:]:
    data=json.load(open(filename)); seen={}
    for result in data['results']:
        mode=result['mode'];v=k+seen.get(mode,0);seen[mode]=v-k+1
        ids=list(range(k))+[v];M=B[np.ix_(ids,ids)].copy();M[:k,:k]=A
        assert phi(M)==result['initial_phi']
        for step in result['records']:
            assert phi(M)==step['before']
            if step['edges'] is not None:
                for i,j in step['edges']:M[i,j]*=-1;M[j,i]*=-1
                assert step['after']<step['before']
            assert phi(M)==step['after'];checks+=1
        assert phi(M)==result['final_phi']
        assert M[np.triu_indices(len(M),1)].tolist()==result['final_upper']
        assert result['reached']==(phi(M)<=result['target'])
    assert all(count==10 for count in seen.values())
print(json.dumps({'status':'PASS','replayed_stages':checks,'spin_states_per_norm':len(x)}))
