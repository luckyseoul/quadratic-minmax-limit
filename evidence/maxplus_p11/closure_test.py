"""Completeness test for the assembled Max+ set.

Max+ = {y in {+-1}^n : Cy = py} is invariant under the automorphism group of C.
So applying a dilation/Frobenius element gamma to any solution must land back
in the set. If the k=6 gauge-orbit enumeration dropped states, closure fails.

Membership is tested by sorting the row-bytes once and using searchsorted.
"""
import numpy as np, sys, time
sys.path.insert(0, '/tmp/e1work')
from dilation import build_group

p = 11
q = p * p
n = q + 1

t0 = time.time()
Y = np.load('/tmp/e1work/maxplus_p11_eps1.npy')   # (Nh, n) int8, col 0 = y_inf
Nh = len(Y)
print(f"loaded {Y.shape}  ({time.time()-t0:.0f}s)", flush=True)

def keys_of(A):
    Ac = np.ascontiguousarray(A)
    return Ac.view(np.dtype((np.void, Ac.dtype.itemsize * Ac.shape[1]))).reshape(-1)

t1 = time.time()
K = keys_of(Y)
order = np.argsort(K, kind='stable')
Ksorted = K[order]
print(f"sorted keys ({time.time()-t1:.0f}s)", flush=True)

def in_set(rows):
    kk = keys_of(rows)
    idx = np.searchsorted(Ksorted, kk)
    idx = np.clip(idx, 0, len(Ksorted) - 1)
    return Ksorted[idx] == kk

group, _ = build_group(p)
print(f"group size {len(group)}", flush=True)

rng = np.random.default_rng(0)
sample_idx = rng.choice(Nh, 20000, replace=False)
S = Y[sample_idx]

# sanity: the sample itself must be in the set
base_ok = in_set(S).all()
print(f"sanity (sample in set): {base_ok}", flush=True)

fails = []
for gi in range(0, len(group), max(1, len(group) // 12)):
    perm, dirperm, alphas = group[gi]
    # act on finite coords only; y_inf (col 0) is fixed
    T = np.empty_like(S)
    T[:, 0] = S[:, 0]
    T[:, 1:] = S[:, 1:][:, perm]
    ok = in_set(T)
    nfail = int((~ok).sum())
    fails.append((gi, nfail))
    print(f"  gamma[{gi:3d}]: {len(S)-nfail}/{len(S)} images in set"
          f"{'  <-- CLOSURE FAILURE' if nfail else ''}", flush=True)

tot = sum(f for _, f in fails)
print(f"\nTOTAL closure failures: {tot} over {len(fails)} group elements x {len(S)} samples")
print("VERDICT:", "CLOSED (set consistent with completeness)" if tot == 0
      else "NOT CLOSED -> enumeration is missing solutions")
