"""GPU chunked rewrite of moments.py.

The original materializes Q = Nh x npair float32 (37.4M x 7381 = 1.1 TB) and
dies. Here M = Q^T Q / Nh is accumulated in chunks on the V100.

Exactness: entries of Q are +/-1, so a per-chunk Q_c^T Q_c has integer entries
bounded by the chunk size (50k << 2^24), making each fp32 GEMM exact. V100 has
no TF32 path, so cuBLAS fp32 really is fp32. Chunk results accumulate in float64.
"""
import numpy as np, sys, time, itertools
sys.path.insert(0, '/tmp/e1work')
sys.path.insert(0, '/tmp/claude-1000/-claude/e3c53bf5-2bab-41d0-9b30-ce6e3a2d6316/scratchpad/qml-verify/src')
from minmax_quadratic import paley_conference_prime_power
from fractions import Fraction
import cupy as cp

p = int(sys.argv[1]) if len(sys.argv) > 1 else 11
q = p * p
n = q + 1

if p == 11:
    Yfull = np.load('/tmp/e1work/maxplus_p11_eps1.npy')          # int8, (Nh, n)
else:
    Y = np.load(f'/tmp/maxplus_p{p}.npy')
    Yfull = Y[Y[:, 0] == 1]
Nh = len(Yfull)
print(f"p={p}: Nh={Nh}  n={n}", flush=True)

C = paley_conference_prime_power(p).astype(np.int64)
iu = np.triu_indices(n, 1)
npair = len(iu[0])
print(f"npair={npair}", flush=True)

iu0_g = cp.asarray(iu[0].astype(np.int32))
iu1_g = cp.asarray(iu[1].astype(np.int32))

# ---- accumulate M = E[q_a q_b] over all solutions, chunked on GPU ----
t0 = time.time()
M_acc = np.zeros((npair, npair), dtype=np.float64)
CH = 50_000
nch = (Nh + CH - 1) // CH
for ci, lo in enumerate(range(0, Nh, CH)):
    Yc = cp.asarray(Yfull[lo:lo + CH])                  # (b, n) int8
    Qc = (Yc[:, iu0_g] * Yc[:, iu1_g]).astype(cp.float32)  # (b, npair), +/-1
    G = Qc.T @ Qc                                        # exact in fp32
    M_acc += cp.asnumpy(G).astype(np.float64)
    del Yc, Qc, G
    if (ci + 1) % 100 == 0 or ci + 1 == nch:
        cp.get_default_memory_pool().free_all_blocks()
        print(f"  chunk {ci+1}/{nch}  {time.time()-t0:.0f}s", flush=True)
M = M_acc / Nh
print(f"M4 built {M.shape} in {time.time()-t0:.1f}s", flush=True)

# ---- four-point moment consistency + mu bounds ----
pid = -np.ones((n, n), dtype=np.int64)
pid[iu[0], iu[1]] = np.arange(npair)
pid[iu[1], iu[0]] = pid[iu[0], iu[1]]

quads = np.array(list(itertools.combinations(range(n), 4)), dtype=np.int64)
i, j, k, l = quads.T
print(f"quads: {len(quads):,}", flush=True)

m1 = M[pid[i, j], pid[k, l]]
m2 = M[pid[i, k], pid[j, l]]
m3 = M[pid[i, l], pid[j, k]]
consist = max(np.abs(m1 - m2).max(), np.abs(m1 - m3).max())
print(f"pairing consistency max diff: {consist:.2e}", flush=True)

m4 = m1
kap = C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k]
L = Fraction(p - 2, 2 * p * p)
T = Fraction(p - 2, p * (2 * p - 1))
sel = np.abs(kap) == 1
mu = np.abs(m4[sel]).max()
print(f"max |mu| over |kappa|=1 four-sets: {mu:.9f}   (mu*Nh = {mu*Nh:.3f})", flush=True)
print(f"L=(p-2)/2p^2 = {float(L):.9f}   |mu| <= L ? {mu <= float(L)+1e-12}", flush=True)
print(f"|T|=(p-2)/(p(2p-1)) = {float(T):.9f}   |mu| < |T| ? {mu < float(T)-1e-12}", flush=True)
print(f"max |m4| over ALL four-sets: {np.abs(m4).max():.9f}", flush=True)

selE = sel & (i == 0) & (j == 1)
Gm = C[k[selE], l[selE]] * m4[selE]
print(f"through-e |kappa|=1 sets: {selE.sum()}   min G = {Gm.min():.9f}   > -|T| ? {Gm.min() > -float(T)+1e-12}", flush=True)

np.save(f'/tmp/e1work/m4diag_p{p}.npy', np.diag(M))
print("done", flush=True)
