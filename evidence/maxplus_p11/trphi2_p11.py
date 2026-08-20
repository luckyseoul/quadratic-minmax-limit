"""Exact tr(Phi^2) at p=11 via the integer pair-moment Gram.

tr(Phi^2) = 4*||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2,   M = G/Nh, G = Q^T Q integer.
Verified exactly against known spectra at p=5,7.

G entries are bounded by Nh=3.7457e7 so G_ab^2 <= 1.4e15; a full row sum can
reach 1.04e19 which overflows int64, so rows are summed in halves and combined
as Python ints (exact).
"""
import numpy as np, time
from fractions import Fraction as F
import cupy as cp

p=11; q=p*p; n=q+1
Y=np.load('/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy')
Nh=len(Y)
iu=np.triu_indices(n,1); npair=len(iu[0])
print(f"p={p} Nh={Nh} n={n} npair={npair}",flush=True)

iu0=cp.asarray(iu[0].astype(np.int32)); iu1=cp.asarray(iu[1].astype(np.int32))
G=np.zeros((npair,npair),dtype=np.int64)
CH=50_000; t0=time.time(); nch=(Nh+CH-1)//CH
for ci,lo in enumerate(range(0,Nh,CH)):
    Yc=cp.asarray(Y[lo:lo+CH])
    Qc=(Yc[:,iu0]*Yc[:,iu1]).astype(cp.float32)   # +/-1, chunk<=50k => GEMM exact
    G+=cp.asnumpy(Qc.T@Qc).astype(np.int64)
    del Yc,Qc
    if (ci+1)%150==0 or ci+1==nch:
        cp.get_default_memory_pool().free_all_blocks()
        print(f"  {ci+1}/{nch}  {time.time()-t0:.0f}s",flush=True)

assert (np.diag(G)==Nh).all(), "diag(G) must equal Nh"
half=npair//2
tot=0
for r in range(npair):
    row=G[r]
    tot+=int(row[:half]@row[:half])+int(row[half:]@row[half:])
print(f"||G||_F^2 = {tot}",flush=True)

normM2=F(tot,Nh*Nh)
trPhi2=4*normM2-3*n**2+F(2*n**2*(n-1),p**2)
N=2*Nh
trGhat2=trPhi2*N*N
print(f"tr(Phi^2)  = {trPhi2} = {float(trPhi2):.6f}")
print(f"tr(Ghat^2) = {trGhat2}  integer={trGhat2.denominator==1}")
D=Nh//(2*p)
print(f"D={D}  A_p = D^2*tr(Phi^2) = {trPhi2*D*D}")
print(f"check tr(Ghat^2)==16p^2*A_p : {trGhat2==16*p*p*trPhi2*D*D}")
np.save('/mnt/storage/e1work/maxplus_p11/G_pairmoment_p11.npy',G)
