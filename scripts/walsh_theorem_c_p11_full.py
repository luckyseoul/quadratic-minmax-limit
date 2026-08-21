#!/usr/bin/env python3
"""FULL-ensemble Theorem C check at p=11, properly vectorized.
Rows packed as (lo:uint64, hi:uint64) pairs (122 bits). Elimination is
sequential (inherent to Gaussian elimination) but O(basis_size) uint64 XORs
per row via plain Python ints from numpy scalars -- cheap. Consistency
checks (pass 2) use np.bitwise_count for whole-chunk vectorized popcount."""
import numpy as np, sys, time
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/src')
sys.path.insert(0,'/home/nick/quadratic-minmax-limit/scripts/maxplus_profile_enum')
from kgen import field_ctx

p=11; q=p*p; n=q+1
q_,mul,chi,tr=field_ctx(p)
def sub(u,v): return (u%p-v%p)%p + ((u//p-v//p)%p)*p
ONE=1
C=np.zeros((n,n),dtype=np.int64); C[0,1:]=1; C[1:,0]=1
for e1 in range(q):
    for e2 in range(q):
        if e1!=e2: C[1+e1,1+e2]=chi(sub(e1,e2))
def order_of(e):
    x,o=e,1
    while x!=ONE: x=mul(x,e); o+=1
    return o
gen=next(e for e in range(2,q) if order_of(e)==q-1)
pi=np.zeros(n,dtype=np.int64); pi[0]=0
for e in range(q): pi[1+e]=1+mul(e,gen)
d=np.zeros(n,dtype=np.int64); d[0]=1
d[1:]=-C[pi[0],pi[1:]]*C[0,1:]
ok=(-d[:,None]*d[None,:]*C==C[np.ix_(pi,pi)]); np.fill_diagonal(ok,True); assert ok.all()

A=np.load('/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy',mmap_mode='r')
Ntot=A.shape[0]
i,j=0,1
CH=1_000_000
print(f"p=11: streaming {Ntot} rows, n={n}, edge=({i},{j})",flush=True)

# bit-plane weights for packing 122 bits -> (lo64, hi64)
w_lo=(np.uint64(1)<<np.arange(64,dtype=np.uint64))
w_hi=(np.uint64(1)<<np.arange(n-64,dtype=np.uint64))

def to_minus_bits_lohi(chunk_i64):
    Ym=(d[None,:]*chunk_i64[:,pi])
    B=((1-Ym)//2).astype(np.uint64)                  # 0/1, shape (m,n)
    fe=C[i,j]*Ym[:,i]*Ym[:,j]
    lo=(B[:,:64]*w_lo[None,:]).sum(axis=1,dtype=np.uint64)
    hi=(B[:,64:]*w_hi[None,:]).sum(axis=1,dtype=np.uint64)
    return lo,hi,fe<0

MASKN=np.uint64((1<<n)-1) if n<64 else None  # n=122>64, full 128 bits used

t0=time.time()
# ---- Pass 1: build echelon basis of B_U (pivot bit -> (lo,hi) python ints)
basis={}   # pivot -> (lo,hi)
def pivot_bit(lo,hi):
    if hi: return 64+int(hi).bit_length()-1
    if lo: return int(lo).bit_length()-1
    return -1
def add_row(lo,hi):
    while True:
        pb=pivot_bit(lo,hi)
        if pb<0: return False
        if pb in basis:
            blo,bhi=basis[pb]; lo^=blo; hi^=bhi
        else:
            basis[pb]=(lo,hi); return True

nU=nUc=0
for lo0 in range(0,Ntot,CH):
    chunk=A[lo0:lo0+CH].astype(np.int64)
    lo,hi,mask=to_minus_bits_lohi(chunk)
    loU,hiU=lo[mask],hi[mask]
    nU+=int(mask.sum()); nUc+=int((~mask).sum())
    for a,b in zip(loU.tolist(),hiU.tolist()):
        add_row(int(a),int(b))
    if (lo0//CH)%10==0:
        print(f"  pass1 {lo0}/{Ntot}  rank={len(basis)}  t={time.time()-t0:.0f}s",flush=True)
rankU=len(basis)
print(f"pass1 done: |U|={nU} |Uc|={nUc} rank(B_U)={rankU} (n={n}, n/2={n//2}) [{time.time()-t0:.0f}s]",flush=True)

# ---- RREF the basis (reduce pivot rows against each other) -> nullspace + solve x1
piv=sorted(basis.keys(),reverse=True)
rows=[basis[b] for b in piv]
for a in range(len(rows)):
    alo,ahi=rows[a]
    for bidx in range(len(rows)):
        if bidx==a: continue
        blo,bhi=rows[bidx]
        pb=piv[a]
        bit = (bhi>>(pb-64))&1 if pb>=64 else (blo>>pb)&1
        if bit:
            rows[bidx]=(blo^alo, bhi^ahi)
pivset=set(piv)
prow={piv[k]:k for k in range(len(piv))}
def getbit(lo,hi,c):
    return (hi>>(c-64))&1 if c>=64 else (lo>>c)&1
free=[c for c in range(n) if c not in pivset]
null_basis=[]           # list of (lo,hi) python-int pairs
for f in free:
    vlo = 1<<f if f<64 else 0
    vhi = 1<<(f-64) if f>=64 else 0
    for pv in piv:
        rlo,rhi=rows[prow[pv]]
        if getbit(rlo,rhi,f):
            if pv>=64: vhi|=1<<(pv-64)
            else: vlo|=1<<pv
    null_basis.append((vlo,vhi))
print(f"nullspace dim = {len(null_basis)} (expect n-rank = {n-rankU}) [{time.time()-t0:.0f}s]",flush=True)

# ---- solve B_U x1 = all-ones via augmented elimination (bit n = RHS)
AUGBIT=n
basis2={}
def add_row_aug(lo,hi,rhs):
    v_lo,v_hi=lo,hi
    aug=rhs
    while True:
        pb=pivot_bit(v_lo,v_hi)
        if pb<0:
            if aug: return False   # inconsistent row: 0 = 1
            return True
        if pb in basis2:
            blo,bhi,baug=basis2[pb]
            v_lo^=blo; v_hi^=bhi; aug^=baug
        else:
            basis2[pb]=(v_lo,v_hi,aug); return True
consistent=True
for lo0 in range(0,Ntot,CH):
    chunk=A[lo0:lo0+CH].astype(np.int64)
    lo,hi,mask=to_minus_bits_lohi(chunk)
    loU,hiU=lo[mask],hi[mask]
    for a,b in zip(loU.tolist(),hiU.tolist()):
        if not add_row_aug(int(a),int(b),1):
            consistent=False
print(f"pass1b (augmented, solvability) done: consistent={consistent} [{time.time()-t0:.0f}s]",flush=True)

x1=(0,0)
if consistent:
    piv2=sorted(basis2.keys(),reverse=True)
    rows2=[basis2[b] for b in piv2]
    for a in range(len(rows2)):
        alo,ahi,aaug=rows2[a]
        for bidx in range(len(rows2)):
            if bidx==a: continue
            blo,bhi,baug=rows2[bidx]
            pb=piv2[a]
            bit=(bhi>>(pb-64))&1 if pb>=64 else (blo>>pb)&1
            if bit:
                rows2[bidx]=(blo^alo,bhi^ahi,baug^aaug)
    x1lo=x1hi=0
    for k,pv in enumerate(piv2):
        rlo,rhi,raug=rows2[k]
        if raug:
            if pv>=64: x1hi|=1<<(pv-64)
            else: x1lo|=1<<pv
    x1=(x1lo,x1hi)
print(f"x1 = {x1}  [{time.time()-t0:.0f}s]",flush=True)

# ---- Pass 2: vectorized check on U^c: for each check vector, dot-product
# parity must be CONSTANT across all rows of U^c.
check_vecs = null_basis + ([x1] if consistent else [])
seen = [set() for _ in check_vecs]
vlo_arr=np.array([v[0] for v in check_vecs],dtype=np.uint64)
vhi_arr=np.array([v[1] for v in check_vecs],dtype=np.uint64)
for lo0 in range(0,Ntot,CH):
    chunk=A[lo0:lo0+CH].astype(np.int64)
    lo,hi,mask=to_minus_bits_lohi(chunk)
    loC,hiC=lo[~mask],hi[~mask]
    if len(loC)==0: continue
    for vi in range(len(check_vecs)):
        andlo = loC & vlo_arr[vi]
        andhi = hiC & vhi_arr[vi]
        par = (np.bitwise_count(andlo) + np.bitwise_count(andhi)) & 1
        vals=set(np.unique(par).tolist())
        seen[vi]|=vals
        if len(seen[vi])>1:
            pass   # already mixed; keep scanning others but could early-exit per-vec
    if (lo0//CH)%10==0:
        print(f"  pass2 {lo0}/{Ntot}  t={time.time()-t0:.0f}s",flush=True)

ker_mixed=sum(1 for s in seen[:len(null_basis)] if len(s)>1)
aff_mixed=1 if (consistent and len(seen[-1])>1) else 0
closed = consistent and ker_mixed==0 and aff_mixed==0
print(f"\np=11 FULL-ENSEMBLE Theorem C (all {Ntot} points, exact):")
print(f"  |U|={nU} |Uc|={nUc}  rank(B_U)={rankU}  ker_dim={len(null_basis)}")
print(f"  solvable={consistent}  ker_mixed={ker_mixed}  aff_mixed={aff_mixed}")
print(f"  CLOSED (Theorem C holds at p=11, full ensemble): {closed}")
print(f"  total time {time.time()-t0:.0f}s")
