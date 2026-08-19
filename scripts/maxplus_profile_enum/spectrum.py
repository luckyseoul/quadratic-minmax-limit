import numpy as np,sys,time
sys.path.insert(0,'/tmp/e1work')
sys.path.insert(0,'/tmp/claude-1000/-claude/e3c53bf5-2bab-41d0-9b30-ce6e3a2d6316/scratchpad/qml-verify/src')
from minmax_quadratic import paley_conference_prime_power
from fractions import Fraction

p=int(sys.argv[1]) if len(sys.argv)>1 else 11
q=p*p;n=q+1
if p==11:
    Yfull=np.load('/tmp/e1work/maxplus_p11_eps1.npy').astype(np.float64)
else:
    Y=np.load(f'/tmp/maxplus_p{p}.npy').astype(np.float64)
    Yfull=Y[Y[:,0]==1]
Nh=len(Yfull); N=2*Nh
print(f"p={p}: Nh={Nh} N={N} D=N/(4p)={Fraction(N,4*p)}")
C=paley_conference_prime_power(p)
ev,EV=np.linalg.eigh(C)
Vp=EV[:,ev>1e-8]   # n x n/2
d=Vp.shape[1]; assert d==n//2
Z61=Yfull@Vp        # Nh x d ;  z z^T features
# svec index
iu=np.triu_indices(d)
sq2=np.sqrt(2.0)
def svec_batch(Zb):
    # rows: svec(z z^T) with off-diag *sqrt2
    O=Zb[:,:,None]*Zb[:,None,:]
    F=O[:,iu[0],iu[1]]
    F[:, iu[0]!=iu[1]]*=sq2
    return F
# diag constraints: svec(v_x v_x^T)
Vc=np.zeros((n,len(iu[0])))
for x in range(n):
    vx=Vp[x]
    O=np.outer(vx,vx)
    r=O[iu]
    r[iu[0]!=iu[1]]*=sq2
    Vc[x]=r
u,s,vt=np.linalg.svd(Vc,full_matrices=True)
rank=(s>1e-9).sum()
Kb=vt[rank:].T     # svecdim x dimZ
dimZ=Kb.shape[1]
print(f"dimZ from constraints: {dimZ}  formula n(n-6)/8={n*(n-6)//8}")
# accumulate Phi_full = E[svec svec^T]
t0=time.time()
svd_dim=len(iu[0])
Phi_full=np.zeros((svd_dim,svd_dim))
step=4000
for lo in range(0,Nh,step):
    F=svec_batch(Z61[lo:lo+step])
    Phi_full+=F.T@F
Phi_full/=Nh
print(f"Phi_full built in {time.time()-t0:.1f}s")
PhiZ=Kb.T@Phi_full@Kb
lam=np.linalg.eigvalsh(PhiZ)
D=Fraction(N,4*p)
print(f"lambda_min={lam[0]:.9f}  lambda_max={lam[-1]:.9f}")
print(f"6-floor: lambda_min>=6 ? {lam[0]>=6}   lambda_max(K)=8-lmin={8-lam[0]:.9f} <=2 ? {8-lam[0]<=2}")
print(f"candidate 8(n-6)/n = {8*(n-6)/n:.9f}")
# cluster eigenvalues
vals=[]
tolc=1e-6
cur=[lam[0]]
for x in lam[1:]:
    if x-cur[-1]<tolc: cur.append(x)
    else: vals.append((np.mean(cur),len(cur))); cur=[x]
vals.append((np.mean(cur),len(cur)))
print("spectrum (value, mult, value*D):")
for v,mult in vals:
    vd=v*float(D)
    print(f"  {v:.9f}  x{mult}   D*lam={vd:.6f}  (round {round(vd)})  err={abs(vd-round(vd)):.2e}")
np.save(f'/tmp/e1work/phiZ_p{p}.npy',PhiZ)
