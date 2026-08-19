import numpy as np,sys,time,itertools
sys.path.insert(0,'/tmp/e1work')
sys.path.insert(0,'/tmp/claude-1000/-claude/e3c53bf5-2bab-41d0-9b30-ce6e3a2d6316/scratchpad/qml-verify/src')
from minmax_quadratic import paley_conference_prime_power
from fractions import Fraction

p=int(sys.argv[1]) if len(sys.argv)>1 else 11
q=p*p;n=q+1
if p==11:
    Yfull=np.load('/tmp/e1work/maxplus_p11_eps1.npy').astype(np.int8)
else:
    Y=np.load(f'/tmp/maxplus_p{p}.npy').astype(np.int8)
    Yfull=Y[Y[:,0]==1]
Nh=len(Yfull)
C=paley_conference_prime_power(p).astype(np.int64)
iu=np.triu_indices(n,1)
npair=len(iu[0])
pid=-np.ones((n,n),dtype=np.int64)
pid[iu[0],iu[1]]=np.arange(npair); pid[iu[1],iu[0]]=pid[iu[0],iu[1]]
t0=time.time()
Q=(Yfull[:,iu[0]]*Yfull[:,iu[1]]).astype(np.float32)
M=(Q.T@Q)/np.float32(Nh)   # npair x npair, 4-point moments (and 2-point on diag)
M=M.astype(np.float64)
print(f"M4 built {M.shape} in {time.time()-t0:.1f}s")
quads=np.array(list(itertools.combinations(range(n),4)),dtype=np.int64)
i,j,k,l=quads.T
m1=M[pid[i,j],pid[k,l]]; m2=M[pid[i,k],pid[j,l]]; m3=M[pid[i,l],pid[j,k]]
consist=max(np.abs(m1-m2).max(),np.abs(m1-m3).max())
print(f"pairing consistency max diff: {consist:.2e}")
m4=m1
kap=C[i,j]*C[k,l]+C[i,k]*C[j,l]+C[i,l]*C[j,k]
L=Fraction(p-2,2*p*p); T=Fraction(p-2,p*(2*p-1))
sel=np.abs(kap)==1
mu=np.abs(m4[sel]).max()
print(f"max |mu| over |kappa|=1 four-sets: {mu:.9f}  = {Fraction(round(mu*Nh)),Nh} check {mu*Nh:.3f}")
print(f"L=(p-2)/2p^2={float(L):.9f}  |mu|<=L ? {mu<=float(L)+1e-12}")
print(f"|T|=(p-2)/(p(2p-1))={float(T):.9f}  |mu|<|T| ? {mu<float(T)-1e-12}")
print(f"max |m4| over ALL four-sets: {np.abs(m4).max():.9f}")
# G = chi(u-v)*mu > T check on through-e |kappa|=1 (e = {inf,0} = indices 0,1):
selE=sel&(i==0)&(j==1)
Gm=(C[k[selE],l[selE]]*m4[selE])
print(f"through-e |kappa|=1 sets: {selE.sum()}  min G = {Gm.min():.9f}  > -|T| ? {Gm.min()>-float(T)+1e-12}")
np.save(f'/tmp/e1work/m4diag_p{p}.npy',np.diag(M))
