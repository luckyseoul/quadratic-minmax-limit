import numpy as np
rng = np.random.default_rng(42)
for n in [16,20,24]:
    N=2**n; idx=np.arange(N)[:,None]; X=((idx>>np.arange(n)[None,:])&1*2-1).astype(np.float32)
    Cs=[]; 
    iters=40 if n<=20 else 20
    for _ in range(iters):
        Z=rng.choice([-1,1],size=(n,n)).astype(np.float32)
        mx=np.abs(X@Z).sum(axis=1).max()
        Cs.append(mx)
    Cs=np.array(Cs)
    print(f'n={n}: C_cross=mean(max||xZ||_1/n^{3/2})={Cs.mean()/n**1.5:.4f} p5={np.percentile(Cs,5)/n**1.5:.4f} max={Cs.max()/n**1.5:.4f}')
