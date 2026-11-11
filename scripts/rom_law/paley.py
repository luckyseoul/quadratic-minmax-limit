import numpy as np, struct, subprocess, sys
def paley(q):
    # prime q = 1 mod 4: conference matrix of order q+1
    chi=np.zeros(q,dtype=int)
    for a in range(1,q): chi[(a*a)%q]=1
    chi[chi==0]=-1; chi[0]=0
    C=np.zeros((q+1,q+1))
    C[0,1:]=1; C[1:,0]=1
    for i in range(q):
        for j in range(q):
            C[i+1,j+1]=chi[(i-j)%q]
    return C
def run(J,name='tmp.bin'):
    n=J.shape[0]
    with open(name,'wb') as f:
        f.write(struct.pack('i',n)); f.write(np.ascontiguousarray(J,dtype=np.float64).tobytes())
    out=subprocess.check_output(['./gs2',name]).decode().split()
    return float(out[0]),float(out[1])
if __name__=='__main__':
    for q in (5,13,17):
        C=paley(q); n=q+1
        assert np.allclose(C@C,(n-1)*np.eye(n))
        mx,mn=run(C)
        print(q,n,'Phi=',max(mx,-mn)/2,'spectral',n*np.sqrt(n-1)/2, 'rho',max(mx,-mn)/(n*np.sqrt(n-1)))
