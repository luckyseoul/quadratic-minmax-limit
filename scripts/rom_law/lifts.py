import numpy as np, struct, subprocess, sys
from paley import paley
def write(J,name):
    n=J.shape[0]
    with open(name,'wb') as f:
        f.write(struct.pack('i',n)); f.write(np.ascontiguousarray(J,dtype=np.float64).tobytes())
def phi(A,name='lift.bin'):
    write(A,name); o=subprocess.check_output(['./gs2',name]).decode().split()
    return max(float(o[0]),-float(o[1]))/2
def lift(C,d):
    n=C.shape[0]; B=C+d*np.eye(n)
    return np.block([[C,B],[B,-C]])
