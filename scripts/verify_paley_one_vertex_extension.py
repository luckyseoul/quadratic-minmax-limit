"""Independent direct replay of the selected order-14 Paley extension."""
import json, sys
import numpy as np

result = json.load(open(sys.argv[1]))
q, n = 13, 14
chi = np.zeros(q, dtype=np.int8)
for a in range(1, q): chi[(a*a) % q] = 1
chi[chi == 0] = -1
C = np.zeros((n, n), dtype=np.int8)
C[0, 1:] = C[1:, 0] = 1
for i in range(q):
 for j in range(i+1, q): C[i+1,j+1] = C[j+1,i+1] = chi[(j-i) % q]
assert np.array_equal(C @ C, (n-1) * np.eye(n, dtype=np.int8))
b = np.ones(n, dtype=np.int8)
for i in range(1, n): b[i] = 1 if (result['b_mask'] >> (i-1)) & 1 else -1
values=[]
for mask in range(1 << n):
 x=np.array([1 if (mask>>i)&1 else -1 for i in range(n)], dtype=np.int8)
 values.append(abs(sum(int(C[i,j])*int(x[i])*int(x[j]) for i in range(n) for j in range(i+1,n))) + abs(int(b@x)))
assert max(values) == result['extension_phi']
print(json.dumps({'status':'PASS','source_phi':result['source_phi'],'extension_phi':max(values),'b_mask':result['b_mask']}))
