\\ Export floating rows for cone/LP searches; exact certification remains in GP.
default(parisize, 12G);
default(realprecision, 80);
X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
H = read("/home/nick/p11_half_rational_250.gpbin");
Ctarget = H[2][16,];
S = mfsturm(mf);
M = mfcoefs(mf, S) * K;
path = "/home/nick/p11_functional_float_rows_20260827_v2.txt";
for(i = 1, S + 1, write(path, vector(matsize(M)[2], j, 1.0 * M[i,j])));
write(path, vector(#Ctarget, j, 1.0 * Ctarget[j]));
print("WROTE=", path, " ROWS=", S + 2, " COLS=", matsize(M)[2]);
quit;
