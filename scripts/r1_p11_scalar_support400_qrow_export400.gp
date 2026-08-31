\\ Export every coefficient through 400 from the exact 60D support kernel.
default(parisize, 8G);
X = read("/home/nick/p11_scalar_support400_mf_k_reduced_v1_20260827.gpbin");
mf = X[1];
K = X[2];
E = 400;
gettime();
Q = mfcoefs(mf, E) * K;
print("COEFFICIENT_MS=", gettime());
path = "/home/nick/p11_scalar_support400_qrows_exact400_v1_20260827.txt";
system(Str("rm -f -- ", path));
for(i = 1, E + 1, write(path, Q[i,]));
print("QROWS_WRITTEN=", path, " ROWS=", E + 1, " COLS=", matsize(K)[2]);
quit;
