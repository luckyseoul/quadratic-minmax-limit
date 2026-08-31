\\ Measure stabilization of the p=11 scalar 0,3 mod 4 support subspace.
default(parisize, 12G);
Y = read("/home/nick/p11_scalar_full_mf_k_reduced.gpbin");
X = Y[#Y];
mf = X[1];
K = X[2];
H = 400;
gettime();
M = mfcoefs(mf, H) * K;
print("COEFFICIENT_MS=", gettime(), " ROWS=", H + 1, " BASE_DIM=", matsize(K)[2]);
cutoffs = [180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400];
for(j = 1, #cutoffs, C = cutoffs[j]; forbidden = List(); for(e = 0, C, if(e % 4 == 1 || e % 4 == 2, listput(forbidden, e + 1))); forbidden = Vec(forbidden); F = vecextract(M, forbidden, Str("1..", matsize(M)[2])); print("CUTOFF=", C, " SUPPORT_DIM=", matsize(matker(F))[2], " FORBIDDEN_ROWS=", #forbidden));
quit;
