\\ Exact one-column smoke test for the scalar p=11 support space at cusp 1/2.
default(parisize, 8G);
default(realprecision, 250);
X = read("/home/nick/p11_scalar_support400_mf_k_reduced_v1_20260827.gpbin");
mf = X[1];
K = X[2];
F = mflinear(mf, K[,1]);
params = 0;
gettime();
c = mfslashexpansion(mf, F, [1,0;2,1], 15, 1, &params);
print("BATCH_MS=", gettime());
print("CUSP_HALF_PARAMS=", params);
print("COEFFICIENT_COUNT=", #c);
print("COEFFICIENTS=", c);
quit;
