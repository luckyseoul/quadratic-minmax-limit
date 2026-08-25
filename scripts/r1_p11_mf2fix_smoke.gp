\\ Smoke test the guarded PARI half-integral slash patch on one genuine
\\ p=11 Kohnen basis vector at the formerly failing cusp.
default(parisize, 12G);
default(realprecision, 250);

X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
F = mflinear(mf, K[,1]);
params = 0;
c = mfslashexpansion(mf, F, [1,0;2,1], 15, 1, &params);
print("CUSP_HALF_PARAMS=", params);
print("COEFFICIENT_COUNT=", #c);
print("CUSP_HALF_FIRST=", c[1]);
print("CUSP_HALF_TARGET_INDEX_15=", c[16]);
print("EXACT_ZERO_COUNT=", sum(j = 1, #c, c[j] == 0));
writebin("/home/nick/p11_half_reduced_col1_exact_250.gpbin", c);
print("CACHE_WRITTEN=/home/nick/p11_half_reduced_col1_exact_250.gpbin");
quit;
