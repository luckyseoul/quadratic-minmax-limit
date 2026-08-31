\\ Export the exact coefficient-extraction functional for the p=11 R1
\\ half-cusp target in the Fourier-pivot-normalized Kohnen basis.
default(parisize, 12G);
X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
RI = X[3];
H = read("/home/nick/p11_half_rational_250.gpbin");
Ctarget = H[2][16,];
base_target = H[3][16];
S = mfsturm(mf);
Q = mfcoefs(mf, S) * K;
print("STURM=", S);
print("PIVOT_INDICES_ZERO_BASED=", vector(#RI, i, RI[i] - 1));
print("PIVOT_IDENTITY=", vecextract(Q, RI, Str("1..", #RI)) == matid(#RI));
print("TARGET_FUNCTIONAL=", Ctarget);
print("TARGET_NONZERO=", sum(i = 1, #Ctarget, Ctarget[i] != 0));
print("TARGET_SIGNS_POS=", sum(i = 1, #Ctarget, Ctarget[i] > 0));
print("TARGET_SIGNS_NEG=", sum(i = 1, #Ctarget, Ctarget[i] < 0));
print("BASE_TARGET_EXACT=", base_target);
print("BASE_TARGET_NORM_EXACT=", norm(base_target));
quit;
