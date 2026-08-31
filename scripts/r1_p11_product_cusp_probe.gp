\\ Probe PARI's half-integral slash failure through the integral-weight
\\ product F*theta.  The first three quotient coefficients at cusp 1/2
\\ vanish identically; this records the numerical debris left by PARI.
default(parisize, 12G);
default(realprecision, 250);

X = read("/home/nick/p11_mf_k.gpbin")[1];
mf = X[1];
K = X[2];
DK = matsize(K)[2];
anchor_coordinates = vector(DK, j, j)~;
F = mflinear(mf, K * anchor_coordinates);
FT = mfmul(F, mfTheta());
mfi = mfinit(FT);
print("HALF_PARAMS=", mfparams(mf));
print("PRODUCT_PARAMS=", mfparams(FT));
print("PRODUCT_SPACE_PARAMS=", mfparams(mfi));
print("PRODUCT_SPACE_DIM=", mfdim(mfi));

params = 0;
V = mfslashexpansion(mfi, FT, [1,0;2,1], 18, 0, &params);
print("PRODUCT_CUSP_PARAMS=", params);
for (j = 1, 8, print("PRODUCT_INDEX=", j - 1, " EXPONENT=", exponent(V[j]), " VALUE=", V[j]));
print("LEADING_SEPARATIONS_BITS=", vector(3, j, exponent(V[4]) - exponent(V[j])));
quit;
