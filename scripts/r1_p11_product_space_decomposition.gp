\\ Exact decomposition of the p=11 scalar-theta product subspace into the
\\ Eisenstein, old-cuspidal, and new-cuspidal parts of weight 31, level 44.
\\ This is a diagnostic for a stable integral-weight cusp-zero transform.
default(parisize, 8G);

scalar_cache = getenv("P11_SCALAR_CACHE");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
X = read(scalar_cache);
mfhalf = X[1];
K = X[2];
D = matsize(K)[2];

mf = mfinit([44, 31, -11], 4);
mfnew = mfinit([44, 31, -11], 0);
mfold = mfinit([44, 31, -11], 2);
mfeis = mfinit([44, 31, -11], 3);
S = mfsturm(mf);
if(S != 183, error("unexpected integral Sturm bound", S));

\\ Fourier matrix of the 45-dimensional half-integral input, multiplied
\\ coefficientwise by theta via the exact Cauchy product.
H = mfcoefs(mfhalf, S) * K;
P = matrix(S + 1, D);
for(n = 0, S, for(a = 0, sqrtint(n), P[n + 1,] = P[n + 1,] + if(a, 2, 1) * H[n - a^2 + 1,]));

B = mfcoefs(mf, S);
idx = matindexrank(B);
RI = idx[1];
CJ = idx[2];
if(#RI != mfdim(mf), error("full integral basis has deficient Fourier rank"));
R = vecextract(B, RI, CJ);
to_full_coordinates(A) =
  R^-1 * vecextract(A, RI, Str("1..", matsize(A)[2]));

Cprod = to_full_coordinates(P);
Ceis = to_full_coordinates(mfcoefs(mfeis, S));
Cold = to_full_coordinates(mfcoefs(mfold, S));
Cnew = to_full_coordinates(mfcoefs(mfnew, S));
Cparts = matconcat([Ceis, Cold, Cnew]);
if(matsize(Cparts) != [mfdim(mf), mfdim(mf)], error("component dimensions do not sum to the full space"));
if(matrank(Cparts) != mfdim(mf), error("component coordinates are singular"));

parts = Cparts^-1 * Cprod;
de = mfdim(mfeis);
do = mfdim(mfold);
dn = mfdim(mfnew);
Peis = vecextract(parts, Str("1..", de), Str("1..", D));
Pold = vecextract(parts, Str(de + 1, "..", de + do), Str("1..", D));
Pnew = vecextract(parts, Str(de + do + 1, "..", de + do + dn), Str("1..", D));

print("PRODUCT_DIM=", matrank(Cprod));
print("EISENSTEIN_DIM=", de, " PRODUCT_PROJECTION_RANK=", matrank(Peis));
print("OLD_DIM=", do, " PRODUCT_PROJECTION_RANK=", matrank(Pold));
print("NEW_DIM=", dn, " PRODUCT_PROJECTION_RANK=", matrank(Pnew));
print("STACKED_COMPONENT_RANK=", matrank(parts));
quit;
