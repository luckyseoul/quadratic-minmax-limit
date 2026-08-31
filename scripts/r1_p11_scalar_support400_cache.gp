\\ Build a pivot-normalized p=11 scalar subspace satisfying the lattice
\\ support condition n == 0 or 3 (mod 4) through exponent 400.
\\ Environment: P11_SCALAR_CACHE, P11_OUTPUT_DIRECTORY.
default(parisize, 12G);

scalar_cache = getenv("P11_SCALAR_CACHE");
output_directory = getenv("P11_OUTPUT_DIRECTORY");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(output_directory) != "t_STR", error("missing P11_OUTPUT_DIRECTORY"));

Y = read(scalar_cache);
X = Y[#Y];
mf = X[1];
K0 = X[2];
S = mfsturm(mf);
H = 400;
E = 200;

gettime();
M = mfcoefs(mf, H) * K0;
print("COEFFICIENT_MS=", gettime());
forbidden = List();
for(e = 0, H, if(e % 4 == 1 || e % 4 == 2, listput(forbidden, e + 1)));
forbidden = Vec(forbidden);
J = matker(vecextract(M, forbidden, Str("1..", matsize(M)[2])));
print("SUPPORT_DIM=", matsize(J)[2], " SUPPORT_CUTOFF=", H);
if(matsize(J)[2] != 60, error("unexpected support-kernel dimension"));

Q = M * J;
Qsturm = vecextract(Q, Str("1..", S + 1), Str("1..", matsize(Q)[2]));
idx = matindexrank(Qsturm);
RI = idx[1];
CJ = idx[2];
if(#RI != 60 || #CJ != 60, error("support subspace lacks full Sturm rank"));
R = vecextract(Qsturm, RI, CJ);
Jred = vecextract(J, Str("1..", matsize(J)[1]), CJ) * R^-1;
K = K0 * Jred;
Qred = M * Jred;
if(vecextract(Qred, RI, Str("1..", 60)) != matid(60), error("support pivot normalization failed"));
if(vecextract(Qred, forbidden, Str("1..", 60)) != matrix(#forbidden, 60), error("forbidden support survives"));
print("PIVOT_INDICES=", RI);

cache = Str(output_directory, "/p11_scalar_support400_mf_k_reduced_v1_20260827.gpbin");
writebin(cache, [mf, K, RI, H]);
path = Str(output_directory, "/p11_scalar_support400_qrows_exact_v1_20260827.txt");
for(i = 1, E + 1, write(path, Qred[i,]));
print("CACHE_WRITTEN=", cache);
print("QROWS_WRITTEN=", path, " ROWS=", E + 1, " COLS=", matsize(K)[2]);
quit;
