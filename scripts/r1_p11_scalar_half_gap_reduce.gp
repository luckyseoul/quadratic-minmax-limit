\\ Restrict the exact scalar support space by the cusp-1/2 gap and export
\\ a Fourier-pivot-normalized basis through exponent 400.
\\ Environment: P11_SCALAR_CACHE, P11_HALF_GAP_CACHE,
\\              P11_HALF_REDUCED_OUTPUT_DIRECTORY.
default(parisize, 8G);
scalar_cache = getenv("P11_SCALAR_CACHE");
gap_cache = getenv("P11_HALF_GAP_CACHE");
output_directory = getenv("P11_HALF_REDUCED_OUTPUT_DIRECTORY");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(gap_cache) != "t_STR", error("missing P11_HALF_GAP_CACHE"));
if(type(output_directory) != "t_STR", error("missing P11_HALF_REDUCED_OUTPUT_DIRECTORY"));

X = read(scalar_cache);
mf = X[1];
K0 = X[2];
H = read(gap_cache);
G = H[1];
Ctarget = H[2];
if(matsize(G) != [150, 60], error("unexpected scalar half-gap matrix size"));
if(matrank(G) != 15, error("unexpected scalar half-gap rank"));
J = matker(G);
if(matsize(J) != [60, 45], error("unexpected scalar half-gap kernel size"));
Kgap = K0 * J;

S = mfsturm(mf);
E = 400;
gettime();
Qgap = mfcoefs(mf, E) * Kgap;
print("COEFFICIENT_MS=", gettime());
Qsturm = vecextract(Qgap, Str("1..", S + 1), Str("1..", matsize(Qgap)[2]));
idx = matindexrank(Qsturm);
RI = idx[1];
CJ = idx[2];
if(#RI != 45 || #CJ != 45, error("half-gap scalar space lacks full Sturm rank"));
R = vecextract(Qsturm, RI, CJ);
T = vecextract(J, Str("1..", 60), CJ) * R^-1;
K = K0 * T;
Q = mfcoefs(mf, E) * K;
if(G * T != matrix(150, 45), error("half-gap kernel reduction failed"));
if(vecextract(Q, RI, Str("1..", 45)) != matid(45), error("half-gap pivot normalization failed"));
Ctarget_reduced = Ctarget * T;

cache = Str(output_directory, "/p11_scalar_halfgap_support400_mf_k_reduced_v1_20260828.gpbin");
system(Str("rm -f -- ", cache));
writebin(cache, [mf, K, RI, E, Ctarget_reduced]);
path = Str(output_directory, "/p11_scalar_halfgap_support400_qrows_exact400_v1_20260828.txt");
system(Str("rm -f -- ", path));
for(i = 1, E + 1, write(path, Q[i,]));
target_path = Str(output_directory, "/p11_scalar_halfgap_target_exact_v1_20260828.txt");
system(Str("rm -f -- ", target_path));
for(i = 1, matsize(Ctarget_reduced)[1], write(target_path, Ctarget_reduced[i,]));
print("HALF_GAP_RANK=15");
print("REDUCED_DIMENSION=", matsize(K)[2]);
print("PIVOT_INDICES=", RI);
print("CACHE_WRITTEN=", cache);
print("QROWS_WRITTEN=", path);
print("TARGET_ROWS_WRITTEN=", target_path);
quit;
