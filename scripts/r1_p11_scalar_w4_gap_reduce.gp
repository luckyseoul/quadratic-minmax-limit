\\ Restrict the p=11 half-gap scalar space by the exact W_4 cusp-1/11 gap.
default(parisize, 8G);

scalar_cache = getenv("P11_SCALAR_CACHE");
w4_cache = getenv("P11_W4_CACHE");
output_directory = getenv("P11_W4_REDUCED_OUTPUT_DIRECTORY");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(w4_cache) != "t_STR", error("missing P11_W4_CACHE"));
if(type(output_directory) != "t_STR", error("missing P11_W4_REDUCED_OUTPUT_DIRECTORY"));

X = read(scalar_cache);
mf = X[1];
K0 = X[2];
Ctarget0 = X[5];
D0 = matsize(K0)[2];
W = read(w4_cache);
if(type(W[1]) == "t_VEC" && #W[1] == 5 && type(W[1][1]) == "t_MAT", W = W[#W]);
W4q = W[1];
if(matsize(W4q) != [184, D0], error("unexpected W4 q-row matrix"));

\\ Theta|W_4 is a nonzero scalar multiple of theta.  If the lattice theta
\\ factor has zero coefficients 1..4, its product with theta obeys these
\\ four rational relations after eliminating the unknown constant term.
G = matrix(4, D0);
G[1,] = W4q[2,] - 2 * W4q[1,];
G[2,] = W4q[3,];
G[3,] = W4q[4,];
G[4,] = W4q[5,] - 2 * W4q[1,];
if(matrank(G) != 4, error("unexpected exact W4 cusp-gap rank"));
J = matker(G);
if(matsize(J) != [D0, D0 - 4], error("unexpected W4 kernel dimension"));
Kgap = K0 * J;

S = mfsturm(mf);
E = 400;
gettime();
Qgap = mfcoefs(mf, E) * Kgap;
coefficient_ms = gettime();
Qsturm = vecextract(Qgap, Str("1..", S + 1), Str("1..", matsize(Qgap)[2]));
idx = matindexrank(Qsturm);
RI = idx[1];
CJ = idx[2];
D = matsize(J)[2];
if(#RI != D || #CJ != D, error("W4-reduced scalar space lacks full Sturm rank"));
R = vecextract(Qsturm, RI, CJ);
T = vecextract(J, Str("1..", D0), CJ) * R^-1;
K = K0 * T;
Q = mfcoefs(mf, E) * K;
if(G * T != matrix(4, D), error("exact W4 kernel reduction failed"));
if(vecextract(Q, RI, Str("1..", D)) != matid(D), error("W4 pivot normalization failed"));
Ctarget = Ctarget0 * T;

cache = Str(output_directory, "/p11_scalar_w4gap_support400_mf_k_reduced_v1_20260828.gpbin");
path = Str(output_directory, "/p11_scalar_w4gap_support400_qrows_exact400_v1_20260828.txt");
target_path = Str(output_directory, "/p11_scalar_w4gap_target_exact_v1_20260828.txt");
gap_path = Str(output_directory, "/p11_scalar_w4gap_rows_exact_v1_20260828.txt");
file_exists(candidate) = {
  my(fd);
  iferr(fd = fileopen(candidate, "r"), E, return(0));
  fileclose(fd);
  return(1);
};
if(file_exists(cache) || file_exists(path) || file_exists(target_path) || file_exists(gap_path), error("refusing to append to existing W4 output"));
writebin(cache, [mf, K, RI, E, Ctarget]);
for(i = 1, E + 1, write(path, Q[i,]));
for(i = 1, matsize(Ctarget)[1], write(target_path, Ctarget[i,]));
for(i = 1, 4, write(gap_path, G[i,]));
print("EXACT_W4_GAP_RANK=4");
print("REDUCED_DIMENSION=", D);
print("PIVOT_INDICES=", RI);
print("COEFFICIENT_MS=", coefficient_ms);
print("CACHE_WRITTEN=", cache);
print("QROWS_WRITTEN=", path);
print("TARGET_ROWS_WRITTEN=", target_path);
print("GAP_ROWS_WRITTEN=", gap_path);
quit;
