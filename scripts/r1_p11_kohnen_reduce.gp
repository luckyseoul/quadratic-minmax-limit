\\ Replace PARI's badly conditioned kernel basis for the p=11 Kohnen space
\\ by an exact Fourier-pivot-normalized basis of the same column space.
default(parisize, 12G);

X = read("/home/nick/p11_mf_k.gpbin")[1];
mf = X[1];
K = X[2];
S = mfsturm(mf);
Q = mfcoefs(mf, S) * K;
idx = matindexrank(Q);
RI = idx[1];
CJ = idx[2];
if (#RI != matsize(K)[2] || #CJ != matsize(K)[2], error("unexpected Kohnen rank", [#RI,#CJ]));
all_rows = Str("1..", matsize(K)[1]);
all_cols = Str("1..", matsize(K)[2]);
R = vecextract(Q, RI, CJ);
T = R^-1;
Kred = vecextract(K, all_rows, CJ) * T;
Qred = mfcoefs(mf, S) * Kred;

old_height = vecmax(vector(matsize(K)[2], j, vecmax(vector(matsize(K)[1], i, exponent(K[i,j])))));
new_height = vecmax(vector(matsize(Kred)[2], j, vecmax(vector(matsize(Kred)[1], i, exponent(Kred[i,j])))));
q_height = vecmax(vector(matsize(Qred)[2], j, vecmax(vector(matsize(Qred)[1], i, exponent(Qred[i,j])))));
print("KOHNEN_DIM=", matsize(Kred)[2]);
print("PIVOT_INDICES=", RI);
print("PIVOT_IDENTITY=", vecextract(Qred, RI, all_cols) == matid(#RI));
print("OLD_COORDINATE_MAX_EXPONENT=", old_height);
print("REDUCED_COORDINATE_MAX_EXPONENT=", new_height);
print("REDUCED_Q_MAX_EXPONENT=", q_height);
print("REDUCED_Q_RANK=", matrank(Qred));
writebin("/home/nick/p11_mf_k_reduced.gpbin", [mf, Kred, RI, CJ]);
print("CACHE_WRITTEN=/home/nick/p11_mf_k_reduced.gpbin");
quit;
