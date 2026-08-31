\\ Export exact infinity-coefficient rows of the pivot-normalized p=11
\\ ordinary-theta Kohnen space for the joint raw-shell cone.
default(parisize, 12G);
Y = read("/home/nick/p11_scalar_full_mf_k_reduced.gpbin");
X = Y[#Y];
mf = X[1];
K = X[2];
S = mfsturm(mf);
H = 200; \\ harmonic-theta Sturm bound; retain all coupled shell rows
if (H < S, error("export bound is below the scalar Sturm bound"));
M = mfcoefs(mf, H) * K;
path = "/home/nick/p11_scalar_full_qrows_exact_v3_20260827.txt";
for(i = 1, H + 1, write(path, M[i,]));
print("WROTE=", path, " ROWS=", H + 1, " COLS=", matsize(M)[2], " SCALAR_STURM=", S);
quit;
