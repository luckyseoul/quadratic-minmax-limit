\\ Build a Fourier-pivot-normalized Kohnen basis for the p=11 ordinary
\\ theta space.  This supports a joint modular cone in which the unknown
\\ ordinary theta coefficients replace direct lattice enumeration.
\\ Environment: P11_OUTPUT_DIRECTORY.
default(parisize, 12G);

output_directory = getenv("P11_OUTPUT_DIRECTORY");
if(type(output_directory) != "t_STR", error("missing P11_OUTPUT_DIRECTORY"));

gettime();
mf = mfinit([44, 61/2, 44], 4);
print("MFINIT_MS=", gettime());
S = mfsturm(mf);
print("FULL_DIM=", mfdim(mf));
print("FULL_STURM=", S);

gettime();
B = mfcoefs(mf, S);
forbidden = List();
for(e = 0, S, if(e % 4 == 1 || e % 4 == 2, listput(forbidden, e + 1)));
forbidden = Vec(forbidden);
all_columns = Str("1..", matsize(B)[2]);
K = matker(vecextract(B, forbidden, all_columns));
print("KOHNEN_KERNEL_MS=", gettime());
print("KOHNEN_DIM=", matsize(K)[2]);

gettime();
Q = B * K;
if(vecextract(Q, forbidden, Str("1..", matsize(Q)[2])) != matrix(#forbidden, matsize(Q)[2]), error("forbidden scalar coefficients survive"));
idx = matindexrank(Q);
RI = idx[1];
CJ = idx[2];
if (#RI != matsize(K)[2] || #CJ != matsize(K)[2], error("unexpected scalar Kohnen rank", [#RI, #CJ]));
all_rows = Str("1..", matsize(K)[1]);
R = vecextract(Q, RI, CJ);
Kred = vecextract(K, all_rows, CJ) * R^-1;
Qred = mfcoefs(mf, S) * Kred;
if (vecextract(Qred, RI, Str("1..", #RI)) != matid(#RI), error("scalar pivot normalization failed"));
print("REDUCTION_MS=", gettime());
print("PIVOT_INDICES=", RI);
cache = Str(output_directory, "/p11_scalar_full_mf_k_reduced.gpbin");
writebin(cache, [mf, Kred, RI, CJ]);
print("CACHE_WRITTEN=", cache);
quit;
