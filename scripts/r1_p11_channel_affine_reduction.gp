\\ Exact affine reductions for the three p=11 square-circle channels.
\\ They combine every proved cusp gap with the complete dual coefficients
\\ through exponent 28.  The kernel is then normalized by independent later
\\ Fourier coefficients.  Thus the free variables are actual q-coefficients,
\\ rather than the enormous coordinates returned directly by matker().
\\ Floating reduced coefficient rows are exported for reconnaissance and the
\\ same affine reductions are retained exactly for QSopt_ex.
default(parisize, 4G);
default(realprecision, 80);

p = 11; n = p^2 + 1; d = n/2;
archive = getenv("P11_MODULAR_DIRECTORY");
output_directory = getenv("P11_OUTPUT_DIRECTORY");
if(type(archive) != "t_STR", archive = "/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25");
if(type(output_directory) != "t_STR", error("missing P11_OUTPUT_DIRECTORY"));
H = read(Str(archive, "/p11_half_rational_250.gpbin"));
Z = read(Str(archive, "/p11_cusp_zero_exact_250.gpbin"));
Q = read(Str(archive, "/p11_cusp_quarter_exact_250.gpbin"));
R = read(Str(archive, "/p11_cusp_p_exact_250.gpbin"));
Minf = H[1];
Ctarget = H[2][16,];
Xmf = read(Str(archive, "/p11_mf_k_reduced.gpbin"))[1];
mf = Xmf[1];
Kmf = Xmf[2];
Mfull = mfcoefs(mf, mfsturm(mf))*Kmf;
if(Mfull[1..matsize(Minf)[1],] != Minf, error("infinity coefficient cache mismatch"));

stack(A, B) = matconcat([A~, B~])~;
file_exists(path) = {
  my(fd);
  iferr(fd = fileopen(path, "r"), E, return(0));
  fileclose(fd);
  return(1);
};
Acusp = H[2][1..15,];
Acusp = stack(Acusp, Z[1]);
Acusp = stack(Acusp, Q[1]);
Acusp = stack(Acusp, R[1]);

a2_shadow = 1/4 * (1 - (p - 1)^2/(d + 2)) - (p - 1)^2/(4*p*(d + 2));
\\ Props. 15.633--15.635 state the e=20 and e=24 operators at
\\ H(u/2)=H(u)/16.  Infinity q-rows are the unscaled H(u) theta series.
\\ The v2 reduction omitted this factor and therefore fixed two false rows.
second = 16 * [a2_shadow, a2_shadow + (p - 1)/(8*p), a2_shadow + (p + 1)/(8*p)];
first = -2/(d + 2);
third = 16 * 1/4 * (1 - (p + 1)^2/(d + 2));
fourth = [-(p^4 + 2*p^3 - 69*p^2 + 136*p + 26)/(p^2 + 5), (p^4 - 14*p^3 + 89*p^2 - 196*p + 24)/(p^2 + 5), (p^4 - 10*p^3 + 69*p^2 - 176*p - 76)/(p^2 + 5)];
labels = ["circle-kernel", "circle-low", "circle-high"];

reductions = vector(3);
for(c = 1, 3, {
  A = stack(Acusp, Minf[1..29,]);
  b = concat(vector(matsize(Acusp)[1]), vector(29))~;
  b[matsize(Acusp)[1] + 11 + 1] = first;
  b[matsize(Acusp)[1] + 20 + 1] = second[c];
  b[matsize(Acusp)[1] + 24 + 1] = third;
  b[matsize(Acusp)[1] + 27 + 1] = fourth[c];
  x0 = matinverseimage(A, b);
  if(x0 == [] || A*x0 != b, error("inconsistent channel", labels[c]));
  K = matker(A);
  B = Mfull*K;
  idx = matindexrank(B);
  RI = idx[1]; CJ = idx[2];
  if(#RI != matsize(K)[2] || #CJ != matsize(K)[2],
    error("unexpected free rank", [labels[c], #RI, #CJ, matsize(K), matsize(B), matrank(B)]));
  all_k_rows = Str("1..", matsize(K)[1]);
  all_k_cols = Str("1..", matsize(K)[2]);
  J = vecextract(B, RI, CJ);
  Kq = vecextract(K, all_k_rows, CJ) * J^-1;
  coeff0 = Mfull*x0;
  x0q = x0 - Kq*vecextract(coeff0, RI);
  Mq = Mfull*Kq;
  baseq = Mfull*x0q;
  if(A*x0q != b || A*Kq != 0, error("affine normalization failed", labels[c]));
  if(vecextract(Mq, RI, all_k_cols) != matid(#RI), error("q pivots not identity", labels[c]));
  if(vecextract(baseq, RI) != 0, error("q pivot base not zero", labels[c]));
  reductions[c] = [
    labels[c], x0q, Kq, matrank(A), RI,
    Ctarget*x0q, Ctarget*Kq
  ];

  \\ Each line is [constant, 32 free-coordinate coefficients].  The last
  \\ line is the target functional; preceding lines are exponents 0,1,... .
  path = Str(output_directory, "/p11_affine_", labels[c], "_qrows_float_unscaled_v3_20260828.txt");
  if(file_exists(path), error("refusing to append to existing output", path));
  for(i = 1, matsize(Mq)[1],
    write(path, concat([1.0*baseq[i]], vector(matsize(Mq)[2], j, 1.0*Mq[i,j]))));
  write(path, concat([1.0*(Ctarget*x0q)], vector(matsize(Kq)[2], j, 1.0*(Ctarget*Kq)[j])));
  max_q_exp = vecmax(vector(matsize(Mq)[2], j, vecmax(vector(matsize(Mq)[1], i, exponent(Mq[i,j])))));
  print(labels[c], " RANK=", matrank(A), " FREE=", matsize(Kq)[2],
        " QPIVOTS_ZERO_BASED=", vector(#RI, i, RI[i] - 1),
        " QMAXEXP=", max_q_exp, " WROTE=", path);
});
exact_path = Str(output_directory, "/p11_channel_affine_qreductions_unscaled_v3_20260828.gpbin");
if(file_exists(exact_path), error("refusing to append to existing output", exact_path));
writebin(exact_path, reductions);
print("EXACT_WROTE=", exact_path);
quit;
