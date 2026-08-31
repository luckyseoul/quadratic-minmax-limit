\\ Extend the exact p=11 channel affine q-row exports through exponent 400.
\\ The final line remains the half-cusp target functional.
default(parisize, 10G);

H = read("/home/nick/p11_half_rational_250.gpbin");
Ctarget = H[2][16,];
Xmf = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = Xmf[1];
Kmf = Xmf[2];
E = 400;
gettime();
Mfull = mfcoefs(mf, E) * Kmf;
print("COEFFICIENT_MS=", gettime());
reductions = read("/home/nick/p11_channel_affine_qreductions_v2_20260827.gpbin");

for(c = 1, #reductions, {
  label = reductions[c][1];
  x0 = reductions[c][2];
  K = reductions[c][3];
  base = Mfull*x0;
  M = Mfull*K;
  path = Str("/home/nick/p11_affine_", label, "_qrows_exact400_20260827.txt");
  system(Str("rm -f -- ", path));
  for(i = 1, matsize(M)[1],
    write(path, concat([base[i]], vector(matsize(M)[2], j, M[i,j]))));
  write(path, concat([Ctarget*x0], vector(matsize(K)[2], j, (Ctarget*K)[j])));
  print(label, " ROWS=", matsize(M)[1], " FREE=", matsize(M)[2], " WROTE=", path);
});
quit;
