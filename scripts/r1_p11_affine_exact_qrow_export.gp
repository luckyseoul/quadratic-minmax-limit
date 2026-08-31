\\ Export the q-pivot-normalized p=11 channel affine spaces over Q.
\\ Every coefficient line is [constant, 32 free-coordinate entries].
\\ The final line is the half-cusp target functional in the same coordinates.
default(parisize, 4G);

H = read("/home/nick/p11_half_rational_250.gpbin");
Ctarget = H[2][16,];
Xmf = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = Xmf[1];
Kmf = Xmf[2];
Mfull = mfcoefs(mf, mfsturm(mf))*Kmf;
reductions = read("/home/nick/p11_channel_affine_qreductions_v2_20260827.gpbin");

for(c = 1, #reductions, {
  label = reductions[c][1];
  x0 = reductions[c][2];
  K = reductions[c][3];
  base = Mfull*x0;
  M = Mfull*K;
  path = Str("/home/nick/p11_affine_", label, "_qrows_exact_20260827.txt");
  system(Str("rm -f -- ", path));
  for(i = 1, matsize(M)[1],
    write(path, concat([base[i]], vector(matsize(M)[2], j, M[i,j]))));
  write(path, concat([Ctarget*x0], vector(matsize(K)[2], j, (Ctarget*K)[j])));
  print(label, " ROWS=", matsize(M)[1], " FREE=", matsize(M)[2], " WROTE=", path);
});
quit;
