\\ Exact scalar p=11 support-space block at cusp 1/2.
\\ The first 15 transformed coefficients vanish for the odd-coset theta;
\\ Every fixed-index coefficient spans one cyclotomic phase line across the
\\ support basis.  Strip that phase numerically and reconstruct the exact
\\ rational row with high-precision clearance, as in the harmonic cache.
default(parisize, 12G);
default(realprecision, 250);
X = read("/home/nick/p11_scalar_support400_mf_k_reduced_v1_20260827.gpbin");
mf = X[1];
K = X[2];
D = matsize(K)[2];
forms = vector(D, j, mflinear(mf, K[,j]));
params = 0;
gettime();
batch = mfslashexpansion(mf, forms, [1,0;2,1], 15, 0, &params);
print("BATCH_MS=", gettime());
print("CUSP_HALF_PARAMS=", params);
if(type(batch) != "t_VEC" || #batch != D, error("unexpected scalar half batch"));
refs = vector(16);
phases = vector(16);
for(m = 1, 16, {
  best = 0;
  ref = 0;
  for(j = 1, D, if(abs(batch[j][m]) > best, {best = abs(batch[j][m]); ref = j;}));
  refs[m] = ref;
  phases[m] = if(ref, batch[ref][m], 1);
});
max_imag_exponent = -1000000;
max_residual_exponent = -1000000;
recover_ratio(m, j) =
{
  if(!refs[m], return(0));
  my(ratio = batch[j][m] / phases[m], rat = bestappr(real(ratio)));
  max_imag_exponent = max(max_imag_exponent, exponent(imag(ratio)));
  max_residual_exponent = max(max_residual_exponent, exponent(ratio - rat));
  return(rat);
}
Cexact = matrix(16, D, m, j, recover_ratio(m, j));
print("MAX_IMAG_EXPONENT=", max_imag_exponent);
print("MAX_RESIDUAL_EXPONENT=", max_residual_exponent);
if(max_imag_exponent > -400 || max_residual_exponent > -400,
  error("insufficient scalar half rational reconstruction clearance",
    [max_imag_exponent, max_residual_exponent]));
gap = vecextract(Cexact, Str("1..15"), Str("1..", D));
target = vecextract(Cexact, Str("16..16"), Str("1..", D));
print("EXACT_ROWS=", matsize(Cexact)[1], " COLS=", matsize(Cexact)[2]);
print("GAP_RANK=", matrank(gap));
print("TARGET_RANK=", matrank(target));
cache = "/home/nick/p11_scalar_half_gap_exact_v1_20260827.gpbin";
writebin(cache, [gap, target, refs, params]);
path = "/home/nick/p11_scalar_half_gap_exact_v1_20260827.txt";
system(Str("rm -f -- ", path));
for(i = 1, matsize(gap)[1], write(path, gap[i,]));
print("CACHE_WRITTEN=", cache);
print("GAP_ROWS_WRITTEN=", path);
quit;
