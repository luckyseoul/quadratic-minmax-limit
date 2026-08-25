\\ Batch all p=11 Kohnen forms at cusp 1/2, strip the exact fixed-index
\\ cyclotomic phases, and recover the resulting rational constraint matrix.
default(parisize, 12G);
precision_text = getenv("P11_REALPRECISION");
real_precision = if(precision_text, eval(precision_text), 250);
default(realprecision, real_precision);

X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
DK = matsize(K)[2];
base = read("/home/nick/p11_half_reduced_col1_exact_250.gpbin");
z = exp(2 * Pi * I / 11);
phase = vector(16, m, subst(lift(base[m]), 't, z));
forms = vector(DK, j, mflinear(mf, K[,j]));

gettime();
params = 0;
batch = mfslashexpansion(mf, forms, [1,0;2,1], 15, 0, &params);
batch_time = gettime();
print("REAL_PRECISION=", real_precision);
print("KOHNEN_DIM=", DK);
print("CUSP_HALF_PARAMS=", params);
print("BATCH_TIME_MS=", batch_time);

max_imag_exponent = -1000000;
max_residual_exponent = -1000000;
max_numerator_exponent = -1000000;
max_denominator_exponent = -1000000;
recover_ratio(m, j) =
{
  my(ratio = batch[j][m] / phase[m], rat = bestappr(real(ratio)));
  max_imag_exponent = max(max_imag_exponent, exponent(imag(ratio)));
  max_residual_exponent = max(max_residual_exponent, exponent(ratio - rat));
  max_numerator_exponent = max(max_numerator_exponent, exponent(numerator(rat)));
  max_denominator_exponent = max(max_denominator_exponent, exponent(denominator(rat)));
  return(rat);
}
Crat = matrix(16, DK, m, j, recover_ratio(m, j));
print("MAX_IMAG_EXPONENT=", max_imag_exponent);
print("MAX_RESIDUAL_EXPONENT=", max_residual_exponent);
print("MAX_NUMERATOR_EXPONENT=", max_numerator_exponent);
print("MAX_DENOMINATOR_EXPONENT=", max_denominator_exponent);
if (max_imag_exponent > -400 || max_residual_exponent > -400, error("insufficient phase/rational reconstruction clearance", [max_imag_exponent, max_residual_exponent]));

Minf = mfcoefs(mf, 32) * K;
Chalf = Crat[1..15,];
Ctarget = Crat[16..16,];
Ahalf = matconcat([Minf~, Chalf~])~;
rank_inf = matrank(Minf);
rank_half = matrank(Chalf);
rank_combined = matrank(Ahalf);
target_variation = matrank(Ctarget * matker(Ahalf));
print("INFINITY_RANK=", rank_inf);
print("CUSP_HALF_RATIONAL_RANK=", rank_half);
print("INFINITY_HALF_RANK=", rank_combined);
print("INFINITY_HALF_DIMENSION=", DK - rank_combined);
print("TARGET_VARIATION_RANK=", target_variation);

cache_path = Str("/home/nick/p11_half_rational_", real_precision, ".gpbin");
writebin(cache_path, [Minf, Crat, base, params]);
print("CACHE_WRITTEN=", cache_path);
quit;
