\\ Test the one-dimensional cyclotomic phase law at cusp 1/2.  If every
\\ fixed-index coefficient is a rational multiple of the exact first
\\ reduced-basis column, rational reconstruction replaces number-field LLL.
default(parisize, 12G);
default(realprecision, 250);

X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
base = read("/home/nick/p11_half_reduced_col1_exact_250.gpbin");
z = exp(2 * Pi * I / 11);
phase = vector(16, m, subst(lift(base[m]), 't, z));
indices = [1, 2, 3, 8, 16];

report_one(j, c, h) =
{
  my(m = indices[h], ratio, rat);
  ratio = c[m] / phase[m];
  rat = bestappr(real(ratio));
  print("COLUMN=", j, " INDEX=", m - 1, " RATIONAL=", rat, " IMAG_EXPONENT=", exponent(imag(ratio)), " RESIDUAL_EXPONENT=", exponent(ratio - rat));
}
report_batch_column(j) =
{
  for (h = 1, #indices, report_one(j, batch[j], h));
  print("PHASE_COLUMNS_DONE=", j);
}
forms = vector(6, j, mflinear(mf, K[,j]));
params = 0;
batch = mfslashexpansion(mf, forms, [1,0;2,1], 15, 0, &params);
print("CUSP_HALF_PARAMS=", params);
for (j = 1, 6, report_batch_column(j));
quit;
