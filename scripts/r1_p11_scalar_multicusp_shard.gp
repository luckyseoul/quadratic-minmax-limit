\\ Exact per-column scalar transformations at the three remaining cusps.
\\ Environment: P11_START, P11_END, P11_SCALAR_CACHE, P11_CUSP_OUTPUT_DIRECTORY.
default(parisize, 8G);
default(realprecision, 250);
start_column = eval(getenv("P11_START"));
end_column = eval(getenv("P11_END"));
scalar_cache = getenv("P11_SCALAR_CACHE");
output_directory = getenv("P11_CUSP_OUTPUT_DIRECTORY");
if(type(start_column) != "t_INT", error("invalid P11_START"));
if(type(end_column) != "t_INT", error("invalid P11_END"));
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(output_directory) != "t_STR", error("missing P11_CUSP_OUTPUT_DIRECTORY"));
X = read(scalar_cache);
mf = X[1];
K = X[2];
if(start_column < 1, error("P11_START is below one"));
if(end_column > matsize(K)[2], error("P11_END exceeds the scalar basis dimension"));
if(start_column > end_column, error("empty scalar multicusp shard range"));

cyclotomic_column(c, count, degree) = vector(count * degree, i, polcoef(lift(c[((i - 1) \ degree) + 1]), (i - 1) % degree));

for(j = start_column, end_column, {
  F = mflinear(mf, K[,j]);
  params_zero = 0;
  params_quarter = 0;
  params_p = 0;
  gettime();
  c_zero = mfslashexpansion(mf, F, [0,-1;1,0], 5, 1, &params_zero);
  zero_ms = gettime();
  gettime();
  c_quarter = mfslashexpansion(mf, F, [1,0;4,1], 23, 1, &params_quarter);
  quarter_ms = gettime();
  gettime();
  c_p = mfslashexpansion(mf, F, [1,0;11,1], 4, 1, &params_p);
  p_ms = gettime();
  if(params_zero != [0,44,[1,0;0,1]], error("unexpected cusp-zero parameters"));
  if(params_quarter != [0,11,[1,0;0,1]], error("unexpected cusp-quarter parameters"));
  if(params_p != [0,4,[1,0;0,1]], error("unexpected cusp-p parameters"));
  column_zero = cyclotomic_column(c_zero, 6, 20);
  column_quarter = cyclotomic_column(c_quarter, 24, 10);
  column_p = cyclotomic_column(c_p, 5, 2);
  path = Str(output_directory, "/column_", j, ".gpbin");
  writebin(path, [j, column_zero, column_quarter, column_p]);
  print("COLUMN=", j, " ZERO_MS=", zero_ms, " QUARTER_MS=", quarter_ms, " P_MS=", p_ms, " WROTE=", path);
});
quit;
