\\ Exact per-column scalar support transformations at cusp 1/2.
\\ Environment: P11_START, P11_END, P11_SCALAR_CACHE,
\\              P11_HALF_OUTPUT_DIRECTORY.
default(parisize, 5G);
default(realprecision, 250);
start_column = eval(getenv("P11_START"));
end_column = eval(getenv("P11_END"));
output_directory = getenv("P11_HALF_OUTPUT_DIRECTORY");
cache_path = getenv("P11_SCALAR_CACHE");
if(type(start_column) != "t_INT", error("invalid P11_START"));
if(type(end_column) != "t_INT", error("invalid P11_END"));
if(type(output_directory) != "t_STR", error("missing P11_HALF_OUTPUT_DIRECTORY"));
if(type(cache_path) != "t_STR", error("missing P11_SCALAR_CACHE"));
X = read(cache_path);
mf = X[1];
K = X[2];
if(start_column < 1, error("P11_START is below one"));
if(end_column > matsize(K)[2], error("P11_END exceeds the scalar basis dimension"));
if(start_column > end_column, error("empty scalar-half shard range"));
for(j = start_column, end_column, {
  F = mflinear(mf, K[,j]);
  params = 0;
  gettime();
  c = mfslashexpansion(mf, F, [1,0;2,1], 15, 1, &params);
  elapsed = gettime();
  if(#c != 16, error("unexpected scalar-half coefficient count", [j, #c]));
  column = vector(160, i,
    polcoef(lift(c[((i - 1) \ 10) + 1]), (i - 1) % 10));
  path = Str(output_directory, "/column_", j, ".gpbin");
  writebin(path, [j, column, params]);
  print("COLUMN=", j, " MS=", elapsed, " PARAMS=", params, " WROTE=", path);
});
quit;
