\\ Assemble the 60 exact scalar cusp-1/2 columns and export the gap matrix.
\\ Environment: P11_HALF_INPUT_DIRECTORY.
default(parisize, 4G);
input_directory = getenv("P11_HALF_INPUT_DIRECTORY");
if(type(input_directory) != "t_STR", error("missing P11_HALF_INPUT_DIRECTORY"));
M = matrix(160, 60);
params = 0;
for(j = 1, 60, {
  X = read(Str(input_directory, "/column_", j, ".gpbin"));
  if(X[1] != j || #X[2] != 160, error("bad scalar-half column", j));
  if(j == 1, params = X[3], if(X[3] != params, error("cusp parameter mismatch", j)));
  M[,j] = X[2]~;
});
gap = vecextract(M, "1..150", "1..60");
target = vecextract(M, "151..160", "1..60");
print("GAP_RANK=", matrank(gap));
print("TARGET_RANK=", matrank(target));
print("CUSP_HALF_PARAMS=", params);
cache = Str(input_directory, "/p11_scalar_half_gap_exact_v1_20260827.gpbin");
writebin(cache, [gap, target, params]);
path = Str(input_directory, "/p11_scalar_half_gap_exact_v1_20260827.txt");
system(Str("rm -f -- ", path));
for(i = 1, matsize(gap)[1], write(path, gap[i,]));
print("CACHE_WRITTEN=", cache);
print("GAP_ROWS_WRITTEN=", path);
quit;
