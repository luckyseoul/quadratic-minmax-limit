\\ Extend the exact p=11 harmonic channel affine q-row exports.
\\
\\ Environment:
\\   P11_MF_CACHE              exact reduced Kohnen modular-form cache
\\   P11_CHANNEL_REDUCTIONS    exact three-channel affine reductions
\\   P11_TARGET_SOURCE_PREFIX  prefix of prior exports containing target rows
\\   P11_OUTPUT_DIRECTORY      fresh output directory
\\   P11_EXPONENT_LIMIT        largest infinity exponent to export
default(parisize, 8G);

mf_cache = getenv("P11_MF_CACHE");
reductions_cache = getenv("P11_CHANNEL_REDUCTIONS");
target_prefix = getenv("P11_TARGET_SOURCE_PREFIX");
output_directory = getenv("P11_OUTPUT_DIRECTORY");
exponent_text = getenv("P11_EXPONENT_LIMIT");
if(type(mf_cache) != "t_STR", error("missing P11_MF_CACHE"));
if(type(reductions_cache) != "t_STR", error("missing P11_CHANNEL_REDUCTIONS"));
if(type(output_directory) != "t_STR", error("missing P11_OUTPUT_DIRECTORY"));
if(type(exponent_text) != "t_STR", error("missing P11_EXPONENT_LIMIT"));
E = eval(exponent_text);
if(type(E) != "t_INT" || E < 1, error("invalid P11_EXPONENT_LIMIT"));

Xmf = read(mf_cache)[1];
mf = Xmf[1];
Kmf = Xmf[2];
gettime();
Mfull = mfcoefs(mf, E) * Kmf;
coefficient_ms = gettime();
reductions = read(reductions_cache);

for(c = 1, #reductions, {
  label = reductions[c][1];
  x0 = reductions[c][2];
  K = reductions[c][3];
  base = Mfull * x0;
  M = Mfull * K;
  if(#reductions[c] >= 7,
    target = concat([reductions[c][6]], vector(matsize(M)[2], j, reductions[c][7][j])),
    if(type(target_prefix) != "t_STR", error("missing P11_TARGET_SOURCE_PREFIX for legacy reduction"));
    prior = readvec(Str(target_prefix, label, "_qrows_exact_20260827.txt"));
    target = prior[#prior]
  );
  if(#target != matsize(M)[2] + 1, error("target width mismatch", label));
  path = Str(output_directory, "/p11_affine_", label, "_qrows_exact_e", E, "_20260828.txt");
  for(i = 1, matsize(M)[1],
    write(path, concat([base[i]], vector(matsize(M)[2], j, M[i,j]))));
  write(path, target);
  print(label, " ROWS=", matsize(M)[1], " FREE=", matsize(M)[2], " WROTE=", path);
});
print("EXPONENT_LIMIT=", E);
print("COEFFICIENT_MS=", coefficient_ms);
quit;
