\\ Extend the exact p=11 W4-gap-reduced scalar q-row export.
\\
\\ Environment:
\\   P11_SCALAR_W4_CACHE   reduced scalar cache [mf,K,pivots,E,target]
\\   P11_OUTPUT_DIRECTORY  fresh output directory
\\   P11_EXPONENT_LIMIT    largest infinity exponent to export
default(parisize, 8G);

scalar_cache = getenv("P11_SCALAR_W4_CACHE");
output_directory = getenv("P11_OUTPUT_DIRECTORY");
exponent_text = getenv("P11_EXPONENT_LIMIT");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_W4_CACHE"));
if(type(output_directory) != "t_STR", error("missing P11_OUTPUT_DIRECTORY"));
if(type(exponent_text) != "t_STR", error("missing P11_EXPONENT_LIMIT"));
E = eval(exponent_text);
if(type(E) != "t_INT" || E < 1, error("invalid P11_EXPONENT_LIMIT"));

X = read(scalar_cache);
mf = X[1];
K = X[2];
Ctarget = X[5];
gettime();
Q = mfcoefs(mf, E) * K;
coefficient_ms = gettime();
path = Str(output_directory, "/p11_scalar_w4gap_qrows_exact_e", E, "_20260828.txt");
target_path = Str(output_directory, "/p11_scalar_w4gap_target_exact_e", E, "_20260828.txt");
for(i = 1, matsize(Q)[1], write(path, Q[i,]));
for(i = 1, matsize(Ctarget)[1], write(target_path, Ctarget[i,]));
print("EXPONENT_LIMIT=", E);
print("REDUCED_DIMENSION=", matsize(K)[2]);
print("COEFFICIENT_MS=", coefficient_ms);
print("QROWS_WRITTEN=", path);
print("TARGET_ROWS_WRITTEN=", target_path);
quit;
