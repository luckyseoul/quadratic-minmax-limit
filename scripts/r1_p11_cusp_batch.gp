\\ Batch one p=11 cusp over the reduced Kohnen basis.  The custom PARI
\\ batch path reuses the expensive transformed Eisenstein basis, then
\\ applies PARI's exact cyclotomic reconstruction to every form.
default(parisize, 12G);
label = getenv("P11_CUSP");
precision_text = getenv("P11_REALPRECISION");
real_precision = if(precision_text, eval(precision_text), 250);
default(realprecision, real_precision);

slash_matrix = 0;
gap = 0;
if (label == "zero", slash_matrix = [0,-1;1,0]; gap = 6; field_degree = 20; expected_alpha = 0; expected_width = 44);
if (label == "quarter", slash_matrix = [1,0;4,1]; gap = 24; field_degree = 10; expected_alpha = 0; expected_width = 11);
if (label == "p", slash_matrix = [1,0;11,1]; gap = 5; field_degree = 2; expected_alpha = 0; expected_width = 4);
if (label == "2p", slash_matrix = [1,0;22,1]; gap = 1; field_degree = 2; expected_alpha = 3/4; expected_width = 1);
if (type(slash_matrix) != "t_MAT", error("unknown P11_CUSP label", label));

X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
mf = X[1];
K = X[2];
DK = matsize(K)[2];
forms = vector(DK, j, mflinear(mf, K[,j]));

gettime();
params = 0;
batch = mfslashexpansion(mf, forms, slash_matrix, gap - 1, 1, &params);
batch_time = gettime();
if (type(batch) != "t_VEC", print("BATCH_FAILED"); quit(1));
print("CUSP_LABEL=", label);
print("REAL_PRECISION=", real_precision);
print("KOHNEN_DIM=", DK);
print("CUSP_PARAMS=", params);
print("BATCH_TIME_MS=", batch_time);
print("FIELD_DEGREE=", field_degree);
if (params[1] != expected_alpha || params[2] != expected_width, error("unexpected cusp parameters", params));

raw_cache_path = Str("/home/nick/p11_cusp_", label, "_exact_batch_", real_precision, ".gpbin");
writebin(raw_cache_path, [batch, params]);
print("RAW_CACHE_WRITTEN=", raw_cache_path);
Cexact = matrix(gap * field_degree, DK, i, j, polcoef(lift(batch[j][((i - 1) \ field_degree) + 1]), (i - 1) % field_degree));
print("EXACT_MATRIX_ROWS=", matsize(Cexact)[1]);
print("EXACT_ROW_RANK=", matrank(Cexact));

cache_path = Str("/home/nick/p11_cusp_", label, "_exact_", real_precision, ".gpbin");
writebin(cache_path, [Cexact, params]);
print("CACHE_WRITTEN=", cache_path);
quit;
