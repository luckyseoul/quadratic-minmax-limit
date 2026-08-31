\\ Audit candidate scalar cusp gaps against the certified coefficients 0..28.
default(parisize, 6G);
scalar_cache = getenv("P11_SCALAR_CACHE");
multicusp_cache = getenv("P11_MULTICUSP_CACHE");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(multicusp_cache) != "t_STR", error("missing P11_MULTICUSP_CACHE"));
X = read(scalar_cache);
mf = X[1];
K = X[2];
Y = read(multicusp_cache);
Gzero = Y[4];
Gquarter = Y[5];
Gp = Y[6];
Ainf = mfcoefs(mf, 28) * K;
b = vector(29);
b[1] = 1;
b[12] = 244;
b[21] = 16104;
b[25] = 14762;
b[28] = 442860;

audit(label, G) =
{
  my(A = matconcat([Ainf~, G~])~, rhs = concat(b, vector(matsize(G)[1]))~, r, ra);
  r = matrank(A);
  ra = matrank(matconcat([A, rhs]));
  print("AUDIT=", label, " ROWS=", matsize(G)[1], " RANK=", r,
        " AUGMENTED_RANK=", ra, " CONSISTENT=", r == ra,
        " AFFINE_DIMENSION=", matsize(K)[2] - r);
};

audit("infinity", matrix(0, 45));
for(k = 1, 5, audit(Str("zero_1_", k), vecextract(Gzero, Str("1..", 20 * k), "1..45")));
for(k = 1, 23, audit(Str("quarter_1_", k), vecextract(Gquarter, Str("1..", 10 * k), "1..45")));
for(k = 1, 4, audit(Str("p_1_", k), vecextract(Gp, Str("1..", 2 * k), "1..45")));
audit("zero_all", Gzero);
audit("quarter_all", Gquarter);
audit("p_all", Gp);
audit("quarter_plus_p", matconcat([Gquarter~, Gp~])~);
audit("all_three", matconcat([Gzero~, Gquarter~, Gp~])~);
quit;
