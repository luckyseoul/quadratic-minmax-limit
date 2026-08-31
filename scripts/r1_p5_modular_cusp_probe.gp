\\ Three-cusp rank calibration for the p=5 exceptional harmonic channel.
\\ Coefficients are 600 times the unphased exceptional dual shell scalars.
default(parisize, 2G);
default(realprecision, 200);

p = 5;
\\ mffromqf(4p G^-1) identifies the character as trivial at p=5.
mf = mfinit([4 * p, 21/2, 1], 1);
B = mfbasis(mf);
D = #B;
S = mfsturm(mf);
M = mfcoefs(mf, S);
print("DIM=", D);
print("STURM=", S);

known_bound = 13;
V = vector(S + 1);
V[5 + 1] = -5;
V[8 + 1] = 48;
V[9 + 1] = 45;
V[12 + 1] = -654;
V[13 + 1] = 6;

\\ The p=5 scaled dual norm is supported on 0 or 1 modulo 4.
rows = List();
for (n = 0, known_bound, listput(rows, n));
for (n = known_bound + 1, S, if (n % 4 == 2 || n % 4 == 3, listput(rows, n)));
rows = Vec(rows);
Ainf = matrix(#rows, D, i, j, M[rows[i] + 1, j]);
binf = vector(#rows, i, V[rows[i] + 1])~;
print("INFINITY_ROWS=", #rows);
print("INFINITY_RANK=", matrank(Ainf));
print("INFINITY_AUGMENTED_RANK=", matrank(matconcat([Ainf, binf])));

\\ Cusp 1/2: q^(1/20) sum c_m q^(m/5).  The odd coset minimum
\\ 26/40 = 1/20 + 3/5, so m=0,1,2 vanish and m=3 is the target.
gap = 3;
field_degree = 4;
Acusp = matrix(gap * field_degree, D);
fill_cusp_row(j, m, c) =
{
  for (r = 0, field_degree - 1,
    Acusp[m * field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_cusp_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [1,0;2,1], gap - 1, 1, &params);
  if (params[1] != 1/20 || params[2] != 5,
    error("unexpected cusp-half parameters", params));
  for (m = 0, gap - 1, fill_cusp_row(j, m, c));
}
for (j = 1, D, fill_cusp_column(j));
print("CUSP_HALF_ROWS=", matsize(Acusp)[1]);
print("CUSP_HALF_RANK=", matrank(Acusp));

A = matrix(#rows + matsize(Acusp)[1], D, i, j, if (i <= #rows, Ainf[i,j], Acusp[i - #rows,j]));
b = concat(binf~, vector(matsize(Acusp)[1]))~;
print("HALF_COMBINED_RANK=", matrank(A));
print("HALF_COMBINED_AUGMENTED_RANK=", matrank(matconcat([A, b])));

\\ Cusp 0: min(L)=p+1=6 gives exponent 6/(8p)=3/20,
\\ so m=0,1,2 vanish in q^(m/20).
zero_gap = 3;
zero_field_degree = 8;
Czero = matrix(zero_gap * zero_field_degree, D);
fill_zero_row(j, m, c) =
{
  for (r = 0, zero_field_degree - 1,
    Czero[m * zero_field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_zero_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [0,-1;1,0], zero_gap - 1, 1, &params);
  if (params[1] != 0 || params[2] != 20,
    error("unexpected cusp-zero parameters", params));
  for (m = 0, zero_gap - 1, fill_zero_row(j, m, c));
}
for (j = 1, D, fill_zero_column(j));
print("CUSP_ZERO_ROWS=", matsize(Czero)[1]);
print("CUSP_ZERO_RANK=", matrank(Czero));

A2 = matrix(matsize(A)[1] + matsize(Czero)[1], D, i, j, if (i <= matsize(A)[1], A[i,j], Czero[i - matsize(A)[1],j]));
b2 = concat(b~, vector(matsize(Czero)[1]))~;
print("ZERO_COMBINED_RANK=", matrank(A2));
print("ZERO_COMBINED_AUGMENTED_RANK=", matrank(matconcat([A2, b2])));

\\ Cusp 1/4: the phase i^s is a combination of the unshifted and parity
\\ branches.  The unshifted minimum gives 16*6/(8p)=12/5, hence m=12.
quarter_gap = 12;
quarter_field_degree = 4;
Cquarter = matrix(quarter_gap * quarter_field_degree, D);
fill_quarter_row(j, m, c) =
{
  for (r = 0, quarter_field_degree - 1,
    Cquarter[m * quarter_field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_quarter_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [1,0;4,1], quarter_gap - 1, 1, &params);
  if (params[1] != 0 || params[2] != 5,
    error("unexpected cusp-quarter parameters", params));
  for (m = 0, quarter_gap - 1, fill_quarter_row(j, m, c));
}
for (j = 1, D, fill_quarter_column(j));
print("CUSP_QUARTER_ROWS=", matsize(Cquarter)[1]);
print("CUSP_QUARTER_RANK=", matrank(Cquarter));

A3 = matrix(matsize(A2)[1] + matsize(Cquarter)[1], D, i, j, if (i <= matsize(A2)[1], A2[i,j], Cquarter[i - matsize(A2)[1],j]));
b3 = concat(b2~, vector(matsize(Cquarter)[1]))~;
rankA3 = matrank(A3);
rankAug3 = matrank(matconcat([A3, b3]));
print("THREECUSP_RANK=", rankA3);
print("THREECUSP_AUGMENTED_RANK=", rankAug3);
print("THREECUSP_SOLUTION_DIMENSION=", D - rankA3);

if (rankA3 == rankAug3, {
  my(x3 = matinverseimage(A3, b3), F3, params = 0, target, cp, c2p, first,
     bij, plus_coordinates, sh);
  F3 = mflinear(mf, x3);
  target = mfslashexpansion(mf, F3, [1,0;2,1], gap, 1, &params);
  print("THREECUSP_TARGET_EXACT=", target[gap + 1]);
  print("THREECUSP_TARGET_NUMERIC=", mfslashexpansion(mf, F3, [1,0;2,1], gap, 0)[gap + 1]);
  params = 0;
  cp = mfslashexpansion(mf, F3, [1,0;p,1], 30, 1, &params);
  first = 0; for (m = 0, 30, if (cp[m + 1] != 0, first = m; break));
  print("CUSP_ONE_OVER_P_PARAMS=", params);
  print("CUSP_ONE_OVER_P_FIRST_INDEX=", first);
  print("CUSP_ONE_OVER_P_FIRST_EXACT=", cp[first + 1]);
  params = 0;
  c2p = mfslashexpansion(mf, F3, [1,0;2 * p,1], 30, 1, &params);
  first = 0; for (m = 0, 30, if (c2p[m + 1] != 0, first = m; break));
  print("CUSP_ONE_OVER_2P_PARAMS=", params);
  print("CUSP_ONE_OVER_2P_FIRST_INDEX=", first);
  print("CUSP_ONE_OVER_2P_FIRST_EXACT=", c2p[first + 1]);
  bij = mfkohnenbijection(mf);
  plus_coordinates = matinverseimage(bij[3], x3);
  print("KOHNEN_PLUS_DIM=", matsize(bij[3])[2]);
  print("KOHNEN_PLUS_CONTAINS_FORM=", plus_coordinates != []);
  sh = mfshimura(mf, F3, 5);
  print("SHIMURA_PARAMS=", mfparams(sh[2]));
  print("SHIMURA_0_20=", mfcoefs(sh[2], 20));
  print("SHIMURA_NEW_DECOMPOSITION=", mftonew(sh[1], sh[2]));
});

quit;
