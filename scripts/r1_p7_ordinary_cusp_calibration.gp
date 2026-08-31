\\ Control for r1_p7_modular_cusp_probe.gp.  This uses the ordinary p=7
\\ dual theta series (no harmonic polynomial and no exceptional block).

default(parisize, 4G);
default(realprecision, 200);
mf = mfinit([28, 25/2, 28], 4);
B = mfbasis(mf);
D = #B;
S = mfsturm(mf);
print("DIM=", D);
print("STURM=", S);

V = vector(S + 1);
V[1] = 1;
V[7 + 1] = 100;
V[12 + 1] = 2800;
V[15 + 1] = 36260;
V[16 + 1] = 63700;
V[19 + 1] = 436800;
V[20 + 1] = 825720;

rows = List();
for (n = 0, 20, listput(rows, n));
for (n = 21, S, if (n % 4 == 1 || n % 4 == 2, listput(rows, n)));
rows = Vec(rows);
M = mfcoefs(mf, S);
Ainf = matrix(#rows, D, i, j, M[rows[i] + 1, j]);
binf = vector(#rows, i, V[rows[i] + 1])~;
print("INFINITY_ROWS=", #rows);
print("INFINITY_RANK=", matrank(Ainf));
print("INFINITY_AUGMENTED_RANK=", matrank(matconcat([Ainf, binf])));

gap = 6;
field_degree = 6;
Acusp = matrix(gap * field_degree, D);
fill_cusp_row(j, m, c) =
{
  for (r = 0, field_degree - 1, Acusp[m * field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_cusp_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [1,0;2,1], gap - 1, 1, &params);
  if (params[1] != 1/28 || params[2] != 7, error("unexpected cusp parameters", params));
  for (m = 0, gap - 1, fill_cusp_row(j, m, c));
}
for (j = 1, D, fill_cusp_column(j));
print("CUSP_RANK=", matrank(Acusp));

gap_rank(g) =
{
  my(A, b, rg, raug);
  A = matrix(#rows + g * field_degree, D, i, j, if (i <= #rows, Ainf[i,j], Acusp[i - #rows,j]));
  b = concat(binf~, vector(g * field_degree))~;
  rg = matrank(A);
  raug = matrank(matconcat([A, b]));
  print("GAP_TERMS=", g, " RANK=", rg, " AUGMENTED_RANK=", raug, " CONSISTENT=", rg == raug);
}
for (g = 0, gap, gap_rank(g));

quit;
