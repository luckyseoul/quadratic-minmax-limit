\\ Exact modular-space rank probe for the p=7 exceptional harmonic channel.
\\ Laboratory computation only: it tests whether the known infinity-cusp
\\ coefficients, global mod-4 support, and the odd-coset gap at cusp 1/2
\\ uniquely determine a form in S_{33/2}(Gamma0(28), chi_28).

default(parisize, 4G);
default(realprecision, 200);
mf = mfinit([28, 33/2, 28], 1);
B = mfbasis(mf);
D = #B;
S = mfsturm(mf);
print("DIM=", D);
print("STURM=", S);

\\ Coefficients of 1512 times the unphased exceptional dual theta series.
V = vector(S + 1);
V[7 + 1] = -7;
V[12 + 1] = -36;
V[15 + 1] = 441;
V[16 + 1] = 1274;
V[19 + 1] = 7152;
V[20 + 1] = -5172;

\\ Every coefficient through 20 is known.  Globally the quadratic form has
\\ support only at exponents n == 0 or 3 (mod 4).
rows = List();
for (n = 0, 20, listput(rows, n));
for (n = 21, S, if (n % 4 == 1 || n % 4 == 2, listput(rows, n)));
rows = Vec(rows);
M = mfcoefs(mf, S);
Ainf = matrix(#rows, D, i, j, M[rows[i] + 1, j]);
binf = vector(#rows, i, V[rows[i] + 1])~;
print("INFINITY_ROWS=", #rows);
print("INFINITY_RANK=", matrank(Ainf));

\\ At cusp 1/2 PARI returns q^(1/28) sum_m c_m q^(m/7).
\\ The odd coset begins at norm 50, exponent 50/56=25/28, hence m=6;
\\ c_0,...,c_5 vanish.  Expand each Q(zeta_7) equation in the power basis.
gap = 6;
field_degree = 6;
Acusp = matrix((gap + 1) * field_degree, D);
fill_cusp_row(j, m, c) =
{
  for (r = 0, field_degree - 1, Acusp[m * field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_cusp_column(j) =
{
  my(params = 0, c);
  params = 0;
  c = mfslashexpansion(mf, B[j], [1,0;2,1], gap, 1, &params);
  if (params[1] != 1/28 || params[2] != 7, error("unexpected cusp parameters", params));
  for (m = 0, gap, fill_cusp_row(j, m, c));
}
for (j = 1, D, fill_cusp_column(j));
print("CUSP_ROWS=", matsize(Acusp)[1]);
print("CUSP_RANK=", matrank(Acusp));

gap_rank(g) =
{
  my(Ag, bg, rg, raug);
  Ag = matrix(#rows + g * field_degree, D, i, j, if (i <= #rows, Ainf[i,j], Acusp[i - #rows,j]));
  bg = concat(binf~, vector(g * field_degree))~;
  rg = matrank(Ag);
  raug = matrank(matconcat([Ag, bg]));
  print("GAP_TERMS=", g, " RANK=", rg, " AUGMENTED_RANK=", raug, " CONSISTENT=", rg == raug);
}
for (g = 0, gap, gap_rank(g));

gap_rank_low_only(g) =
{
  my(Alow, blow, rg, raug);
  Alow = matrix(21 + g * field_degree, D, i, j, if (i <= 21, M[i,j], Acusp[i - 21,j]));
  blow = concat(vector(21, i, V[i]), vector(g * field_degree))~;
  rg = matrank(Alow);
  raug = matrank(matconcat([Alow, blow]));
  print("LOW_ONLY_GAP_TERMS=", g, " RANK=", rg, " AUGMENTED_RANK=", raug, " CONSISTENT=", rg == raug);
}
for (g = 0, gap, gap_rank_low_only(g));

A = matrix(#rows + gap * field_degree, D, i, j, if (i <= #rows, Ainf[i,j], Acusp[i - #rows,j]));
b = concat(binf~, vector(gap * field_degree))~;
rankA = matrank(A);
rankAug = matrank(matconcat([A, b]));
print("COMBINED_RANK=", rankA);
print("AUGMENTED_RANK=", rankAug);
print("SOLUTION_DIMENSION=", D - rankA);

target_report() =
{
  my(xA, KA, Ctarget, variation, target);
  xA = matinverseimage(A, b);
  KA = matker(A);
  Ctarget = Acusp[(gap * field_degree + 1)..((gap + 1) * field_degree),];
  variation = Ctarget * KA;
  target = Ctarget * xA;
  print("TARGET_CUSP_COMPONENTS=", target~);
  print("TARGET_VARIATION_RANK=", matrank(variation));
}
if (rankA == rankAug, target_report());

\\ Cusp 0 is the theta transform to the primal lattice L.  At p=7 the
\\ exact lattice minimum is 8.  Under S its exponent is 8/(8p)=1/7,
\\ i.e. m=4 in q^(m/28), so m=0,...,3 vanish.
zero_gap = 4;
zero_field_degree = 12;
Czero = matrix(zero_gap * zero_field_degree, D);
fill_zero_row(j, m, c) =
{
  for (r = 0, zero_field_degree - 1, Czero[m * zero_field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_zero_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [0,-1;1,0], zero_gap - 1, 1, &params);
  if (params[1] != 0 || params[2] != 28, error("unexpected cusp-zero parameters", params));
  for (m = 0, zero_gap - 1, fill_zero_row(j, m, c));
}
for (j = 1, D, fill_zero_column(j));
print("CUSP_ZERO_ROWS=", matsize(Czero)[1]);
print("CUSP_ZERO_RANK=", matrank(Czero));
A2 = matrix(matsize(A)[1] + matsize(Czero)[1], D, i, j, if (i <= matsize(A)[1], A[i,j], Czero[i - matsize(A)[1],j]));
b2 = concat(b~, vector(matsize(Czero)[1]))~;
rankA2 = matrank(A2);
rankAug2 = matrank(matconcat([A2, b2]));
print("MULTICUSP_RANK=", rankA2);
print("MULTICUSP_AUGMENTED_RANK=", rankAug2);
print("MULTICUSP_SOLUTION_DIMENSION=", D - rankA2);
multicusp_target_report() =
{
  my(x2, K2, Ctarget, variation, target);
  x2 = matinverseimage(A2, b2);
  K2 = matker(A2);
  Ctarget = Acusp[(gap * field_degree + 1)..((gap + 1) * field_degree),];
  variation = Ctarget * K2;
  target = Ctarget * x2;
  print("MULTICUSP_TARGET_COMPONENTS=", target~);
  print("MULTICUSP_TARGET_VARIATION_RANK=", matrank(variation));
}
if (rankA2 == rankAug2, multicusp_target_report());

\\ At cusp 1/4, the phase i^s is a linear combination of the unshifted
\\ and parity-character transforms because s is 0 or 3 modulo 4.  Its
\\ shortest possible unshifted primal contribution again has norm 8.  The
\\ c=4 scaling gives exponent 16*8/(8p)=16/7, i.e. m=16 in q^(m/7).
\\ Hence m=0,...,15 vanish.
quarter_gap = 16;
quarter_field_degree = 6;
Cquarter = matrix(quarter_gap * quarter_field_degree, D);
fill_quarter_row(j, m, c) =
{
  for (r = 0, quarter_field_degree - 1, Cquarter[m * quarter_field_degree + r + 1, j] = polcoef(lift(c[m + 1]), r));
}
fill_quarter_column(j) =
{
  my(params = 0, c);
  c = mfslashexpansion(mf, B[j], [1,0;4,1], quarter_gap - 1, 1, &params);
  if (params[1] != 0 || params[2] != 7, error("unexpected cusp-quarter parameters", params));
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
threecusp_target_report() =
{
  my(x3, K3, Ctarget, variation, target, F3, params = 0, c12, c0, c4,
     target_row, certificate, known_shells, known_weights, pos, cp, c2p,
     first);
  x3 = matinverseimage(A3, b3);
  K3 = matker(A3);
  Ctarget = Acusp[(gap * field_degree + 1)..((gap + 1) * field_degree),];
  variation = Ctarget * K3;
  target = Ctarget * x3;
  print("THREECUSP_TARGET_COMPONENTS=", target~);
  print("THREECUSP_TARGET_VARIATION_RANK=", matrank(variation));
  \\ Recover the uniquely forced form and print the first coefficients at
  \\ the three geometric cusps.  These are needed to calibrate PARI's
  \\ half-integral slash normalization against the exact primal shells.
  F3 = mflinear(mf, x3);
  params = 0;
  c12 = mfslashexpansion(mf, F3, [1,0;2,1], gap, 1, &params);
  print("THREECUSP_TARGET_EXACT=", c12[gap + 1]);
  print("THREECUSP_TARGET_NUMERIC=", mfslashexpansion(mf, F3, [1,0;2,1], gap, 0)[gap + 1]);
  params = 0;
  c0 = mfslashexpansion(mf, F3, [0,-1;1,0], zero_gap, 1, &params);
  print("THREECUSP_ZERO_FIRST_EXACT=", c0[zero_gap + 1]);
  print("THREECUSP_ZERO_FIRST_NUMERIC=", mfslashexpansion(mf, F3, [0,-1;1,0], zero_gap, 0)[zero_gap + 1]);
  params = 0;
  c4 = mfslashexpansion(mf, F3, [1,0;4,1], quarter_gap, 1, &params);
  print("THREECUSP_QUARTER_FIRST_EXACT=", c4[quarter_gap + 1]);
  print("THREECUSP_QUARTER_FIRST_NUMERIC=", mfslashexpansion(mf, F3, [1,0;4,1], quarter_gap, 0)[quarter_gap + 1]);
  \\ Exact dual witness: A3^T certificate is the first power-basis
  \\ component of the target cusp coefficient.  Since every non-shell
  \\ right-hand side is zero, only these six entries contribute.
  target_row = Ctarget[1,];
  certificate = matinverseimage(A3~, target_row~);
  if (certificate == [] || A3~ * certificate != target_row~,
    error("failed to construct target dual certificate"));
  known_shells = [7, 12, 15, 16, 19, 20];
  known_weights = vector(#known_shells);
  for (h = 1, #known_shells,
    pos = 0;
    for (i = 1, #rows, if (rows[i] == known_shells[h], pos = i; break));
    if (pos == 0, error("known shell missing from infinity rows"));
    known_weights[h] = certificate[pos];
  );
  print("THREECUSP_CERTIFICATE_NONZERO_ROWS=", sum(i = 1, #certificate, certificate[i] != 0));
  print("THREECUSP_KNOWN_SHELLS=", known_shells);
  print("THREECUSP_KNOWN_WEIGHTS=", known_weights);
  print("THREECUSP_CERTIFIED_TARGET_COMPONENT=", sum(h = 1, #known_shells, known_weights[h] * V[known_shells[h] + 1]));
  \\ Inspect the two remaining cusp classes for the all-p extension.
  params = 0;
  cp = mfslashexpansion(mf, F3, [1,0;7,1], 30, 1, &params);
  first = 0; for (m = 0, 30, if (cp[m + 1] != 0, first = m; break));
  print("CUSP_ONE_OVER_P_PARAMS=", params);
  print("CUSP_ONE_OVER_P_FIRST_INDEX=", first);
  print("CUSP_ONE_OVER_P_FIRST_EXACT=", cp[first + 1]);
  params = 0;
  c2p = mfslashexpansion(mf, F3, [1,0;14,1], 30, 1, &params);
  first = 0; for (m = 0, 30, if (c2p[m + 1] != 0, first = m; break));
  print("CUSP_ONE_OVER_2P_PARAMS=", params);
  print("CUSP_ONE_OVER_2P_FIRST_INDEX=", first);
  print("CUSP_ONE_OVER_2P_FIRST_EXACT=", c2p[first + 1]);
}
if (rankA3 == rankAug3, threecusp_target_report());

first_nonzero(v) =
{
  for (j = 1, #v, if (v[j] != 0, return(j)));
  return(0);
}
scan_single_coefficient_relations() =
{
  my(xA, KA, Ctarget, tv, nv, j0, factor, target0, coeff0);
  xA = matinverseimage(A, b);
  KA = matker(A);
  Ctarget = Acusp[(gap * field_degree + 1)..((gap + 1) * field_degree),];
  tv = (Ctarget * KA)[1,];
  target0 = (Ctarget * xA)[1];
  for (n = 21, S,
    if (n % 4 == 0 || n % 4 == 3,
      nv = M[n + 1,] * KA;
      j0 = first_nonzero(nv);
      if (j0 != 0,
        factor = tv[j0] / nv[j0];
        if (tv == factor * nv,
          coeff0 = M[n + 1,] * xA;
          print("TARGET_SINGLE_COEFF_RELATION n=", n,
            " target0=", target0, " coeff0=", coeff0,
            " factor=", factor,
            " law: target=target0+factor*(a_n-coeff0)")
        )
      )
    )
  );
}
if (0 && rankA == rankAug, scan_single_coefficient_relations());

scan_prefix_relation() =
{
  my(xA, KA, Ctarget, tv, allowed, R, Rplus, rankR, rankPlus, coeffs, base);
  xA = matinverseimage(A, b);
  KA = matker(A);
  Ctarget = Acusp[(gap * field_degree + 1)..((gap + 1) * field_degree),];
  tv = (Ctarget * KA)[1,];
  allowed = select(n -> n % 4 == 0 || n % 4 == 3, [21..S]);
  for (h = 1, #allowed,
    R = matrix(h, matsize(KA)[2], i, j, (M[allowed[i] + 1,] * KA)[j]);
    Rplus = matrix(h + 1, matsize(KA)[2], i, j, if (i <= h, R[i,j], tv[j]));
    rankR = matrank(R);
    rankPlus = matrank(Rplus);
    if (rankR == rankPlus,
      coeffs = matinverseimage(R~, tv~);
      base = vector(h, i, M[allowed[i] + 1,] * xA);
      print("TARGET_PREFIX_RELATION_CUTOFF=", allowed[h], " SHELLS=", allowed[1..h]);
      print("TARGET_PREFIX_RELATION_WEIGHTS=", coeffs~);
      print("TARGET_PREFIX_RELATION_BASE=", base);
      print("TARGET_PREFIX_RELATION_TARGET0=", (Ctarget * xA)[1]);
      return();
    )
  );
  print("TARGET_PREFIX_RELATION_NOT_FOUND");
}
if (0 && rankA == rankAug, scan_prefix_relation());

recover_form() =
{
  my(x, F, params = 0, cusp);
  x = matinverseimage(A, b);
  if (x == [] || A * x != b, error("failed to recover a preimage"));
  F = mflinear(mf, x);
  print("RECOVERED_INFINITY_0_40=", mfcoefs(F, 40));
  params = 0;
  cusp = mfslashexpansion(mf, F, [1,0;2,1], 12, 1, &params);
  print("RECOVERED_CUSP_PARAMS=", params);
  print("RECOVERED_CUSP_0_12=", cusp);
}
if (0 && rankA == rankAug, recover_form());

\\ Leave the five later nonzero shell coefficients free.  This displays the
\\ affine coefficient space forced by modularity, support, the first shell,
\\ and the full six-term odd-coset gap.
unknown_shells = [12, 15, 16, 19, 20];
is_unknown_shell(n) = sum(i = 1, #unknown_shells, n == unknown_shells[i]) > 0;
rows_rel = List();
for (i = 1, #rows, if (!is_unknown_shell(rows[i]), listput(rows_rel, rows[i])));
rows_rel = Vec(rows_rel);
Ainf_rel = matrix(#rows_rel, D, i, j, M[rows_rel[i] + 1, j]);
binf_rel = vector(#rows_rel, i, V[rows_rel[i] + 1])~;
Arel = matrix(#rows_rel + gap * field_degree, D, i, j, if (i <= #rows_rel, Ainf_rel[i,j], Acusp[i - #rows_rel,j]));
brel = concat(binf_rel~, vector(gap * field_degree))~;
rankRel = matrank(Arel);
rankRelAug = matrank(matconcat([Arel, brel]));
print("RELAXED_RANK=", rankRel);
print("RELAXED_AUGMENTED_RANK=", rankRelAug);
print("RELAXED_SOLUTION_DIMENSION=", D - rankRel);
relaxed_report() =
{
  my(xrel, Krel, coeff_affine);
  xrel = matinverseimage(Arel, brel);
  Krel = matker(Arel);
  coeff_affine = matrix(#unknown_shells, 1 + matsize(Krel)[2], i, j, if (j == 1, M[unknown_shells[i] + 1,] * xrel, M[unknown_shells[i] + 1,] * Krel[,j - 1]));
  print("UNKNOWN_SHELL_AFFINE_COLUMNS=", coeff_affine);
}
if (0 && rankRel == rankRelAug, relaxed_report());

quit;
