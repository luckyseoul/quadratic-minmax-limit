\\ Exact rank test for the p=11 exceptional theta form inside its Kohnen
\\ subspace.  No shell values are needed: this asks whether coefficients
\\ through 3p-1 and the five geometric cusp gaps force the target at 1/2.
default(parisize, 12G);
default(realprecision, 200);

p = 11;
mf = mfinit([44, 69/2, 44], 1);
print("FULL_DIM=", mfdim(mf));
print("FULL_STURM=", mfsturm(mf));
K = mfkohnenbasis(mf);
DK = matsize(K)[2];
print("KOHNEN_DIM=", DK);

\\ The form is known at infinity through scaled norm 3p-1=32.  Work in
\\ Kohnen coordinates, where the global 0,3 mod-4 support is automatic.
known_bound = 32;
Minf = mfcoefs(mf, known_bound) * K;
print("INFINITY_RANK=", matrank(Minf));

\\ Convert the Kohnen coordinate columns to modular forms once.
BK = vector(DK, j, mflinear(mf, K[,j]));

\\ Cusp 1/2: alpha=1/44, width 11.  The odd-coset minimum gives target
\\ index (p^2-1)/8=15, so indices 0,...,14 vanish.
half_gap = 15;
half_degree = 10;
Chalf = matrix(half_gap * half_degree, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [1,0;2,1], half_gap - 1, 1, &params));
  if (j == 1, print("CUSP_HALF_PARAMS=", params));
  for (m = 0, half_gap - 1, for (r = 0, half_degree - 1, Chalf[m * half_degree + r + 1,j] = polcoef(lift(c[m + 1]), r)));
});
print("CUSP_HALF_RANK=", matrank(Chalf));
Ahalf = matconcat([Minf~, Chalf~])~;
print("INFINITY_HALF_RANK=", matrank(Ahalf));

\\ Cusp 0: min(L)=p+1=12 gives index (p+1)/2=6.
zero_gap = 6;
zero_degree = 20;
Czero = matrix(zero_gap * zero_degree, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [0,-1;1,0], zero_gap - 1, 1, &params));
  if (j == 1, print("CUSP_ZERO_PARAMS=", params));
  for (m = 0, zero_gap - 1, for (r = 0, zero_degree - 1, Czero[m * zero_degree + r + 1,j] = polcoef(lift(c[m + 1]), r)));
});
print("CUSP_ZERO_RANK=", matrank(Czero));
Azero = matconcat([Ahalf~, Czero~])~;
print("INFINITY_HALF_ZERO_RANK=", matrank(Azero));

\\ Cusp 1/4: the unshifted L minimum gives index 2(p+1)=24.
quarter_gap = 24;
quarter_degree = 10;
Cquarter = matrix(quarter_gap * quarter_degree, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [1,0;4,1], quarter_gap - 1, 1, &params));
  if (j == 1, print("CUSP_QUARTER_PARAMS=", params));
  for (m = 0, quarter_gap - 1, for (r = 0, quarter_degree - 1, Cquarter[m * quarter_degree + r + 1,j] = polcoef(lift(c[m + 1]), r)));
});
print("CUSP_QUARTER_RANK=", matrank(Cquarter));
Aquarter = matconcat([Azero~, Cquarter~])~;
print("THREECUSP_RANK=", matrank(Aquarter));

\\ Cusp 1/p: width 4, first possible index (p-1)/2=5.
p_gap = 5;
p_degree = 2;
Cp = matrix(p_gap * p_degree, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [1,0;p,1], p_gap - 1, 1, &params));
  if (j == 1, print("CUSP_P_PARAMS=", params));
  for (m = 0, p_gap - 1, for (r = 0, p_degree - 1, Cp[m * p_degree + r + 1,j] = polcoef(lift(c[m + 1]), r)));
});
print("CUSP_P_RANK=", matrank(Cp));
Ap = matconcat([Aquarter~, Cp~])~;
print("FOURCUSP_RANK=", matrank(Ap));

\\ Cusp 1/(2p): width 1, one geometrically forced zero term.
C2p = matrix(1, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [1,0;2 * p,1], 0, 1, &params));
  if (j == 1, print("CUSP_2P_PARAMS=", params));
  C2p[1,j] = polcoef(lift(c[1]), 0);
});
print("CUSP_2P_RANK=", matrank(C2p));
Aall = matconcat([Ap~, C2p~])~;
rank_all = matrank(Aall);
print("ALL_CONSTRAINT_RANK=", rank_all);
print("SOLUTION_DIMENSION=", DK - rank_all);

\\ Does the first nonzero odd-coset coefficient vary in the residual kernel?
Ctarget = matrix(half_degree, DK);
for (j = 1, DK, {
  my(params = 0, c = mfslashexpansion(mf, BK[j], [1,0;2,1], half_gap, 1, &params));
  for (r = 0, half_degree - 1, Ctarget[r + 1,j] = polcoef(lift(c[half_gap + 1]), r));
});
print("TARGET_VARIATION_RANK=", matrank(Ctarget * matker(Aall)));
quit;
