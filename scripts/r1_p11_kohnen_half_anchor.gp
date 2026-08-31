\\ Work around PARI's high-shift mfslashexpansion failure by evaluating a
\\ generic anchored basis at cusp 1/2 and subtracting the anchor exactly.
default(parisize, 12G);
default(realprecision, 200);
X = read("/home/nick/p11_mf_k.gpbin")[1];
mf = X[1];
K = X[2];
DK = matsize(K)[2];
Minf = mfcoefs(mf, 32) * K;
print("KOHNEN_DIM=", DK);
print("INFINITY_RANK=", matrank(Minf));

gap = 15;
degree = 10;
anchor_coordinates = vector(DK, j, j)~;
anchor_scale = 1000003;
Fanchor = mflinear(mf, K * anchor_coordinates);
params = 0;
c_anchor = mfslashexpansion(mf, Fanchor, [1,0;2,1], gap, 1, &params);
print("CUSP_HALF_PARAMS=", params);
print("ANCHOR_FIRST_0_15=", c_anchor);

Chalf = matrix(gap * degree, DK);
Ctarget = matrix(degree, DK);
for (j = 1, DK, {
  my(coords = anchor_scale * anchor_coordinates, F, c, recovered);
  coords[j] = coords[j] + 1;
  F = mflinear(mf, K * coords);
  c = mfslashexpansion(mf, F, [1,0;2,1], gap, 1);
  recovered = c - anchor_scale * c_anchor;
  for (m = 0, gap - 1, for (r = 0, degree - 1, Chalf[m * degree + r + 1,j] = polcoef(lift(recovered[m + 1]), r)));
  for (r = 0, degree - 1, Ctarget[r + 1,j] = polcoef(lift(recovered[gap + 1]), r));
  if (j % 10 == 0 || j == DK, print("HALF_COLUMNS_DONE=", j));
});
Ahalf = matconcat([Minf~, Chalf~])~;
print("CUSP_HALF_RANK=", matrank(Chalf));
print("INFINITY_HALF_RANK=", matrank(Ahalf));
print("INFINITY_HALF_DIMENSION=", DK - matrank(Ahalf));
print("TARGET_VARIATION_RANK=", matrank(Ctarget * matker(Ahalf)));
writebin("/home/nick/p11_half_blocks.gpbin", [Minf, Chalf, Ctarget]);
print("BLOCK_CACHE_WRITTEN=/home/nick/p11_half_blocks.gpbin");
quit;
