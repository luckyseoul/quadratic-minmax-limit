\\ Locate the earliest ordinary-dual coefficient that controls the sole
\\ residual p=11 half-cusp target direction after only proved shell gaps.
default(parisize, 12G);

root = "/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25/";
X = read(Str(root, "p11_mf_k_reduced.gpbin"))[1];
mf = X[1];
Kbasis = X[2];
H = read(Str(root, "p11_half_rational_250.gpbin"));
Z = read(Str(root, "p11_cusp_zero_exact_250.gpbin"));

Crat = H[2];
Chalf = Crat[1..15,];
Ctarget = Crat[16..16,];
Czero = Z[1];
M = mfcoefs(mf, 250) * Kbasis;
Mpre = M[1..20,];                 \\ coefficients 0..19

stack(A, B) = matconcat([A~, B~])~;
Ahalf = stack(Mpre, Chalf);
Abase = stack(Ahalf, Czero);
K = matker(Abase);
Vtarget = Ctarget * K;

print("KOHNEN_DIM=", matsize(Kbasis)[2]);
print("STURM=", mfsturm(mf));
print("BASE_RANK=", matrank(Abase));
print("BASE_RESIDUAL_DIMENSION=", matsize(Abase)[2] - matrank(Abase));
print("BASE_TARGET_VARIATION_RANK=", matrank(Vtarget));

anchors = 0;
checkrow(s) =
{
  my(Crow = M[(s + 1)..(s + 1),], Vrow = Crow * K,
     rowrank = matrank(Vrow), joint = matrank(stack(Crow, Ctarget) * K),
     after = matrank(Ctarget * matker(stack(Abase, Crow))), alpha = 0,
     relation = 0);
  if (rowrank == 1 && joint == 1,
    for (j = 1, matsize(Vrow)[2],
      if (Vrow[1,j] != 0 && !relation,
        alpha = Vtarget[1,j] / Vrow[1,j];
        relation = (Vtarget == alpha * Vrow)
      )
    );
  );
  if (after == 0 || relation,
    anchors++;
    print("ANCHOR_S=", s, " ROW_RANK=", rowrank, " JOINT_RANK=", joint,
          " AFTER=", after, " PROPORTIONAL=", relation,
          if(relation, Str(" TARGET_OVER_ROW=", alpha), ""));
  );
}
for(s = 20, 200, checkrow(s));
print("ANCHOR_COUNT=", anchors);
quit;
