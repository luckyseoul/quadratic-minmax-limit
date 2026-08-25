\\ Exact p=11 modular rank audit after correcting the ordinary-dual shell.
\\
\\ The old laboratory probe constrained every infinity coefficient through
\\ 32 to vanish in the harmonic variation.  Complete qfminim enumeration
\\ now shows that coefficient 20 has a square-circle anisotropic term, and
\\ coefficients above 20 have not been proved complete.  This script keeps
\\ only the rigorous infinity rows 0..19 and tests whether the half-cusp
\\ target is controlled by the newly identified coefficient 20.
default(parisize, 12G);

H = read("/home/nick/p11_half_rational_250.gpbin");
Z = read("/home/nick/p11_cusp_zero_exact_250.gpbin");
Q = read("/home/nick/p11_cusp_quarter_exact_250.gpbin");
R = read("/home/nick/p11_cusp_p_exact_250.gpbin");

Minf = H[1];
Crat = H[2];
Chalf = Crat[1..15,];
Ctarget = Crat[16..16,];
Czero = Z[1];
Cquarter = Q[1];
Cp = R[1];

\\ GP row 1 is q^0, so rows 1..20 are coefficients 0..19 and row 21
\\ is the exact second dual shell at scaled norm 20.
Minf_pre_second = Minf[1..20,];
Csecond = Minf[21..21,];

stack(A, B) = matconcat([A~, B~])~;

report(label, A) =
{
  my(K = matker(A), Vtarget = Ctarget * K, Vsecond = Csecond * K,
     Vjoint = stack(Csecond, Ctarget) * K, alpha = 0, relation = 0,
     second_augmented = stack(A, Csecond));
  print(label, "_RANK=", matrank(A));
  print(label, "_RESIDUAL_DIMENSION=", matsize(A)[2] - matrank(A));
  print(label, "_TARGET_VARIATION_RANK=", matrank(Vtarget));
  print(label, "_SECOND_VARIATION_RANK=", matrank(Vsecond));
  print(label, "_JOINT_VARIATION_RANK=", matrank(Vjoint));
  print(label, "_TARGET_AFTER_SECOND_RANK=", matrank(Ctarget * matker(second_augmented)));
  if (matrank(Vjoint) == 1 && matrank(Vsecond) == 1,
    for (j = 1, matsize(Vsecond)[2],
      if (Vsecond[1,j] != 0 && !relation,
        alpha = Vtarget[1,j] / Vsecond[1,j];
        relation = (Vtarget == alpha * Vsecond)
      )
    );
  );
  print(label, "_TARGET_SECOND_PROPORTIONAL=", relation);
  if (relation, print(label, "_TARGET_OVER_SECOND=", alpha));
};

Ahalf = stack(Minf_pre_second, Chalf);
Azero = stack(Ahalf, Czero);
Aquarter = stack(Azero, Cquarter);
Ap = stack(Aquarter, Cp);

print("KOHNEN_DIM=", matsize(Minf)[2]);
print("INFINITY_PRE_SECOND_RANK=", matrank(Minf_pre_second));
print("SECOND_ROW_ADDS_RANK=", matrank(stack(Minf_pre_second, Csecond)) - matrank(Minf_pre_second));
report("HALF", Ahalf);
report("ZERO", Azero);
report("QUARTER", Aquarter);
report("P", Ap);
quit;
