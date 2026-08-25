\\ Exact p=11 witness that the known modular/cusp constraints, including
\\ the complete second dual shell, do not determine the half-cusp target.
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
Minf_pre_second = Minf[1..20,];
Csecond = Minf[21..21,];

stack(A, B) = matconcat([A~, B~])~;
A = stack(stack(stack(stack(Minf_pre_second, Chalf), Czero), Cquarter), Cp);
Asecond = stack(A, Csecond);
K = matker(Asecond);
j0 = 0;
for (j = 1, matsize(K)[2], if (!j0 && (Ctarget * K[,j])[1] != 0, j0 = j));
if (!j0, error("target unexpectedly vanishes on the exact residual kernel"));
w = K[,j0] / (Ctarget * K[,j0])[1];

print("CONSTRAINT_RANK=", matrank(Asecond));
print("RESIDUAL_DIMENSION=", matsize(Asecond)[2] - matrank(Asecond));
print("WITNESS_CONSTRAINT_ZERO=", Asecond * w == 0);
print("WITNESS_TARGET=", (Ctarget * w)[1]);
print("WITNESS_SECOND_SHELL=", (Csecond * w)[1]);
print("WITNESS_NONZERO_COORDINATES=", sum(i = 1, #w, w[i] != 0));
print("WITNESS_MAX_NUMERATOR_EXPONENT=", vecmax(vector(#w, i, exponent(numerator(w[i])))));
print("WITNESS_MAX_DENOMINATOR_EXPONENT=", vecmax(vector(#w, i, exponent(denominator(w[i])))));
writebin("/home/nick/p11_modular_independence_witness.gpbin", [w, Asecond * w, Ctarget * w]);
print("CACHE_WRITTEN=/home/nick/p11_modular_independence_witness.gpbin");
quit;
