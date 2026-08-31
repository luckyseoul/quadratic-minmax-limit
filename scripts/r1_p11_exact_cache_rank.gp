\\ Combine the exact/reconstructed p=11 cusp caches currently available.
default(parisize, 12G);

H = read("/home/nick/p11_half_rational_250.gpbin");
Z = read("/home/nick/p11_cusp_zero_exact_250.gpbin");
Q = read("/home/nick/p11_cusp_quarter_exact_250.gpbin");
P = read("/home/nick/p11_cusp_p_exact_250.gpbin");
Minf = H[1];
Crat = H[2];
Chalf = Crat[1..15,];
Ctarget = Crat[16..16,];
Czero = Z[1];
Cquarter = Q[1];
Cp = P[1];

Ainfhalf = matconcat([Minf~, Chalf~])~;
Azero = matconcat([Ainfhalf~, Czero~])~;
Aquarter = matconcat([Azero~, Cquarter~])~;
Ap = matconcat([Aquarter~, Cp~])~;
print("INFINITY_RANK=", matrank(Minf));
print("INFINITY_HALF_RANK=", matrank(Ainfhalf));
print("ZERO_BLOCK_RANK=", matrank(Czero));
print("INFINITY_HALF_ZERO_RANK=", matrank(Azero));
print("RESIDUAL_DIMENSION=", matsize(Azero)[2] - matrank(Azero));
print("TARGET_VARIATION_RANK=", matrank(Ctarget * matker(Azero)));
print("QUARTER_BLOCK_RANK=", matrank(Cquarter));
print("INFINITY_HALF_ZERO_QUARTER_RANK=", matrank(Aquarter));
print("QUARTER_RESIDUAL_DIMENSION=", matsize(Aquarter)[2] - matrank(Aquarter));
print("QUARTER_TARGET_VARIATION_RANK=", matrank(Ctarget * matker(Aquarter)));
print("P_BLOCK_RANK=", matrank(Cp));
print("THROUGH_P_RANK=", matrank(Ap));
print("THROUGH_P_RESIDUAL_DIMENSION=", matsize(Ap)[2] - matrank(Ap));
print("THROUGH_P_TARGET_VARIATION_RANK=", matrank(Ctarget * matker(Ap)));
quit;
