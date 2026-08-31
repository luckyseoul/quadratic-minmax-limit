\\ Export all proved non-infinity cusp-gap rows in the 66-coordinate
\\ pivot-normalized p=11 Kohnen basis.  Exact certification uses the source
\\ GP binaries; this decimal file is only for floating LP reconnaissance.
default(parisize, 2G);
default(realprecision, 80);

H = read("/home/nick/p11_half_rational_250.gpbin");
Z = read("/home/nick/p11_cusp_zero_exact_250.gpbin");
Q = read("/home/nick/p11_cusp_quarter_exact_250.gpbin");
R = read("/home/nick/p11_cusp_p_exact_250.gpbin");

stack(A, B) = matconcat([A~, B~])~;
A = H[2][1..15,];
A = stack(A, Z[1]);
A = stack(A, Q[1]);
A = stack(A, R[1]);

path = "/home/nick/p11_cusp_constraints_float_20260827.txt";
for(i = 1, matsize(A)[1], write(path, vector(matsize(A)[2], j, 1.0 * A[i,j])));
print("ROWS=", matsize(A)[1], " COLS=", matsize(A)[2], " RANK=", matrank(A));
print("WROTE=", path);
quit;
