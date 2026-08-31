\\ Export the exact q-pivots and half-cusp target functional needed to
\\ evaluate any p=11 modular form from its infinity coefficients.
default(parisize, 1G);
H = read("/home/nick/p11_half_rational_250.gpbin");
X = read("/home/nick/p11_mf_k_reduced.gpbin")[1];
RI = X[3];
Ctarget = H[2][16,];
path = "/home/nick/p11_trace_target_exact_v2_20260827.txt";
system(Str("rm -f -- ", path));
write(path, vector(#RI, i, RI[i] - 1));
write(path, Ctarget);
print("PIVOTS=", #RI, " TARGET_WIDTH=", #Ctarget, " WROTE=", path);
quit;
