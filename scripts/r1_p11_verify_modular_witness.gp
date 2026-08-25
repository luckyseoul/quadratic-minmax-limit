\\ Verify the exact p=11 modular-independence payload backed up on storage.
default(parisize, 2G);
X = read("/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25/p11_modular_independence_witness.gpbin");
\\ Binary files carrying PARI variable priorities return [symbolic,data].
W = if (#X == 2 && type(X[2]) == "t_VEC" && #X[2] == 3, X[2], X);
print("WITNESS_COORDINATES=", #W[1]);
print("WITNESS_CONSTRAINT_ZERO=", W[2] == 0);
print("WITNESS_TARGET=", W[3]);
if (#W[1] != 66 || W[2] != 0 || W[3] != [1]~, error("invalid modular witness payload"));
quit;
