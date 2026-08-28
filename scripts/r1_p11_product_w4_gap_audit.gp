\\ Audit the exact W_4 product rows against the geometric cusp-1/11 gap.
default(parisize, 4G);

scalar_cache = getenv("P11_SCALAR_CACHE");
w4_cache = getenv("P11_W4_CACHE");
multicusp_cache = getenv("P11_MULTICUSP_CACHE");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(w4_cache) != "t_STR", error("missing P11_W4_CACHE"));

X = read(scalar_cache);
mf = X[1];
K = X[2];
D = matsize(K)[2];
W = read(w4_cache);
print("W4_CACHE_TYPE=", type(W), " LENGTH=", #W, " FIRST_TYPE=", type(W[1]));
if(type(W[1]) == "t_VEC" && #W[1] == 5 && type(W[1][1]) == "t_MAT", W = W[#W]);
W4q = W[1];
if(matsize(W4q) != [184, D], error("unexpected W4 q-row matrix"));

\\ If the half-integral theta series has no terms q^1,...,q^4 at this
\\ cusp, multiplication by theta gives coefficients c*[1,2q,0q^2,0q^3,
\\ 2q^4] through q^4.  Eliminate the unknown constant c.
G = matrix(4, D);
G[1,] = W4q[2,] - 2 * W4q[1,];
G[2,] = W4q[3,];
G[3,] = W4q[4,];
G[4,] = W4q[5,] - 2 * W4q[1,];

Ainf = mfcoefs(mf, 28) * K;
known_indices = [0, 11, 20, 24, 27];
Aknown = matrix(#known_indices, D, i, j, Ainf[known_indices[i] + 1,j]);
bknown = [1, 244, 16104, 14762, 442860]~;
A = matconcat([Aknown~, G~])~;
b = concat(bknown~, vector(4))~;

print("EXACT_W4_GAP_RANK=", matrank(G));
print("KNOWN_RANK=", matrank(Aknown));
print("KNOWN_PLUS_W4_GAP_RANK=", matrank(A));
print("KNOWN_PLUS_W4_GAP_AUGMENTED_RANK=", matrank(matconcat([A, b])));
print("AFFINE_DIMENSION=", D - matrank(A));

if(type(multicusp_cache) == "t_STR", {
  Y = read(multicusp_cache);
  Gp = Y[6];
  print("OLD_RECONSTRUCTED_P_GAP_RANK=", matrank(Gp));
  print("EXACT_PLUS_OLD_P_GAP_RANK=", matrank(matconcat([G~, Gp~])~));
  print("KNOWN_PLUS_OLD_P_GAP_RANK=", matrank(matconcat([Aknown~, Gp~])~));
  print("KNOWN_PLUS_OLD_P_GAP_AUGMENTED_RANK=", matrank(matconcat([matconcat([Aknown~, Gp~])~, concat(bknown~, vector(matsize(Gp)[1]))~])));
});
quit;
