\\ Exact W_4 action on the p=11 scalar-theta product subspace.
\\
\\ PARI's generic full-oldspace mfatkininit path reconstructs a 147 by 147
\\ matrix numerically and is unstable here.  We instead split the oldspace
\\ into the level-11 copies B(1), B(2), B(4) and the level-22 copies B(1),
\\ B(2).  Matrix identities give the level-11 action directly and reduce the
\\ level-22 action to the exact fast-path W_2 matrix at level 22.
\\
\\ Environment:
\\   P11_SCALAR_CACHE  half-gap scalar cache
\\   P11_W4_OUTPUT     optional binary output path
default(parisize, 16G);
default(realprecision, 250);

scalar_cache = getenv("P11_SCALAR_CACHE");
output_path = getenv("P11_W4_OUTPUT");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));

forms_to_matrix(V, E) =
{
  my(M = matrix(E + 1, #V));
  for(j = 1, #V, M[,j] = mfcoefs(V[j], E)~);
  M;
};

X = read(scalar_cache);
mfhalf = X[1];
K = X[2];
D = matsize(K)[2];

mffull = mfinit([44, 31, -11], 4);
mfeis = mfinit([44, 31, -11], 3);
mfnew44 = mfinit([44, 31, -11], 0);
mfnew11 = mfinit([11, 31, -11], 0);
mfnew22 = mfinit([22, 31, -11], 0);
S = mfsturm(mffull);
if(S != 183, error("unexpected integral Sturm bound", S));
if(mfdim(mfeis) != 6, error("unexpected Eisenstein dimension"));
if(mfdim(mfnew44) != 30, error("unexpected level-44 new dimension"));
if(mfdim(mfnew11) != 29, error("unexpected level-11 new dimension"));
if(mfdim(mfnew22) != 30, error("unexpected level-22 new dimension"));

\\ Fourier matrix of the 45-dimensional half-integral input, multiplied by
\\ the unary theta series by exact Cauchy product.
H = mfcoefs(mfhalf, S) * K;
P = matrix(S + 1, D);
for(n = 0, S, for(a = 0, sqrtint(n), P[n + 1,] = P[n + 1,] + if(a, 2, 1) * H[n - a^2 + 1,]));

V11 = mfbasis(mfnew11);
V22 = mfbasis(mfnew22);
Vold = concat(V11, apply(F -> mfbd(F, 2), V11));
Vold = concat(Vold, apply(F -> mfbd(F, 4), V11));
Vold = concat(Vold, V22);
Vold = concat(Vold, apply(F -> mfbd(F, 2), V22));
if(#Vold != 147, error("unexpected custom oldspace basis size"));

Beis = mfcoefs(mfeis, S);
Bold = forms_to_matrix(Vold, S);
Bnew = mfcoefs(mfnew44, S);
Bparts = matconcat([Beis, Bold, Bnew]);
if(matsize(Bparts) != [S + 1, 183], error("bad component Fourier matrix size"));
idx = matindexrank(Bparts);
RI = idx[1];
CJ = idx[2];
if(#RI != 183 || #CJ != 183, error("custom component basis lacks full rank"));
R = vecextract(Bparts, RI, CJ);
parts = R^-1 * vecextract(P, RI, Str("1..", D));
if(Bparts * parts != P, error("product decomposition failed at Sturm bound"));

\\ Exact Atkin blocks.  For level 44 and Q=4 the normalization constant is
\\ rationalized to one.  For level 22 and Q=2, PARI returns sqrt(2) times
\\ the actual W_2 action; the half-integral powers from B(1),B(2) turn this
\\ into the rational factors 2^15 and 2^-16 below.
gettime();
A_eis = mfatkininit(mfeis, 4);
t_eis = gettime();
gettime();
A_new = mfatkininit(mfnew44, 4);
t_new = gettime();
gettime();
A_22 = mfatkininit(mfnew22, 2);
t_22 = gettime();
if(A_eis[3] != 1, error("unexpected Eisenstein W4 normalization", A_eis[3]));
if(A_new[3] != 1, error("unexpected newspace W4 normalization", A_new[3]));
M22 = A_22[2];
if(type(M22[1,1]) != "t_INT" && type(M22[1,1]) != "t_FRAC", error("level-22 W2 matrix is not exact rational"));
if(M22 * M22 != -2 * matid(30), error("unexpected normalized level-22 W2 square"));

Mold = matrix(147, 147);
I29 = matid(29);
\\ Level 11: B1 -> 2^31 B4, B2 -> -B2, B4 -> 2^-31 B1.
for(j = 1, 29, Mold[58 + j,j] = 2^31);
for(j = 1, 29, Mold[29 + j,29 + j] = -1);
for(j = 1, 29, Mold[j,58 + j] = 1 / 2^31);
\\ Level 22: B1 -> -2^15 (sqrt(2) W2) B2 and
\\ B2 -> 2^-16 (sqrt(2) W2) B1.
for(i = 1, 30, for(j = 1, 30, Mold[117 + i,87 + j] = -2^15 * M22[i,j]));
for(i = 1, 30, for(j = 1, 30, Mold[87 + i,117 + j] = M22[i,j] / 2^16));
if(Mold * Mold != matid(147), error("custom oldspace W4 is not an involution"));

Mparts = matrix(183, 183);
for(i = 1, 6, for(j = 1, 6, Mparts[i,j] = A_eis[2][i,j]));
for(i = 1, 147, for(j = 1, 147, Mparts[6 + i,6 + j] = Mold[i,j]));
for(i = 1, 30, for(j = 1, 30, Mparts[153 + i,153 + j] = A_new[2][i,j]));
if(Mparts * Mparts != matid(183), error("full W4 block is not an involution"));

W4parts = Mparts * parts;
W4q = Bparts * W4parts;

print("STURM=", S);
print("PRODUCT_DIM=", matrank(P));
print("CUSTOM_COMPONENT_RANK=", matrank(Bparts));
print("W4_PRODUCT_DIM=", matrank(W4q));
print("W4_FIRST_NONCONSTANT_ROW_RANKS=", vector(16, j, matrank(vecextract(W4q, Str("2..", j + 1), Str("1..", D)))));
print("W4_ZERO_ROWS=", select(n -> W4q[n + 1,] == vector(D), [0..S]));
print("ATKIN_MS_EIS=", t_eis, " NEW44=", t_new, " NEW22_W2=", t_22);

if(type(output_path) == "t_STR", writebin(output_path, [W4q, parts, Mparts, RI, [t_eis, t_new, t_22]]));
quit;
