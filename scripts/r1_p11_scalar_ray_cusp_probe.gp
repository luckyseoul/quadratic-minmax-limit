\\ Test one exact scalar recession direction against every geometric cusp gap.
default(parisize, 8G);
default(realprecision, 250);
scalar_cache = getenv("P11_SCALAR_CACHE");
coordinate_path = getenv("P11_SCALAR_RAY_COORDINATES");
if(type(scalar_cache) != "t_STR", error("missing P11_SCALAR_CACHE"));
if(type(coordinate_path) != "t_STR", error("missing P11_SCALAR_RAY_COORDINATES"));
X = read(scalar_cache);
mf = X[1];
K = X[2];
v = read(coordinate_path);
if(#v != matsize(K)[2], error("ray coordinate width mismatch"));
F = mflinear(mf, K * v~);

probe(label, slash_matrix, last_index, expected_alpha, expected_width) =
{
  my(params = 0, c, nonzero);
  gettime();
  c = mfslashexpansion(mf, F, slash_matrix, last_index, 1, &params);
  print("CUSP=", label, " MS=", gettime(), " PARAMS=", params);
  if(params[1] != expected_alpha || params[2] != expected_width,
    error("unexpected cusp parameters", [label, params]));
  nonzero = select(i -> c[i + 1] != 0, [0..last_index]);
  print("CUSP=", label, " NONZERO_INDICES=", nonzero);
  print("CUSP=", label, " COEFFICIENTS=", c);
}

\\ Include index 15 at the half cusp: the normalized ray must preserve both
\\ its fifteen zero coefficients and the fixed first-shell coefficient.
probe("half", [1,0;2,1], 15, 1/44, 11);
probe("zero", [0,-1;1,0], 5, 0, 44);
probe("quarter", [1,0;4,1], 23, 0, 11);
probe("p", [1,0;11,1], 4, 0, 4);
quit;
