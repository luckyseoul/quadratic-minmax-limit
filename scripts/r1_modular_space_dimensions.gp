\\ Dimension/Sturm census for the scalar harmonic-theta ambient spaces.
default(parisize, 1G);
forprime(p = 3, 19, {
  my(k = (p^2 + 17) / 4, dim, sturm_upper);
  dim = mfdim([4 * p, k, 4 * p], 1);
  sturm_upper = mfsturm([4 * p, k]);
  print("P=", p, " WEIGHT=", k, " DIM=", dim,
    " GENERIC_STURM_UPPER=", sturm_upper);
});
quit;
