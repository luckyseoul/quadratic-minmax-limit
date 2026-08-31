\\ Validate the custom fixed-cusp batch path against ordinary calls.
default(parisize, 4G);
default(realprecision, 120);

check_batch(params, n) =
{
  my(mf = mfinit(params, 1), B = mfbasis(mf), forms = [B[1], B[2]], P = 0,
     batch, c1, c2, exact_batch, e1, e2);
  batch = mfslashexpansion(mf, forms, [1,0;2,1], n, 0, &P);
  c1 = mfslashexpansion(mf, forms[1], [1,0;2,1], n, 0);
  c2 = mfslashexpansion(mf, forms[2], [1,0;2,1], n, 0);
  print("SPACE=", params, " PARAMS=", P, " BATCH_MATCH=", batch == [c1,c2]);
  exact_batch = mfslashexpansion(mf, forms, [1,0;2,1], n, 1);
  e1 = mfslashexpansion(mf, forms[1], [1,0;2,1], n, 1);
  e2 = mfslashexpansion(mf, forms[2], [1,0;2,1], n, 1);
  print("SPACE=", params, " EXACT_BATCH_MATCH=", exact_batch == [e1,e2]);
  if (batch != [c1,c2] || exact_batch != [e1,e2],
    error("batch regression mismatch", params));
}
check_batch([20, 21/2, 1], 3);
check_batch([28, 33/2, 28], 6);
quit;
