\\ Regression fixture for the guarded mf2gaexpansion relative-zero patch.
default(parisize, 4G);
default(realprecision, 200);

mf5 = mfinit([20, 21/2, 1], 1);
F5 = mfbasis(mf5)[1];
P5 = 0;
C5 = mfslashexpansion(mf5, F5, [1,0;2,1], 3, 1, &P5);
print("P5_PARAMS=", P5);
print("P5_CUSP_0_3=", C5);

mf7 = mfinit([28, 33/2, 28], 1);
F7 = mfbasis(mf7)[1];
P7 = 0;
C7 = mfslashexpansion(mf7, F7, [1,0;2,1], 6, 1, &P7);
print("P7_PARAMS=", P7);
print("P7_CUSP_0_6=", C7);
quit;
