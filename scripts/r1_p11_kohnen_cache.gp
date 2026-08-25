\\ Build and cache the expensive p=11 half-integral/Kohnen basis once.
default(parisize, 12G);
mf = mfinit([44, 69/2, 44], 1);
K = mfkohnenbasis(mf);
print("FULL_DIM=", mfdim(mf));
print("FULL_STURM=", mfsturm(mf));
print("KOHNEN_DIM=", matsize(K)[2]);
writebin("/home/nick/p11_mf_k.gpbin", [mf, K]);
print("CACHE_WRITTEN=/home/nick/p11_mf_k.gpbin");
quit;
