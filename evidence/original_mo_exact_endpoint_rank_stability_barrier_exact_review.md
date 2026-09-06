# Independent review: exact endpoint rank-stability barrier

2026-09-06. Reviewer: optimized_profile_exact.

I directly read the ENTIRE 108-line source
`/tmp/original_mo_exact_endpoint_rank_stability_barrier.md`, SHA256
`c32b3d0aac5dd1551e91cc70c1c5755134222118fc5081d2092032cb77414dd4`.
The hash was checked directly. I supplied no derivation, example,
correction, or proof step to this barrier and independently checked
every displayed claim. The note is self-contained; no external theorem,
numeric test, or unproved signing construction is imported.

## Exact rank obstruction

Modulo two every complete zero-diagonal signing is I+J. Its kernel
equation forces x=(sum x_i)1. For even n this has only the zero
solution; for odd n its solutions form the one-dimensional span of 1.
Nonzero minors modulo two remain nonzero over the reals, proving (1.1).

If the nonzero singular values of B are exactly d, WW^T and W^TW
are projections. From H^2+P_L<=I, every vector in ran(P_L) has zero
quadratic form under H^2 and hence is annihilated by H. This gives
HW=0. The right inequality similarly gives H W^T=0, whose transpose
gives WH=0. Thus AB=BA=0 even without assuming B symmetric.

The real kernel bound then gives rank B<=1 and rules out nonzero B
for even n. A complete rank-one B has Frobenius norm n and hence its
single nonzero singular value is n. Therefore d=n for odd n, and the
claimed exact positive-rank-fraction obstruction at sqrt(n) scale is
valid. It supplies no lower spectral gap away from zero for A/sqrt(n).

## Actual family and every normalization

For F=2I_4-J_4, the spectrum is -2 on the all-one vector and +2 on
its orthogonal complement. Thus H=F^{tensor ell} is symmetric with
diagonal one, H^2=tI, and both signs of sqrt(t) in its spectrum.
Consequently A=H tensor J_2-I has exactly zero diagonal and signs
elsewhere, and B=H tensor C_2 has signs in every entry.

Since J_2 C_2=C_2 J_2=0, direct multiplication gives AB=BA=-B.
On the pair-constant subspace A=2H-I and B=0; on the pair-zero-sum
subspace A=-I and B=2H. Therefore B has rank t=n/2 and nonzero
singular value 2sqrt(t). A remains invertible on both subspaces.

The actual paired K has K^2 block-diagonal with A^2+B^2. The largest
eigenvalue is (2sqrt(t)+1)^2, attained on a negative H eigenspace in
the pair-constant part. The other part has eigenvalue 4t+1. Hence
||K||=d=2sqrt(t)+1, and D=dI is indeed feasible with delta=0.

For N=4t the normalized trace is d/sqrt(N)=1+1/(2sqrt(t))<=5/4.
The full weighted cross law is precisely half zero and half at
4t/d^2, tending to the two endpoints. Meanwhile
||(A/d)(B/d)||=||B||/d^2=2sqrt(t)/d^2 tends to zero. This verifies
all the stated generic near-annihilation and cap premises while B has
macroscopic rank supported on A's eigenvalue -1 subspace.

## Original norms and strict scope

The displayed Boolean eigenvectors of F have eigenvalues of opposite
sign. Taking one negative tensor factor and positive factors elsewhere
produces an actual Boolean H eigenvector of eigenvalue -sqrt(t).
Its tensor with the all-one pair vector is an A eigenvector of
eigenvalue -d, attaining Phi(A)=nd/2. The normalized source norm thus
tends to 1/sqrt(2), not 2/5.

A Boolean tensor eigenvector in the pair-zero-sum subspace likewise
attains beta(B)=2n sqrt(t). The general cross-only lower and the actual
K operator upper give beta(B)<=Phi(K)<=nd, so its normalized paired
norm tends to sqrt(2), not 4/3. The cross endpoint mass is 1/2, not
the strengthened formal profile's 9/25.

Verdict: PASS, no required correction. The construction refutes only
the generic upgrade from exact parity/rank obstruction to asymptotic
rank exclusion under the listed scalar-feasibility/cap conditions.
It fails the decisive small original-source cap and is neither an
actual counterexample to that problem nor a new closure milestone.
Its proper role is the narrowly scoped companion barrier described
in the source, preventing reuse of the invalid stability inference.

No mathematical computation, solver, spectral scan, construction search,
or numerical evaluation was run on any host. Tools were used only to
read the source, verify its hash, and write/read/hash this /tmp receipt.
No canonical repository file was edited and no publication was performed.
