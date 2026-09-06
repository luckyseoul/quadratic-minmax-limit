# Near-scalar internal flat-law transfer on the original source

2026-09-06. Supporting analytic lemma for an actual source argument.
No mathematical computation, construction, or search was run.

This note transfers a specified full weighted internal spectral law to
one large complete principal submatrix of the SAME original source A.
It does not replace the paired covariance, cross block, or active field.
The prospective strict-gain theorem to be applied afterward is separate;
no source-profile exclusion or original MO closure is asserted here.

## 1. Hypotheses and conclusion

Let n tend to infinity, and let

    K=[[A,B],[B^T,-A]]

be an actual complete symmetric zero-diagonal signing of order N=2n.
Let D=diag(D_L,D_R)>0 satisfy D+-K>=0, and set

    S=tr D,       dbar=S/(2n),
    delta=S tr(D^(-1))/(2n)^2-1,
    H_L=D_L^(-1/2) A D_L^(-1/2),
    alpha_n=Phi(A)/n^(3/2),
    Phi(A)=max_(x in {+-1}^n)|x^T A x|/2.

Assume delta tends to zero and the empirical eigenvalue law of the FULL
actual H_L converges weakly to

             chi_*=(9/25)delta_0+(8/25)(delta_(3/4)+delta_(-3/4)).

Then there is a common set J of original left/right index labels, with
q=|J| and q/n tending to one, such that the actual complete signing A_J
satisfies

    dbar/sqrt(q) -> 5/3,
    limsup ||A_J||op/sqrt(q) <= 5/3,
    empirical law of A_J/sqrt(q)
             -> (9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4)),
    Phi(A_J)/q^(3/2) <= alpha_n+o(1).                         (1.1)

In particular a separately assumed bound S<=C(2n)^(3/2) suffices, but
is not necessary: the specified nonzero internal second moment itself
forces the needed source scale. Only a bounded operator norm is claimed;
weak empirical convergence does not force its limit to equal 5/4.

The proof specializes the good-coordinate, principal-interlacing, and
congruence mechanism in
`original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.

## 2. One common original principal source

Write t_i=d_i/dbar, i=1,...,2n. The exact identity is

                 (1/(2n))sum_i (t_i-1)^2/t_i=delta.

If delta>0, choose epsilon=delta^(1/3), eventually at most 1/8, and put

    J={i in {1,...,n}: |t_i-1|<=epsilon and |t_(n+i)-1|<=epsilon},
    a=q/n,       b=1-a.

Each excluded index label has at least one bad coordinate. Since
(t-1)^2/t>=epsilon^2/(1+epsilon) outside the good interval,

          0<=b<=2delta(1+epsilon)/epsilon^2 -> 0.              (2.1)

If delta=0, take epsilon=0 and J={1,...,n}; all the same estimates hold
without divisions by epsilon. This convention also covers sequences
containing both zero and nonzero dispersion terms.

Let H_J be the actual principal compression H_L[J,J] and
Q=diag(sqrt(t_i):i in J). Then

    ||H_L||op<=1,       ||H_J||op<=1,
    A_J/dbar=Q H_J Q,
    ||A_J||op<=(1+epsilon)dbar,
    ||A_J/dbar-H_J||op<=3epsilon.                            (2.2)

These follow from principal feasibility and
||Q-I||op<=epsilon, ||Q||op<=sqrt(1+epsilon). Both good sets were
intersected so the very same A_J is available for the right diagonal
half as well; no independent left/right source matrices are substituted.

## 3. Full internal law and completeness determine the scale

Let F_L and F_J be the empirical eigenvalue distribution functions of
H_L and H_J, normalized respectively by n and q. Interlacing gives

             a F_J(x)<=F_L(x)<=a F_J(x)+b,
             sup_x |F_L(x)-F_J(x)|<=b.                       (3.1)

Thus the empirical law of H_J has the same weak limit chi_*. Since
both weighted spectra are in [-1,1], their second moments differ by
at most 2b. One way to see this finite bound is to treat the positive
and negative eigenvalue squares separately: each unnormalized moment
loses between zero and n-q under the compression, and then normalize.

Weyl's ordered-eigenvalue inequality and (2.2) show that A_J/dbar has
the same weak limit. On [-(1+epsilon),1+epsilon], x^2 has Lipschitz
constant at most 3, so the finite second-moment comparison is

    |(q-1)/dbar^2-(1/n)tr H_L^2|<=2b+9epsilon.              (3.2)

Here the left moment is EXACT: A_J is a complete zero-diagonal signing,
so tr A_J^2=q(q-1). Uniform support of H_L makes weak convergence
imply convergence of its second moment to

                  (16/25)(3/4)^2=9/25.

Because q tends to infinity, (3.2) therefore gives
dbar/sqrt(q)->5/3. Rescaling the empirical law of A_J/dbar by this
convergent scalar proves the spectral limit in (1.1), and (2.2) proves
the stated operator bound.

## 4. The original quadratic norm is retained

Fix any Boolean signing on J and extend it to all n labels by independent
unbiased signs. Its conditional expected original quadratic energy is
exactly its A_J energy: every other edge has zero expected contribution.
Taking the maximizing sign and then absolute values yields the exact
inequality

                           Phi(A_J)<=Phi(A).

Consequently

                    Phi(A_J)/q^(3/2)<=a^(-3/2)alpha_n.       (4.1)

The previous section already gives dbar/sqrt(n)->5/3. Feasibility
implies Phi(A)<=tr(D_L)/2<=S/2=n dbar, so alpha_n is bounded.
Since a tends to one, (4.1) is alpha_n+o(1), proving the last part of
(1.1). If alpha_n tends to 2/5, any independently proved strict lower
bound greater than 2/5 for complete sources with the limiting law in
(1.1) and a bounded normalized operator norm would therefore exclude
this particular ACTUAL near-scalar internal-law regime. Such a lower
bound is not proved or presumed in this transfer lemma.

During derivation no canonical repository file was changed and no
publication was performed. The full original weighted law remains the
input, and A_J is used only for an original-source norm lower bound.
