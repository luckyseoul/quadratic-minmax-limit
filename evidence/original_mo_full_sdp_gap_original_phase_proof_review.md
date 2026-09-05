# Proof-worker full review: full-SDP gap and actual original phases

2026-09-05. Complete analytic PASS; no correction requested.
No mathematical computation or source edit was performed by this reviewer.

Reviewed source: all 274 lines of
`/tmp/original_mo_full_sdp_gap_original_phase_bound.md`, SHA256
`1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`.

## Actual phase construction and the original-norm constant

K^2 plus or minus K|K| is twice the square of the corresponding
spectral phase. Its positivity gives |h_i|<=q coordinatewise, including
when K is singular. Adding the stated nonnegative diagonal padding and
using v_i=q+|h_i| gives two genuine unit-diagonal PSD correlations.
The same coordinate weights occur in both phases, which is essential
to the cancellation of K^2 in their difference.

The reviewer checked the factor-two normalizations directly. For each
phase the arcsine error is at most rho/2 times its off-diagonal squared
Frobenius norm. These norms are bounded respectively by 4 tr(K_+^4)/q^2
and 4 tr(K_-^4)/q^2. The linear difference is kappa J, with J summing
over ordered off-diagonal indices. Since each original expected
quadratic energy belongs to [-Phi(K),Phi(K)], dividing their difference
by two gives exactly Phi(K)>=kappa J/2-rho S_4/q^2.

No rectangular-to-quadratic factor is used to establish kappa/2.
The diagonal padding does not enter either trace pairing with K.

## Coordinate normalization and weighted mask

The diagonal of DK-K|K| is exactly -h. Its weighted Frobenius bound
therefore gives both inequalities (3.2), using weighted Cauchy--Schwarz
for the second. It controls the sum of |h_i|, not its maximum.

The interval f_i in [1/sqrt(2),1], the mask interval [0,1/2], and the
bound for A_0 follow from the displayed elementary scalar inequality.
The bound B_0<=S/2 correctly counts exactly q off-diagonal terms in
each row. No diagonal sign-square identity is used.

The arithmetic-harmonic dispersion gives the L1 deviation of d from
its mean. Because its signed sum is zero, its positive part has half
that L1 norm. This justifies the c_*/2 constant in the weighted a_i
bound. Expanding the two terms of m_ij<=a_i+a_j gives a coefficient
1-1/q on sum d_i a_i, so dropping its negative diagonal correction
indeed gives the stated safe bound (3.5).

Substituting sqrt(Sg/q)=eta N sqrt(gamma) and
sqrt(delta)<=2eta sqrt(gamma) yields exactly the coefficient
eta(1+N/q+c_*) in b_*.

## Exact residual identity and asymptotic statement

The identity (4.1) follows by inserting H=DK-R_1 only on off-diagonal
indices into the normalization loss. Its sign and the retained
canonical value S-g are correct. The weighted residual norm, combined
with m_ij^2<=m_ij/2, gives exactly sqrt(2gB_0), rather than an
unjustified unweighted or maximum-diagonal substitution.

The map b+sqrt(2gamma b) is increasing, so replacing B_0/S by b_* is
valid. The fourth moment bound S_4<=L S_3 gives the exact final term
rho(L/q)(1-gamma). These checks establish the complete finite (1.2),
including parameter ranges where its lower bound is negative.

Under the stated fixed norm cap, the cited interpolation bound gives
L/q=O(N^(-1/4)), and the finite Grothendieck inequality only bounds
S=O(N^(3/2)). The mask and mixed terms are O(sqrt(gamma)) and
O(gamma^(3/4)); both, as well as gamma, have the stated safe
O(sqrt(gamma)) bound on [0,1]. This proves (1.3) with the recorded
N^(5/4) finite term.

## Paired application and limitation

The paired application retains c/Phi(K) unless original p=q_A=0 and
c=Phi(K) are separately assumed. With that explicit saturation and
the uniform compatibility theorem, the conclusion c_D/n>=kappa-o(1)
is correct. The construction does not itself evaluate a field width.

For gamma>=1/4, eta>=1 and N/q>1 force b_*=1/2. The main bracket is
1/2-gamma-sqrt(gamma), at most -1/4. Thus the displayed estimate is
indeed vacuous there. This is only a limitation of this estimate,
not a claim that no other positive-gap argument can work.

## Review provenance

The exact worker authored this theorem. Before receiving its source,
the reviewer requested the small-gap consequence for a separate field
evaluation and discussed its zero-gap motivation and the need to avoid
a maximum-coordinate phase normalizer. The reviewer did not author or
edit the coordinate-normalized construction or masked-residual proof.
The reviewer did author the separately cited residual/compatibility
prerequisite and earlier original-phase spectral work. This receipt
records the subsequent full-source analytic check with that dependency
involvement and prior statement discussion disclosed.
