# Independent full-source review: near-scalar actual cross gain

2026-09-05. Verdict: PASS. No mathematical correction requested.

## 1. Frozen source, prerequisites, and independent role

I directly read the complete final 364-line source
`/tmp/original_mo_near_scalar_cross_spectral_gain.md`, SHA256

    ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9

This is the version whose interlacing argument explicitly uses the
monotone function `(max(x,0))^(2k)`. The earlier hash is superseded.
I had no role in deriving this transfer or its finite factors, balancing
rule, moment estimates, or active-state consequence. The source discloses
the root/proofer/exact contributions; this review is independent of them.

The substantive imported theorem is the complete 411-line actual
complete-cross source, already read and independently reviewed in full
by me in this working interval:

    original_mo_complete_cross_flat_spectral_gain.md
    b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20

My independent receipt for that prerequisite is 201 lines, SHA256
`01dacaf0e4d01edaef3f3b85651748ca5475b2a0b0c83b84b73fb4307e721a3f`.
For the present review I also directly reread both complete contextual
sources identified in Section 9 and checked their hashes:

    original_mo_near_scalar_diagonal_spectral_normalization.md
    280 lines
    c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2

    original_mo_full_sdp_gap_weighted_compatibility.md
    303 lines
    3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0

The transfer rederives the needed good-coordinate and uniform-energy
facts. In particular it does not import trace optimality from the
303-line theorem into the present every-feasible-D compatibility bound.
No other unpublished refinement is needed for this proof.

## 2. Actual normalization and balanced complete block

The row-square inequality for the contraction T is
`sum_(j!=i)1/d_j<=d_i`. Summing it proves
`(N-1)tr(D^(-1))<=S`; scalar Cauchy--Schwarz then gives
`dbar^2>=N-1`. The separate trace cap supplies
`dbar^2<=C^2 N`. With N=2n these give precisely
`1/(2C^2)<=m_0<=n/(2n-1)<=1`.

The exact dispersion identity follows from averaging
`(t-1)^2/t=t-2+1/t`. Literal squared cross signs factor the full
weighted second singular moment as `m=m_0 ell h`. The inverse means
sum to `2(1+delta)`, so their product is at most `(1+delta)^2`.
For the lower direction, half-wise Cauchy--Schwarz gives
`ell h>=1/(t_L t_R)>=1`, since `t_L+t_R=2`. Therefore
`m_0<=m<=m_0(1+delta)^2` is exact. Contraction of W gives
`0<v_2<=m<=1`; positivity follows from the nonzero complete cross block.

The bad-coordinate bound uses the global count over both halves.
If their bad counts are b_L and b_R, then
`q=n-max(b_L,b_R)` and `theta=max(b_L,b_R)/n<=2b_0`.
This verifies the factor two, the balancing rule, `a=q/n>=1/2`,
and the integer bound `q>=ceil(n/2)`. Removing surplus good coordinates
from one half does not spoil its retained diagonal interval.

For the exact auxiliary scale and mean,

    d'=(1+eta)dbar,
    m'=q/d'^2=a m_0/(1+eta)^2>=1/(9C^2).

The last constant uses `a>=1/2`, `m_0>=1/(2C^2)`, and
`(1+eta)^2<=9/4`. The original complete B_J genuinely has every
entry in {+1,-1}. Compressing W is contractive, and multiplying its
two retained diagonal square roots bounds only `||B_J||op<=d'`.
No operator bound on the entire unweighted B is inferred.

The unbiased independent extension of signs annihilates every omitted
cross term, including terms with both coordinates removed. It proves
`beta(B)>=beta(B_J)` without making B_J an optimizer or substituting
it into an actual field or cell.

## 3. Full actual second and fourth singular moments

The symmetric dilation of W_J is indeed a principal submatrix of the
dilation of W, with exactly `2(n-q)` deleted coordinates. For a
symmetric contraction, interlacing and the nondecreasing function
`(max(x,0))^(2k)` show that deleting r coordinates loses between zero
and r in its positive-power eigenvalue sum. Both dilation spectra are
paired, so those positive-power sums are precisely the cross singular
powers, counted once. Division by n gives (4.1), with loss `2theta`.
This treats k=1 and k=2, meaning the first two moments of the SQUARED
singular-value law: the second and fourth singular powers themselves.

On retained coordinates, `sqrt((1-eta)/(1+eta))>=1-eta` for
the stipulated eta range. Thus both P factors lie between `1-eta`
and one. Expanding `P_L W_J P_R-W_J` with two perturbation terms
gives operator error at most `2eta`, while both matrices remain
contractive. Singular-value Weyl comparison and the `2k` Lipschitz
constant of `x^(2k)` on [0,1] yield error `4k a eta` after
normalization by n. All factors in (4.3) are therefore correct.

Literal sign squares give `v_1'=m'` exactly. The fourth-moment lower
bound can be replaced by its positive part because `a v_2'>=0`.
The full moments m and v_2 always remain those of the original W,
including every exceptional coordinate and zero singular value.

## 4. Imported gain, finite factors, and explicit constant

For B_J, direct substitution gives `1-epsilon'=v_2'/m'`.
The change of norm scale is exactly
`q d'/(n dbar)=a(1+eta)`. Applying the 411-line theorem and then
the full fourth-moment comparison produces the coefficient
`kappa(1+eta)/m'` in (5.3), with both normalization factors retained.

The key one-sided comparison is valid: `m'<=m_0<=m` and the bracket
is nonnegative after taking its positive part. Hence its coefficient
can safely be reduced to `kappa/m`, and only then can the positive
part be bounded below by its untruncated argument. This yields (6.1)
without multiplying a negative error by an enlarged uncontrolled factor.

The second contribution is exactly `g_kappa a^2 m_0/(1+eta)`.
The directional upper bound on m, followed by `m_0<=1`, gives its
loss at most `g_kappa(2delta+delta^2+2theta+eta)`. The bracket is
nonnegative, and `1-a^2<=2theta`, so the coefficient directions in
this step are correct.

For `eta=delta^(1/3)<=1/8`, the stated good-coordinate condition
holds and `theta<=9eta/4`. Consequently

    2theta+8a eta<=25eta/2,
    2delta+delta^2+2theta+eta<=6eta.

For the latter, division by eta bounds the first two terms by
`2eta^2+eta^5`, which is less than the available extra one-half
when eta<=1/8; the other terms total at most 11/2.
Combining the first bound with `1/m<=2C^2` gives precisely
`25kappa C^2 eta`, not a missing factor of two. The second gives
`6g_kappa eta`. This checks the printed explicit constant in (1.1).

The prerequisite uniform error can be chosen bounded and nonnegative:
its absolute-moment proof has a uniform second-moment bound at each
fixed lower-m parameter. Its tail supremum over integers
`j>=ceil(n/2)` is therefore finite and tends to zero. The factor
3/2 dominates `a(1+eta)`, and q belongs to that tail, so the defined
R_C(n) controls every auxiliary choice and every allowed dispersion.

At delta=0, the exact dispersion identity forces every t_i=1.
Using q=n and eta=theta=0 directly invokes the original theorem with
the full B. The envelope includes j=n and dominates its error too.
Thus the exact scalar endpoint is covered without taking a singular
trimming limit or claiming a finite-n convergence rate.

## 5. Every-feasible-D active ratio and scope

For Q=diag(sqrt(t_i)), the dispersion identity implies
`||(Q-I)z||^2<=N delta` for every Boolean z, while
`||Qz||=||z||=sqrt(N)`. Since `K/dbar=Q T Q` and T is
contractive, the quadratic difference is at most `2N sqrt(delta)`.
The Phi convention divides this by two. Comparing the two full block
energies with one half flipped therefore bounds the cross pairing by
`N sqrt(delta)`, and division by n gives exactly `2sqrt(delta)`.
Neither optimality of D nor a canonical-primal gap is used here.

The inequality `beta(B)<=Phi(K)` follows from the same half-flip
argument. At the separate actual state with p=q_A=0 and c=Phi(K),
the reverse comparison `beta(B)>=c` forces equality. Therefore (7.2)
has the correct active-state error, and its restriction to actual
within-final-cell representatives is appropriate. No ratio c/beta(B)
is silently omitted on a general state.

Near-flatness of the FULL weighted law means `(m-v_2)/m->0`.
The fixed cap keeps m bounded below, so (1.1) and, on the active
face, (7.2) force a uniform positive leading gain above kappa.
The earlier formal flat endpoint is consequently excluded within
this near-scalar fixed-cap actual branch, even without a global
unweighted operator bound. Its relaxation-only certificate statement
is not retracted or converted into an actual signing counterexample.

The auxiliary block never replaces the original W, covariance, source,
or cell. Small dispersion remains a hypothesis; the proof does not
establish it for every optimizer, reduce all active cells to pure-cross,
or evaluate the remaining full-law ellipsoid upper. Original all-cell
inequality and convergence remain open.

No mathematical computation, test, checker, spectral scan, signing
search, or canonical edit was performed by me. This independent
receipt is the only new artifact I wrote for the review. Full-source
analytic and scope verdict at the stated final hash: PASS.
