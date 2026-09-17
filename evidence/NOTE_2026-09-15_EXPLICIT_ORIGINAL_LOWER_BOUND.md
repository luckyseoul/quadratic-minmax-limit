# An explicit unconditional lower bound for the original sequence

2026-09-15. Analytic quantitative refinement of the September 6 strict
lower bound. For the original complete symmetric zero-diagonal signings,

\[
 \liminf_{n\to\infty}\frac{m_n}{n^{3/2}}
 \;>\;\frac1\pi+10^{-6}.                                  \tag{1}
\]

The error estimates below are uniform in the signing. They do not give
an explicit order after which the bound holds. Convergence and the value
of a possible limit remain open.

The new ingredient is a quantitative spectral-sign rounding estimate,
followed by an update whose quadratic cost is bounded by the original
Boolean norm. No operator bound on the original signing is assumed.

## 1. Notation and reused results

Let A be any complete symmetric zero-diagonal signing of order n>=2.
Write

    M=A/sqrt(n),       q=(n-1)/n,       kappa=2/pi,
    Q_M(x)=x^T M x/2,
    alpha=Phi(A)/n^(3/2)=max_(x in {+1,-1}^n)|Q_M(x)|/n,
    ell=tr|M|/n,       d=q+1-2ell.

Completeness gives `tr M^2=nq`, and Cauchy--Schwarz gives
`0<ell<=sqrt(q)`. The exact nuclear lower from
[the all-law source note, Section 4](NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md)
is

    alpha >= kappa q/(2ell).                               (2)

We also reuse the uniform scalar assertion in
[the local Gaussianization lemma, Section 6](NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md).
For any correlation matrix R with `||R||op<=4`, G Gaussian with covariance
R, and `F=sum_j a_j sign(G_j)` with `max|a_j|<=1/sqrt(n)` and
`sum_j a_j^2<=1`, there is one error `e_4(n)>=0`, tending to zero,
such that

    |E|F|-sqrt(kappa)*sqrt(E F^2)| <= e_4(n).                (3)

This lemma does not require an operator cap on the matrix whose row
provides the coefficients. Its proof uses bounded Schur powers, Hermite
contraction estimates, then a uniformly controlled Hermite tail. It is a
scalar statement; no joint limit for all n fields is needed here.

The existing rational enclosure `7/11<kappa<16/25` is reused throughout.

## 2. A bounded Gaussian frame close to the spectral sign

Let S be the symmetric spectral sign of M, assigning +1 on its kernel.
Then

    S^2=I,       ||S||op=1,       MS=|M|,
    ||M-S||_F^2=nd,       0<=d<=2.

Since M has zero diagonal, `sum_i S_ii^2<=nd`. Define

    I={i:|S_ii|<=1/2},       B=I^c,       P=I_n+S.

Here I_n denotes the identity matrix, distinct from the index set I.
On I, normalize the positive semidefinite principal matrix P[I,I] by
its own diagonal `1+S_ii`. Put an independent identity block on B, with
zero cross blocks, and call the resulting correlation matrix R. Thus

    diag R=1,       0<=R<=4I_n.                             (4)

We need the quantitative bound

    ||R-(I_n+S)||_F <= 4sqrt(nd).                           (5)

To prove it, put `U=diag((1+S_ii)^(-1/2):i in I)`. For |t|<=1/2,

    |(1+t)^(-1/2)-1| <= 2(sqrt(2)-1)|t|.

This follows by rationalizing the difference; the denominator
`sqrt(1+t)(1+sqrt(1+t))` is smallest at t=-1/2. Since `||U||op<=sqrt(2)`
and `||P[I,I]||op<=2`,

    ||U P[I,I] U-P[I,I]||_F
      <=2(sqrt(2)+1)||U-I||_F
      <=4sqrt(sum_(i in I) S_ii^2).

On the remaining entries, R-P consists of the two removed cross blocks
of S and the block -S[B,B]. Every row of S has squared length one, so
the squared Frobenius cost is at most `2|B|`, which is at most
`8 sum_(i in B)S_ii^2`. These entries are disjoint from the good principal
block. Adding the squared costs proves (5).

Every original coordinate remains present. S and R only specify the
rounding distribution; the energy matrix is always M.

## 3. The original phase energy

Let `G~N(0,R)`, `X=sign(G)`, and `C_X=E XX^T`. The Gaussian sign
identity and its nonnegative odd Hermite coefficients give

    C_X=kappa R+C_tail,       0<=C_X<=4I_n,
    C_tail=sum_(odd k>=3)c_k^2 R^(circ k),
    sum_(odd k>=3)c_k^2=1-kappa.

For every k>=1, Schur multiplication gives
`0<=R^(circ k)<=4I_n`. For odd k>=3, the complete off-diagonal entries
of M give

    |tr(M R^(circ k))|
      <=n^(-1/2)sum_(i!=j)|R_ij|^k
      <=n^(-1/2)tr R^2 <=4sqrt(n).

This estimate is uniform in k, so it can be summed over the whole
Hermite tail. By (5), `tr M=0`, `tr(MS)=nell`, and `||M||_F=sqrt(nq)`,
the actual phase energy `e=E Q_M(X)/n` satisfies

    e >= kappa ell/2-2kappa sqrt(d)
                         -2(1-kappa)/sqrt(n).             (6)

This is the baseline of the actual X being updated below. It is not
replaced by the distinct maximized nuclear bound (2).

## 4. Quantitative local-field variances without a source operator cap

Set `a_n=1+n^(-1/2)+n^(-1)`. Off the diagonal,

    |(C_tail)_ij| <= (1-kappa)|R_ij|^3.

For 0<=r<=1 and h=n^(-1/2),

    r^3 <= (1+h+h^2)|r-h|+h^3.

If r<=h this is immediate; otherwise factor `r^3-h^3` and bound
`r^2+rh+h^2` by `1+h+h^2`. Also
`||R_off-M||_F<=5sqrt(nd)` by (5), because M has zero diagonal.
The reverse triangle inequality gives
`|| |R_off|-|M|_entry ||_F<=5sqrt(nd)`, where `|M|_entry` here
means entrywise absolute value off the diagonal, not the spectral |M|.
Since those entries are all h, the preceding scalar inequality yields

    ||C_X-(I_n+kappa S)||_F
      <=[4kappa+5(1-kappa)a_n]sqrt(nd)
                                      +(1-kappa)/sqrt(n). (7)

The diagonal part is accounted for by
`C_X=kappa R+(1-kappa)I_n+H`, where H has zero diagonal.

Write `v_i=E(SX)_i^2=(S C_X S)_ii`. Conjugation by S preserves the
Frobenius norm and fixes `I_n+kappa S`. Taking diagonal parts in (7),
and using `sum_i|S_ii|<=n sqrt(d)`, gives

    (1/n)sum_i|v_i-1|
      <=[5+5(1-kappa)(n^(-1/2)+n^(-1))]sqrt(d)
                                             +(1-kappa)/n.

For v>=0, `sqrt(v)>=1-|v-1|`. For the actual local fields
`F_i=(MX)_i`, the L2 triangle inequality and `C_X<=4I_n` give

    (1/n)sum_i |sqrt(E F_i^2)-sqrt(v_i)|
      <=sqrt(E||(M-S)X||^2/n) <=2sqrt(d).

Consequently, with

    b_n=[7+5(1-kappa)(n^(-1/2)+n^(-1))]sqrt(d)
                                               +(1-kappa)/n,

the uniform Gaussianization (3), applied to the actual rows of M, proves

    f:=(1/n)sum_i E|F_i| >= sqrt(kappa)(1-b_n)-e_4(n).

Since d<=2, this is uniformly

    f >= sqrt(kappa)(1-7sqrt(d))-o(1).                     (8)

The estimates involving SX compare second moments only. Gaussianization
is applied to MX, whose coefficients satisfy (3), never to the possibly
localized rows of S.

## 5. A Boolean update bounded by the original norm

For each realization of X, put `Y=sign(MX)`, choosing +1 at zero, and
`z=(1-p)X+pY` for 0<=p<=1. Multilinearity and zero diagonal show that
`|Q_M(z)|<=n alpha`: independently round its coordinates to signs with
means z_i and take the expectation. Expanding the same quadratic gives

    Q_M(z)=(1-p)^2 Q_M(X)+p(1-p) X^T M Y+p^2 Q_M(Y).

Symmetry implies `X^T M Y=||MX||_1`, and `Q_M(Y)>=-n alpha`.
Taking expectations therefore proves the exact inequality

    (1+p^2)alpha >= (1-p)^2 e+p(1-p)f.                    (9)

It uses the original Boolean norm to pay for the new state's quadratic
energy. No source operator bound or regularity of X -> Y is needed.

Take p=1/10 in (9), then insert (6) and (8). There is a nonnegative
error `r(n)->0`, independent of A, such that

    101alpha >= (81/2)kappa ell+9sqrt(kappa)
                    -(162kappa+63sqrt(kappa))sqrt(d)-r(n). (10)

For example r(n) is obtained directly by multiplying (6)'s displayed
finite error by 81 and (8)'s finite error by 9, using d<=2.
Equation (10) is a uniform all-source bound retaining the actual nuclear
moment, not only a statement about a hypothetical saturation sequence.

## 6. An explicit gap

Put epsilon=10^(-6). Suppose that some sequence of actual signings has
orders tending to infinity and `alpha -> a <= kappa/2+epsilon`.
By (2),

    liminf ell >= kappa/(kappa+2epsilon),
    limsup d <= 4epsilon/(kappa+2epsilon) < 7epsilon.

Here `4/kappa<44/7<7`. In particular `d<=7epsilon` eventually. Replacing
ell by `(q+1-d)/2=1-(d+1/n)/2` in (10) gives

    101(alpha-kappa/2)
      >=9sqrt(kappa)-10kappa-(81/4)kappa d
                    -(162kappa+63sqrt(kappa))sqrt(d)-o(1).

The already available rational enclosure for kappa proves

    sqrt(kappa)>797/1000,       sqrt(kappa)<4/5,
    9sqrt(kappa)-10kappa >773/1000,
    (81/4)kappa<13,       162kappa+63sqrt(kappa)<155,
    sqrt(7epsilon)<27/10000.

Thus the hypothetical limit must satisfy

    101(a-kappa/2)
      >=773/1000-13(7epsilon)-155(27/10000).

But the right side exceeds 101epsilon by exactly

    773/1000-13(7/10^6)-155(27/10000)-101/10^6
      =88577/250000 >0.                                   (11)

This is a contradiction. Apply it to a convergent subsequence attaining
the liminf of the bounded original minima to prove (1). No limiting
spectral law or bound on the original operator norm was imposed.

## 7. Provenance and verification scope

Reused inputs are (2), the complete proof of the scalar lemma underlying
(3), and the existing rational enclosure for kappa. The new argument is
the fixed-threshold spectral-sign frame estimate (5), its quantitative
covariance/field transfer (7)--(8), and the original-norm update (9).
Together these give (10) and the explicit contradiction (11).

The duplication check covered main, all nine other local branch heads,
the three linked worktrees, and the unpublished analytic frontier in the
September 6 campaign archive. The older strict lower is qualitative;
the older adaptive-update penalty uses an actual operator cutoff. No
previous version of (5)--(11) was found. No construction family or
signing census is used.

The author is the root assistant. Whole-proof self-review and the mesh
check are recorded separately under `explicit_original_lower_20260915/`.
The scalar computation checks arithmetic and update coefficients; it
does not formalize the matrix estimates or the imported Gaussian lemma.
No independent human review or formal verification is claimed.

## Addendum, 2026-09-17: independent verification and sharpened constant

The 2026-09-17 operating session re-derived every display above (the two
imported inputs and (4)--(11)), reproduced the margin (11) exactly from the
stated constants, checked (4)--(9) numerically on all signings of orders
2..7, on Paley conference families through order 158, on perturbed and
random signings, with Monte Carlo for (8) and exact update-bound checks for
(9), and replayed the checks on two further mesh nodes. No error was found.
Records: `explicit_original_lower_20260915/SELF_REVIEW.md` (line-by-line
review), `README.md` (results and hashes), `result_*.json` and
`replay_*.json` (receipts).

Sharpened constant, same proof: with the tighter rational enclosure of
`kappa` from the frozen enclosure `31415926/10^7 < pi < 31415927/10^7`
(recorded in `NOTE_2026-09-05_SOURCE_CROSS_NUCLEAR_TRACE_BOUNDARY.md`) and
the exact `d`-bound `d <= 4 eps/(kappa + 2 eps) <= 6.3 eps` in place of the
`7 eps` used above, the same contradiction certifies

    liminf_{n->infty} m_n / n^{3/2} > 1/pi + 4 * 10^{-6},

with exact rational margins for every `eps = 10^{-6}, 2*10^{-6}, 3*10^{-6},
3.5*10^{-6}, 4*10^{-6}` and with `4.5*10^{-6}, 5*10^{-6}` failing as
controls. The same-method ceiling is `eps = 305/68378763 = 4.460449...e-6`
(certified). Claim (1) above is therefore conservative; no analytic estimate
is changed and convergence and the value of a possible limit remain open.
Certificate: `explicit_original_lower_20260915/sharpen_constant.py`,
`result_sharpen.json`.
