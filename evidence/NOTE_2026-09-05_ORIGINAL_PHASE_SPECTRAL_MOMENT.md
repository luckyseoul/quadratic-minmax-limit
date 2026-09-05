# Original quadratic-norm spectral moments from actual Gaussian phases

2026-09-05. Analytic theorem for the ORIGINAL quadratic norm of a
complete signing. The positive and negative Gaussian phases are actual
positive-semidefinite covariances. Their generally unequal diagonal
normalizations are retained. No spectral symmetry, source optimality,
or comparison with a doubled-order optimizer is assumed.

The construction was checked against the existing Gaussian lower bound
in `CORE.md`, Section 4, and the source-preserving whole-edge note.
This is a different same-order constraint, not an order-transport upper.
No computational experiment is used.

## 1. Definitions and exact phase bounds

Let `n>=2`, let `A` be a complete symmetric zero-diagonal signing, and
write
\[
 q=n-1,\qquad
 \Phi(A)=\max_{x\in\{-1,1\}^n}\left|\tfrac12x^TAx\right|,
 \qquad \kappa=2/\pi,\quad\delta=1-\kappa.
\]
The positive and negative spectral parts are
\[
 A_+=(|A|+A)/2,\qquad A_-=(|A|-A)/2,
\]
and put
\[
 a=\lambda_{\max}(A)>0,\qquad b=-\lambda_{\min}(A)>0,
 \qquad S_j^\pm=\operatorname{tr}(A_\pm^j),
\]
\[
 T=A|A|,\qquad h_i=T_{ii}.
\]
Both extremes are positive since `A` is nonzero and has trace zero.

Choose any positive numbers satisfying
\[
           v_+\ge q+\max_i h_i,\qquad
           v_-\ge q-\min_i h_i.                          \tag{1}
\]
Then the following exact inequalities hold:
\[
 \boxed{\quad
 v_+\Phi(A)\ge\kappa S_3^+-{2\delta S_4^+\over v_+},
 \qquad
 v_-\Phi(A)\ge\kappa S_3^--{2\delta S_4^-\over v_-}.
 \quad}                                                  \tag{2}
\]
In particular,
\[
 \boxed{\quad
 (v_++v_-)\Phi(A)
 \ge\kappa\operatorname{tr}|A|^3
      -2\delta\left({S_4^+\over v_+}+{S_4^-\over v_-}\right)
 \ge\kappa\operatorname{tr}|A|^3-\delta n(a^2+b^2).
 \quad}                                                  \tag{3}
\]

### Proof by actual Gaussian phases

Define
\[
 C_+={A^2+T+\operatorname{diag}(v_+-q-h_i)\over v_+},
 \qquad
 C_-={A^2-T+\operatorname{diag}(v_--q+h_i)\over v_-}.       \tag{4}
\]
The matrices `A^2+T=2A_+^2` and `A^2-T=2A_-^2` are PSD.
The added diagonals are nonnegative by (1). Because `(A^2)_{ii}=q`,
both `C_+` and `C_-` have diagonal one. They are therefore genuine
Gaussian correlation matrices, allowing singularity.

Let `g_+` and `g_-` have these covariances and put
`x_+=sign(g_+)`, `x_-=sign(g_-)`. They may be sampled separately;
no particular joint coupling is required. The Gaussian sign identity
from `CORE.md` gives
\[
 \mathbb E Q_A(x_\pm)
       =\kappa\sum_{i<j}A_{ij}\arcsin(C_\pm)_{ij}.        \tag{5}
\]
For every `z in [-1,1]`,
\[
                 |\arcsin z-z|\le(\pi/2-1)z^2.           \tag{6}
\]
Indeed the positive Taylor coefficients after the linear term have
sum `pi/2-1`, and each remaining power of `|z|` is at most `z^2`.
The endpoints follow by continuity.

Since `A_{ii}=0`, its trace pairing with each padding diagonal is zero.
Thus the linear terms in (5) are respectively
\[
      {\kappa\over2}\operatorname{tr}(AC_+)
          ={\kappa S_3^+\over v_+},\qquad
      {\kappa\over2}\operatorname{tr}(AC_-)
          =-{\kappa S_3^-\over v_-}.                    \tag{7}
\]
The absolute remainder in either expectation is at most
\[
 \delta\sum_{i<j}(C_\pm)_{ij}^2
 \le {\delta\over2v_\pm^2}\|2A_\pm^2\|_F^2
 ={2\delta S_4^\pm\over v_\pm^2}.                       \tag{8}
\]
In (8), the off-diagonal part of the covariance is the off-diagonal
part of `2A_pm^2/v_pm`; the diagonal padding is not included in its
Frobenius estimate.

Every realization satisfies `|Q_A(x_pm)|<=Phi(A)`. Apply this to the
positive phase and to the negative of the negative-phase expectation
in (7)--(8), then multiply by its OWN variance `v_pm`. This proves
(2), and adding proves the first inequality in (3). In particular the
two phases are never silently assigned the same normalization.

Finally, `S_4^+<=a^2 S_2^+` and
\[
              2S_2^+=\sum_i(q+h_i)\le n v_+.
\]
Therefore `2S_4^+/v_+<=na^2`. The corresponding negative bound is
`2S_4^-/v_-<=nb^2`, proving the final inequality in (3).

## 2. The exact diagonal-imbalance denominator

The smallest permissible choices in (1) are strictly positive, because
`q+h_i=2(A_+^2)_{ii}` and `q-h_i=2(A_-^2)_{ii}`, and both spectral
parts are nonzero. With these choices,
\[
 v_++v_-=2q+\operatorname{osc}_i h_i,
 \qquad \operatorname{osc}_i h_i=\max_i h_i-\min_i h_i.
\]
Consequently
\[
 \boxed{\quad
 [2(n-1)+\operatorname{osc}_i(A|A|)_{ii}]\Phi(A)
 \ge\kappa\operatorname{tr}|A|^3-\delta n(a^2+b^2).
 \quad}                                                  \tag{9}
\]
In particular, a constant diagonal of `A|A|` yields denominator
`2(n-1)`, even if that constant is nonzero and the spectrum is not
symmetric. This is not an assumption available for every signing.

Under `a,b<=K sqrt(n)` with fixed K, the last error in (9) is
`O_K(n^2)`; after division by its denominator it is `O_K(n)`.
It is subleading on the original `n^(3/2)` norm scale. With no such
operator control, (9) remains exact, but that asymptotic interpretation
of the error is not asserted.

## 3. A sharp row-moment bound from zero matrix diagonal

For every i,
\[
 \boxed{\quad
 -{q(b^2-q)\over b^2+q}\le h_i
                 \le{q(a^2-q)\over a^2+q}.
 \quad}                                                  \tag{10}
\]
The quantities on the right or left may themselves have either sign;
the assertion does not assume both `a,b >= sqrt(q)`.

To prove the upper bound, use a spectral resolution at coordinate i.
This supplies a finitely supported random variable `X` in `[-b,a]`
with
\[
          \mathbb EX=0,\quad\mathbb EX^2=q,
                       \quad\mathbb E[X|X|]=h_i.
\]
Let `p_+,p_-` be the probabilities of positive and negative X. Write
\[
 t=\mathbb EX_+=\mathbb EX_->0,
 \qquad P_2=\mathbb EX_+^2,\quad N_2=\mathbb EX_-^2.
\]
Then `P_2+N_2=q`, `h_i=P_2-N_2`, and `P_2<=a t`.
Also `p_+>=t/a`, while Cauchy--Schwarz gives `p_->=t^2/N_2`.
Because `p_++p_-<=1`,
\[
 1\ge {t\over a}+{t^2\over N_2}
   \ge {P_2\over a^2}+{P_2^2\over a^2N_2}
    ={qP_2\over a^2N_2}.
\]
It follows that `N_2>=q^2/(a^2+q)`, proving the upper bound in (10).
Apply the same argument to `-X` for the lower bound. Equality in the
one-sided scalar bound is attained by the two-point distribution on
`a` and `-q/a` when that support is permitted; no matrix equality
case is presumed.

Thus universally valid choices in (1) are
\[
           v_+={2qa^2\over a^2+q},\qquad
           v_-={2qb^2\over b^2+q}.                       \tag{11}
\]
These are upper bounds on the actual largest variances, not assumptions
that those variances are attained at every coordinate.

Combining (3) with (11) yields the fully spectral constraint
\[
 \boxed{\quad
 2q\left({a^2\over a^2+q}+{b^2\over b^2+q}\right)\Phi(A)
 \ge\kappa\operatorname{tr}|A|^3-\delta n(a^2+b^2).
 \quad}                                                  \tag{12}
\]
The local identity `EX=0, EX^2=q` also gives `ab>=q`, by taking
expectations in `(a-X)(X+b)>=0`. Consequently the factor in
parentheses in (12) is between one and two. Its lower endpoint is
attained when `ab=q`; no such endpoint condition is imposed generally.

If `a,b<=K sqrt(q)`, with necessarily `K>=1`, (12) in particular gives
\[
 \operatorname{tr}|A|^3
 \le {4qK^2\over\kappa(K^2+1)}\Phi(A)
                         +{2\delta\over\kappa}nK^2q.    \tag{13}
\]
At fixed K the leading coefficient is strictly better than the
`4q/kappa` obtained by first bounding the bilinear norm by `4 Phi(A)`.
This is a bounded-source improvement, not a uniform replacement of that
bilinear inequality for all symmetric signings.

## 4. A complementary nuclear-norm constraint

There is also an exact refinement of the uniform lower bound in CORE:
\[
 \boxed{\quad
 \Phi(A)\ge\kappa\sum_{i<j}{1\over\sqrt{d_i d_j}}
          \ge{\kappa n^2(n-1)\over2\operatorname{tr}|A|},
 \qquad d_i=|A|_{ii}>0.
 \quad}                                                  \tag{14}
\]
For the first inequality, the PSD covariances `|A|+A` and `|A|-A`
have the SAME diagonal `d_i`, since `A_{ii}=0`. Normalize them by
that common diagonal. Their correlations `r_ij^+,r_ij^-` satisfy
\[
                   A_{ij}(r_{ij}^+-r_{ij}^-)
                           ={2\over\sqrt{d_i d_j}}.
\]
Oddness of arcsine and the elementary slope bound
`arcsin(u)-arcsin(v)>=u-v` for `u>=v` in `[-1,1]` show that the
difference of their expected ORIGINAL quadratic energies is at least
`2 kappa sum_{i<j}(d_i d_j)^(-1/2)`. Both expectations lie in
`[-Phi(A),Phi(A)]`, proving the claim. The normalization is legitimate:
if `|A|_{ii}=0`, positivity would give a zero row of `|A|`, hence
`(A^2)_{ii}=0`, contrary to q>0.

For the second inequality, use
`1/sqrt(d_i d_j)>=2/(d_i+d_j)` and Cauchy--Schwarz on all unordered
pairs. Their sum of denominators is `(n-1) tr|A|`, which yields
the printed constant. Finally `d_i<=sqrt(q)` by Cauchy--Schwarz in
the local spectral measure. Thus (14) recovers
`Phi(A)>=n sqrt(n-1)/pi` and retains additional diagonal/nuclear
information when `|A|` is not spectrally flat.

## 5. Scope for the live upper comparison

Equations (9), (12), and (14) concern every actual source signing.
They may be imposed on exact or near-original minimizers without
transferring an optimality property to a new random posterior.
For bounded-operator near-minimizers, (12) supplies a leading cubic
moment constraint with an explicit subleading error. Equation (9)
is sharper if the actual diagonal odd-spectral imbalance is controlled.

None of these inequalities evaluates the maximum over attainable
joint shells of the paired Gaussian proposal. In particular they do
not prove `F_A^*<=2 sqrt(2) Phi(A)+o(n^(3/2))`, nor the original
all-orders convergence assertion. The source and cross spectral
constraints are necessary information for that ongoing evaluation,
not sufficient closure statements.
