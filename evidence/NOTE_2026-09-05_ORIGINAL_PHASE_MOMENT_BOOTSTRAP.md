# Norm-only bootstrap of the original phase spectral moment

2026-09-05. This supplements, without changing, the reviewed 262-line
`original_mo_original_phase_spectral_moment.md`, whose SHA-256 is
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
It removes a fixed conference-scale operator cap from the LEADING
cubic-moment interpretation, including the actual diagonal-imbalance
denominator. All assertions remain same-order source constraints.

## 1. Exact phase coefficients

Retain the source note's notation: `q=n-1`, `kappa=2/pi`,
`delta=1-kappa`, extreme eigenvalues `a,-b`, `L=max(a,b)`,
`S_j^pm=tr(A_pm^j)`, `S_j=tr|A|^j`, and `T=A|A|`.

The source note proves, for any admissible positive phase variances,
\[
 v_\pm\Phi(A)\ge\kappa S_3^\pm
                              -{2\delta S_4^\pm\over v_\pm}.
                                                               \tag{1}
\]
Since `S_4^+<=a S_3^+` and `S_4^-<=b S_3^-`, this gives the EXACT
coefficients
\[
 \boxed{\quad
 v_+\Phi(A)\ge\left(\kappa-{2\delta a\over v_+}\right)S_3^+,
 \qquad
 v_-\Phi(A)\ge\left(\kappa-{2\delta b\over v_-}\right)S_3^-.
 \quad}                                                  \tag{2}
\]
The inequalities hold even when a printed coefficient is nonpositive.
Dividing by such a coefficient is only permitted after proving it
strictly positive.

For the source note's universally admissible spectral choices
\[
             v_+={2qa^2\over a^2+q},\qquad
             v_-={2qb^2\over b^2+q},
\]
the coefficients in (2) become respectively
\[
       \gamma_a=\kappa-\delta\left({a\over q}+{1\over a}\right),
 \qquad\gamma_b=\kappa-\delta\left({b\over q}+{1\over b}\right).
                                                               \tag{3}
\]
Since `ab>=q`, both `a,b` lie in `[q/L,L]`. Thus
\[
             \gamma_a,\gamma_b\ge
                      \gamma_{\rm sp}:=\kappa-{2\delta L\over q}.
                                                               \tag{4}
\]
Whenever `gamma_sp>0`, summing (2) proves
\[
 \boxed{\quad
 S_3\le {2q\over\gamma_{\rm sp}}
       \left({a^2\over a^2+q}+{b^2\over b^2+q}\right)\Phi(A).
 \quad}                                                  \tag{5}
\]
If desired, the stronger individual positive coefficients in (3)
can be retained, giving `S_3<=Phi(A)(v_+/gamma_a+v_-/gamma_b)`.

## 2. Norm-only control of the coefficients

The already reviewed finite-dimensional bound
\[
                         L^2\le8\Phi(A)                 \tag{6}
\]
is equation (14), Section 3, of
`NOTE_2026-09-05_DIRECT_CROSS_COVARIANCE_NORMALIZATION.md`.
Its proof was checked directly: interpolation gives
`||A||^2<=2 beta_R(A)`, and zero-diagonal cube polarization gives
`beta_R(A)<=4 Phi(A)`.

Assume only
\[
                         \Phi(A)\le C n^{3/2},           \tag{7}
\]
where `C>0` is fixed. Then `L<=sqrt(8C)n^(3/4)` and
\[
                     \gamma_{\rm sp}
                                =\kappa-O_C(n^{-1/4}).
\]
It is uniformly positive for all sufficiently large n. Therefore (5)
is a LEADING source cubic-moment constraint under (7), without any
assumption that `L/sqrt(n)` remains bounded. In particular
\[
                         S_3\le(4/\kappa+o_C(1))q\Phi(A).
                                                               \tag{8}
\]
The more informative spectral factor in (5) is not discarded when
it is actually controlled. The uniform factor in (8) alone is not
claimed to improve the limiting constant `4/kappa`.

## 3. Retaining the actual diagonal-imbalance denominator

Let
\[
 v_+^*=q+\max_i T_{ii},\qquad
 v_-^*=q-\min_i T_{ii},\qquad
 V^*=v_+^*+v_-^*=2q+\operatorname{osc}_i T_{ii}.            \tag{9}
\]
These are the smallest admissible uniform variances for the two phases.
We now prove a norm-only bootstrap for this STRONGER denominator.

Trace zero gives `tr A_+=tr A_-=S_1/2`. Cauchy--Schwarz over at most
n eigenvalues and then averaging the phase variances give
\[
 S_2^\pm\ge {S_1^2\over4n},\qquad
 v_\pm^*\ge {2S_2^\pm\over n}
                         \ge {S_1^2\over2n^2}.          \tag{10}
\]
Consequently both exact coefficients in (2), with `v_pm=v_pm^*`,
are at least
\[
                  \gamma_*:=\kappa-{4\delta L n^2\over S_1^2}.
                                                               \tag{11}
\]
Whenever `gamma_*>0`, this proves the exact statement
\[
 \boxed{\qquad
  S_3\le {V^*\over\gamma_*}\Phi(A).
 \qquad}                                                 \tag{12}
\]

The source note's independently proved nuclear-norm inequality is
\[
                         S_1\ge{\kappa n^2q\over2\Phi(A)}.
                                                               \tag{13}
\]
Under (7), it gives
`S_1>=kappa sqrt(n)q/(2C)`. Combining this with (6),
\[
 {4\delta L n^2\over S_1^2}
 \le {16\delta C^2\sqrt{8C}\over\kappa^2}
                         {n^{7/4}\over q^2}
                         =O_C(n^{-1/4}).                 \tag{14}
\]
Thus `gamma_*=kappa-O_C(n^(-1/4))` is uniformly positive for large n,
and the actual diagonal version is
\[
 \boxed{\quad
 [2(n-1)+\operatorname{osc}_i(A|A|)_{ii}]\Phi(A)
            \ge[\kappa-O_C(n^{-1/4})]\operatorname{tr}|A|^3.
 \quad}                                                  \tag{15}
\]
Equivalently, with a constant depending only on C,
\[
 \boxed{\quad
 \Phi(A)\ge
 {\kappa\operatorname{tr}|A|^3
       \over2(n-1)+\operatorname{osc}_i(A|A|)_{ii}}
                       -O_C(n^{5/4}).
 \quad}                                                  \tag{16}
\]
For the equivalence, once `gamma_*>=kappa/2`, (12) bounds
`S_3/V^*<=2 Phi(A)/kappa`. The difference between the right-hand
coefficient `kappa` and `gamma_*` then costs at most
`O_C(n^(-1/4)) Phi(A)=O_C(n^(5/4))`.

## 4. Scope

Every exact original minimizer has a uniform bound of the form (7)
by CORE's elementary all-orders construction; so do the bounded-norm
near-minimizers used in the conditional Gaussian work. No new source
regularization is needed to apply (15)--(16). The source signing is not
replaced, and no fixed operator cap is silently introduced.

The error interpretation requiring a fixed operator cap in the original
source note remains correct for its unbootstrapped additive fourth-moment
estimate. The present supplement supplies an additional argument rather
than changing that previously reviewed assertion.

These constraints do not show that the diagonal oscillation is small,
that a cross SDP majorizer is scalar, or that source and cross singular
vectors have any favorable alignment. They remain same-order necessary
constraints, not an evaluated paired Gaussian upper or a convergence proof.
