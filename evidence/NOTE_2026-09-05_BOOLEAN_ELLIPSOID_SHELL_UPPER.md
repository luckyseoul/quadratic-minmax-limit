# Boolean remainder in a Gaussian ellipsoid shell upper bound

2026-09-05. This is an exact finite-dimensional Gaussian upper theorem.
It strengthens an ellipsoid shell bound by retaining an explicit expected
Boolean-coordinate penalty in the completion-square argument. The final scalar
evaluation improves the bound but does not establish the sharp original
comparison or original convergence.

Throughout this note, `kappa = 2/pi`. This is the Gaussian absolute-mean
constant, not the logarithmic tensor-lift constant used in some other
notes. All matrices and Gaussian vectors are real.

## 1. Exact Boolean remainder theorem

Let `S` be a nonempty subset of `{-1,1}^N`. Let `P` be a deterministic
symmetric positive-definite matrix, and suppose
\[
                       z^TPz=q\quad(z\in S).             \tag{1}
\]
Let `E` be a deterministic nonnegative diagonal matrix with `E <= P`
in positive-semidefinite order. Let `g` be centered Gaussian with
arbitrary positive-semidefinite covariance `C`, allowing singular `C`.
Put
\[
 T=\operatorname{tr}(CP^{-1}),\qquad
 r_i=(P^{-1}CP^{-1})_{ii},\qquad e_i=E_{ii},
\]
\[
 A_0=q-\operatorname{tr}E,\qquad
 B_0=T-\sum_i e_i r_i
    =\operatorname{tr}\bigl(CP^{-1}(P-E)P^{-1}\bigr).     \tag{2}
\]
Then `A_0,B_0 >= 0`, and
\[
 \boxed{\quad
 \mathbb E\max_{z\in S}g^Tz
 \le \sqrt{A_0B_0}+\sqrt\kappa\sum_i e_i\sqrt{r_i}.
 \quad}                                                  \tag{3}
\]
In particular, the following simpler consequence is valid:
\[
 \boxed{\quad
 \mathbb E\max_{z\in S}g^Tz
 \le \sqrt{\bigl(q-(1-\kappa)\operatorname{tr}E\bigr)
                 \operatorname{tr}(CP^{-1})}.
 \quad}                                                  \tag{4}
\]

### Proof

Fix any real `tau>0`, and write `w=P^{-1}g`. Completion of the square
and (1) give, for every `z in S`,
\[
 g^Tz={\tau q\over2}+{g^TP^{-1}g\over2\tau}
       -{\tau\over2}(z-w/\tau)^TP(z-w/\tau).             \tag{5}
\]
Because `P >= E` and `z_i` is a sign,
\[
 (z-w/\tau)^TP(z-w/\tau)
 \ge\sum_i e_i(z_i-w_i/\tau)^2
 \ge\sum_i e_i(|w_i|/\tau-1)^2.                          \tag{6}
\]
The last expression no longer depends on `z`; thus it may be subtracted
after maximizing over the shell. Every `w_i` is centered Gaussian with
variance `r_i`, so
\[
 \mathbb E|w_i|=\sqrt{\kappa r_i},\qquad
 \mathbb Ew_i^2=r_i.
\]
Taking expectations in (5)--(6) yields
\[
 \mathbb E\max_{z\in S}g^Tz
 \le {\tau A_0\over2}+{B_0\over2\tau}
       +\sqrt\kappa\sum_i e_i\sqrt{r_i}.                 \tag{7}
\]
For any Boolean `z in S`, `q=z^TPz >= z^TEz=tr E`, proving `A_0>=0`.
The trace expression in (2) proves `B_0>=0`. The infimum over `tau>0`
of the first two terms in (7) is `sqrt(A_0 B_0)`. This remains true
if either quantity is zero, by the appropriate one-sided limit.
Equation (3) follows. No independence among the coordinates is used.

For completeness, for each `s>=0`,
\[
                 1+s^2-2\sqrt\kappa s\ge1-\kappa.
\]
Applying this to `s=sqrt(r_i)/tau` in the expectation of (6) gives
\[
 \mathbb E\max_{z\in S}g^Tz
 \le{\tau\over2}\bigl(q-(1-\kappa)\operatorname{tr}E\bigr)
       +{T\over2\tau}.
\]
Optimizing `tau` proves (4). Its first factor is nonnegative, indeed
at least `kappa q`, since `0<=tr E<=q`. If `T=0`, then `g=0` almost
surely and both bounds hold directly. This also proves that (3),
which optimizes the exact expression in (7), is at least as strong
as the estimate (4) obtained by weakening (7) first.

### Two limiting checks

If `E=0`, (3) reduces to the usual bound `sqrt(q tr(CP^{-1}))`.
If `P=E` is positive diagonal, then `q=tr P` on the entire cube,
`A_0=B_0=0`, and (3) gives exactly
\[
                  \sum_i\sqrt{\kappa C_{ii}}.
\]
That is the unrestricted cube width, including for dependent Gaussian
coordinates. Thus the Boolean correction has the correct cube limit.

## 2. Diagonal-affine shell specialization

Let `D` be positive diagonal and `H` symmetric, with
\[
                         D-H\succeq0,\qquad D+H\succeq0.
                                                               \tag{8}
\]
Assume `z^THz=h` throughout the shell. Write
\[
 d_0=\operatorname{tr}D,\quad
 T_0=\operatorname{tr}(CD^{-1}),\quad
 u={h\over d_0},\quad
 v={\operatorname{tr}(CD^{-1}HD^{-1})\over T_0}.          \tag{9}
\]
If `T_0=0`, the width is zero; otherwise all quantities in (9) are
defined and `u,v in [-1,1]`. For `u`, use (8) on Boolean vectors.
For `v`, put `J=D^{-1/2}HD^{-1/2}`, so `-I<=J<=I`, and take the
trace against the positive-semidefinite matrix `D^{-1/2}CD^{-1/2}`.

For each `-1<eta<1`, choose
\[
 P_\eta=D-\eta H,\qquad E_\eta=(1-|\eta|)D,
 \qquad q_\eta=d_0(1-\eta u).                           \tag{10}
\]
Then `P_eta>0`, `0<=E_eta<=P_eta`, and the shell satisfies (1).
Scalar functional calculus on the contraction `J` gives
\[
 (I-\eta J)^{-1}\preceq{I+\eta J\over1-\eta^2}.
\]
Indeed the scalar difference at `lambda in [-1,1]` is
`eta^2(1-lambda^2)/((1-eta^2)(1-eta lambda)) >= 0`.
Consequently
\[
             \operatorname{tr}(CP_\eta^{-1})
                    \le T_0{1+\eta v\over1-\eta^2}.     \tag{11}
\]
Applying (4) proves the explicitly optimized affine bound
\[
 \boxed{\quad
 \mathbb E\max_{z\in S}g^Tz
 \le\sqrt{d_0T_0}\inf_{-1<\eta<1}
 \sqrt{\frac{[1-\eta u-(1-\kappa)(1-|\eta|)](1+\eta v)}
                  {1-\eta^2}}.
 \quad}                                                  \tag{12}
\]
The case `eta=0` is the cube bound `sqrt(kappa d_0 T_0)`.
In contrast to a putative uniform multiplication of an ellipsoid
bound by `sqrt(kappa)`, (12) keeps the diagonal remainder available
at each particular metric. No uniform multiplier is assumed.

## 3. Exact minimization of the scalar affine expression

For `u,v in [-1,1]`, define
\[
 d=1-\kappa-u,\qquad b=d+\kappa v,
\]
and define the nonnegative-eta squared optimum by
\[
 \mathcal F_+(u,v)=\inf_{0\le\eta<1}
                 {(\kappa+d\eta)(1+v\eta)\over1-\eta^2}.
                                                               \tag{13}
\]
Then the complete formula is
\[
 \boxed{\quad
 \mathcal F_+(u,v)=
 \begin{cases}
 \kappa,&b\ge0,\\[2mm]
 \displaystyle{\kappa-dv+
       \sqrt{(1-u)(2\kappa-1+u)(1-v^2)}\over2},&b<0.
 \end{cases}
 \quad}                                                  \tag{14}
\]
In the second branch put
\[
 \Delta=(1-u)(2\kappa-1+u)(1-v^2).
\]
An optimizer, or endpoint limiting optimizer, is
\[
              \eta_*={-b\over\kappa+dv+\sqrt\Delta}.     \tag{15}
\]
The full squared factor in (12) is
\[
           \min\{\mathcal F_+(u,v),\mathcal F_+(-u,-v)\}. \tag{16}
\]

### Proof and endpoints

Write `a=kappa`, `c=dv`. Subtracting `kappa` from the function in
(13) gives
\[
                  {\eta[b+(\kappa+dv)\eta]\over1-\eta^2}.
\]
If `b>=0`, the bracket is nonnegative on `[0,1]`: its values at
the two endpoints are `b>=0` and
\[
                b+\kappa+dv=(1-u)(1+v)\ge0.
\]
The minimum is therefore `kappa`, attained at zero.

If `b<0`, then `u>1-kappa+kappa v >= 1-2kappa`. In particular
`Delta>=0`; the square root in (14) is real on this entire branch,
not merely when `u>=0`. Moreover
\[
 (\kappa+dv)^2-b^2
   =(\kappa-d)(\kappa+d)(1-v^2)=\Delta,
\]
and `kappa+dv >= -b > 0`. Differentiating (13) gives a numerator
\[
                     b+2(\kappa+dv)\eta+b\eta^2.
\]
Its first zero in `(0,1]` is (15), and the function decreases up
to that zero and increases thereafter when the zero is interior.
Substitution gives the second line of (14).

When `Delta=0` on the `b<0` branch, `eta_*=1` and the same value
is the limit as `eta` increases to one. This occurs only at the
relevant boundary `u=1` or `v=-1`; the formulas include their
intersection. The infimum need not be attained by a positive-definite
metric, but (12) holds for each `eta<1` and therefore for its limit.
Finally, negative `eta` is positive `-eta` with `(u,v)` replaced
by `(-u,-v)`, which proves (16).

For the favorable cross domain `u>=0,v<=0`, the negative branch
has `b_-=1-kappa+u-kappa v>0` and hence never improves `kappa`.
The precise positive-branch improvement criterion is
\[
                         u-\kappa v>1-\kappa.           \tag{17}
\]

## 4. A faithful form of the stronger remainder

Equation (12) uses the simpler bound (4). The stronger (3) remains
available at every metric (10), with no additional assumptions.
One useful intermediate form retains two covariance traces. Define
\[
 T_\eta=\operatorname{tr}(CP_\eta^{-1}),\qquad
 R_\eta=\operatorname{tr}(DP_\eta^{-1}CP_\eta^{-1}).
\]
Cauchy--Schwarz on the final sum in (3) gives
\[
 \mathbb E\max_{z\in S}g^Tz\le\sqrt{d_0}\left[
 \sqrt{(|\eta|-\eta u)
         [T_\eta-(1-|\eta|)R_\eta]}
       +\sqrt\kappa(1-|\eta|)\sqrt{R_\eta}\right].        \tag{18}
\]
The expression in square brackets inside the first square root is
nonnegative because it equals
`tr(C P_eta^{-1}(P_eta-E_eta)P_eta^{-1})`.
No separate upper substitution for `R_eta` is justified without
checking the combined expression: the same trace appears with a
negative sign in one term and a positive sign in the other.
The unsimplified sum in (3) is at least as strong as (18).

## 5. Actual cross covariance and the remaining numerical gap

For the joint-shell reference field of Sections 1--3 in
`original_mo_mu_joint_shell_extension.md` at zero internal shell values,
write
\[
 H_B=\begin{pmatrix}0&B\\B^T&0\end{pmatrix},\qquad
 C=a_0I_{2n}-tcH_B,
\]
where `c=x^TBy`. In the unbiased, mu-normalized construction,
\[
       t={\kappa\over\mu},\qquad
       a_0=n\left(1+{\kappa\over\mu}\right),\qquad
       \mu\ge n-1.                                     \tag{19}
\]
These expressions include the independent Hermite cushion `1-kappa`.
If `D=ell I` is a valid diagonal majorizer of `H_B`, then
\[
 u={c\over n\ell},\qquad
 v=-\rho u,\qquad
 \rho={tn^2\over a_0}={\kappa n\over\mu+\kappa},\qquad
 \sqrt{d_0T_0}=2n\sqrt{a_0}.                            \tag{20}
\]
Here `tr H_B=0` and `tr H_B^2=2n^2` were used. The limiting
strength is at most `kappa+o(1)`, not one; discarding the Hermite
cushion would incorrectly strengthen the evaluated comparison.
The matrix `C` in this paragraph must actually be positive semidefinite
to invoke the Gaussian theorem. If it is indefinite, the theorem is
instead applied to the positive-semidefinite rank-four repaired
covariance of the cited note. In that case (20) is only the uncorrected
trace algebra until the repair contribution has been explicitly
accounted for; no Gaussian with the indefinite covariance is asserted.

As a scalar diagnostic, let
\[
 K_G={\pi\over2\log(1+\sqrt2)},\qquad
 u=K_G^{-1},\qquad v=-\kappa K_G^{-1}.
\]
Equation (14) gives
\[
 \begin{split}
 u&=0.561099852339180\ldots,\\
 v&=-0.357207260271652\ldots,\\
 \eta_*&=0.334096674476977\ldots,\\
 \mathcal F_+(u,v)&=0.565603376492680\ldots,\\
 \sqrt{\mathcal F_+(u,v)}&=0.752066071893075\ldots.
 \end{split}                                            \tag{21}
\]
This beats the central cube factor `sqrt(2/pi)=0.79788456...`.
Nevertheless the resulting width coefficient `1.50413214...` when
`a_0/n -> 1` remains above `sqrt(2)=1.41421356...`.
Thus this evaluated point does not by itself reach the desired
leading comparison at original norm constant one half.

The scalar point in (21) is not asserted to arise from every actual
source or to be a worst shell. In particular a diagonal-majorizer
trace budget does not imply that a uniform majorizer with that
budget exists. Equations (3), (12), and (18) are valid upper bounds
at valid metrics; converting them into an all-source, all-shell
sharp leading comparison remains open.
