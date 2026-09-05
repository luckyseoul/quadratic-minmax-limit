# Cubic frame alignment crosses the weak-Dirac template threshold

2026-09-05. Analytic additive theorem. The scalar-optimal finite-template
completion satisfies `Gamma(C)>283/200` already when `||C||op>=12/5`.
The new ingredient is a quantitative obstruction to anti-alignment of
the cubic tensor features of the actual optimal frames. No rank bound,
sign-matrix realization, or source compatibility is assumed or inferred.

This strengthens a bound on the completion UPPER CERTIFICATE Gamma.
It does not by itself lower-bound an actual constructed Boolean norm.
Its relevance to the previously isolated weak-Dirac diagnostic is
spelled out separately in Section 6.

## 1. Actual scalar-optimal template and reviewed baseline identities

Let C be a real p by p matrix, write q=||C||op>1, and assume that its
actual real vector SDP has value tau(C)=pq. Define

\[
 \Gamma(C)=\max_{a,b\in[-1,1]^p}\left[
 {a^TCb\over p}
 +\sqrt{\left(1-{\|a\|^2\over p}\right)
        \left(1-{\|b\|^2\over p}\right)}\right].           \tag{1}
\]

Choose optimal unit-row frames U,V. Equality in the operator bound
gives `CV=qU`, `C^TU=qV`, and hence

\[
 M=U^TU=V^TV,\qquad\operatorname{tr}M=p,
 \qquad \mu={\lambda_{\max}(M)\over p},\quad
 s={\operatorname{tr}(M^2)\over p^2},\quad
                         0<s\le\mu\le1.                 \tag{2}
\]

The complete derivations of the baseline identities used below are in
*A quantitative completion bound for scalar-optimal finite templates*,
source `original_mo_scalar_template_gamma_bound.md`, final SHA256
`bd5997203c52895744a078048e206241996c46ef485e8975d7955b73be41f1c6`.
They are summarized with their exact hypotheses here.

For a unit top eigenvector z of M, the cube pair Uz,Vz has equal
normalized squared norm mu and bilinear value q mu. Therefore, with
`x=Gamma(C)-1`,

\[
                  \Gamma(C)\ge1+(q-1)\mu,
                  \qquad\mu\le{x\over q-1}.              \tag{3}
\]

Let f(t)=clip(t,-1,1), let G be standard Gaussian, and put

\[
 P=\Pr\{|G|\le1\},\quad \phi={e^{-1/2}\over\sqrt{2\pi}},
 \quad v=E[f(G)^2]=1-2\phi,
\]
\[
 R=\operatorname{Var}(f(G)^2)-2(P-2\phi)^2,
 \qquad D(q)=q(2P^2-v)-v.                                \tag{4}
\]

The nonnegative number R is the squared Gaussian L2 norm of the even
Hermite remainder after deleting degrees zero and two from f^2.
For a= f(Ug), b=f(Vg), with the same standard Gaussian g, the exact
matched-frame identity cancels the entire degree-two fluctuation of
`||a||^2/p-||b||^2/p`. Even Hermite covariance is nonnegative and
bounded by R times squared correlation. Consequently

\[
 E\sqrt{\left(1-{\|a\|^2\over p}\right)
        \left(1-{\|b\|^2\over p}\right)}
                        \ge1-v-\sqrt{Rs/2}.              \tag{5}
\]

The baseline proof bounded every nonlinear odd feature by the crude
operator lower -q. The next section improves its cubic feature using
the same ACTUAL optimal frames.

## 2. Cubic tensor features cannot be fully anti-aligned

Let T=C/q, a contraction, and let U_3,V_3 have rows `u_i tensor^3`
and `v_j tensor^3`. All their row norms are one. Put

\[
                   j_3={\langle U_3,TV_3\rangle_F\over p}.
\]

Then

\[
                         \boxed{j_3\ge-1+{2s^2\over\mu}.} \tag{6}
\]

To prove this, set F_3=U_3+TV_3. Since `U^TT=V^T`,

\[
 U^TF_3=U^TU_3+V^TV_3
           =\sum_i u_i^{\otimes4}+\sum_j v_j^{\otimes4},  \tag{7}
\]

where the last symmetric tensors are flattened into matrices with
one tensor coordinate in the row index and three in the column index.
Flattening preserves Frobenius norm. The squared norm of (7) equals

\[
 \sum_{i,i'}\langle u_i,u_{i'}\rangle^4
 +\sum_{j,j'}\langle v_j,v_{j'}\rangle^4
 +2\sum_{i,j}\langle u_i,v_j\rangle^4.                    \tag{8}
\]

For EACH of these three p by p arrays, its sum of squared inner
products is `tr(M^2)=p^2s`, by the common frame matrix in (2).
Cauchy--Schwarz applied to its p^2 nonnegative squared entries
therefore gives a sum of fourth powers at least p^2 s^2. Thus

\[
 \|U^TF_3\|_F^2\ge4p^2s^2,
 \qquad
 \|F_3\|_F^2\ge{4p^2s^2\over\|U\|_{\rm op}^2}
                         ={4ps^2\over\mu}.               \tag{9}
\]

On the other hand, T is a contraction, so

\[
 \|F_3\|_F^2
 =p+\|TV_3\|_F^2+2\langle U_3,TV_3\rangle_F
                          \le2p+2pj_3.                   \tag{10}
\]

Combining (9)--(10) proves (6). It uses neither an actual-rank bound
nor the replacement of actual frames by an arbitrary moment measure.

## 3. Retaining this positive correction in Gaussian clipping

Expand the bounded odd f in normalized probabilists' Hermite
polynomials. Write its coefficients as b_1,b_3,b_5,...; then
`sum b_k^2=v`, and `b_1=E[G f(G)]=P`. The cubic coefficient is

\[
 b_3={E[(G^3-3G)f(G)]\over\sqrt6}
                         =-{2\phi\over\sqrt6},\qquad
 b_3^2={2\phi^2\over3}.                                  \tag{11}
\]

For example `E[G^3 f(G)]=3P-2phi` follows by adding the truncated
fourth moment `3P-8phi` to the cubic absolute tail `6phi`; subtracting
3P proves (11).

For every odd k, Gaussian Hermite orthogonality gives the contribution
`b_k^2 <U_k,CV_k>_F`, where U_k,V_k have unit tensor-power rows.
Each normalized contraction objective is at least -1, the linear
one equals 1, and the cubic one obeys (6). Consequently

\[
 {E[a^TCb]\over p}
 \ge q\left[2P^2-v+{2b_3^2s^2\over\mu}\right]
 =q(2P^2-v)+{4q\phi^2\over3\mu}s^2.                      \tag{12}
\]

This identity and bound are first valid for finite Hermite sums and
then pass to Gaussian L2 limits, since p is finite and f is bounded.
No assertion of independent tensor features is needed.

The maximum Gamma dominates the expected admissible objective.
Combining (5) and (12) proves the strengthened finite-template bound

\[
 \boxed{\displaystyle
 x=\Gamma(C)-1
       \ge D(q)+{4q\phi^2\over3\mu}s^2-\sqrt{Rs/2}.}     \tag{13}
\]

The frame quantities s and mu in (13) remain those of the actual
optimal frames. In particular (3) is available simultaneously.

## 4. Exact threshold q>=12/5

The companion proof's exact finite-sum enclosures give

\[
 P>P_0={68268\over100000},\qquad
 \phi>\phi_0={24197\over100000},\qquad
 R<R_0={844\over10000},\qquad 2P^2-v>0.                  \tag{14}
\]

Assume for contradiction that q>=12/5 and
`Gamma(C)<=283/200`, so `x<=83/200=0.415`. Equations (3) and (14)
give the following exact rational lower bounds:

\[
 D(q)\ge D(12/5)
 >{12\over5}(2P_0^2-1+2\phi_0)-1+2\phi_0
                 =0.48244551552>0.4824,                  \tag{15}
\]
\[
 {4q\phi^2\over3\mu}
 \ge{4q(q-1)\phi^2\over3(83/200)}
 >{896\over83}\phi_0^2>0.63.                            \tag{16}
\]

Every denominator is positive: q>1, mu>0, and the cap 83/200 is
positive. The function q(q-1) is increasing on q>=12/5.

For positive R_0 and b_0=0.63, elementary one-variable differentiation
gives

\[
 \left(\max_{s\ge0}
         [\sqrt{R_0s/2}-b_0s^2]\right)^3
                         ={27R_0^2\over1024b_0}
                         <(0.067)^3.                    \tag{17}
\]

To check the maximization, put t=sqrt(s). The function is
`sqrt(R_0/2)t-b_0t^4`, whose unique positive critical point satisfies
`4b_0t^3=sqrt(R_0/2)`. Substitution gives exactly (17).
The strict rational comparison in (17), after multiplication by
1024b_0, is just

\[
             0.19233072<0.19402822656.                   \tag{18}
\]

Using (15)--(17) in (13) now forces

\[
              x>0.4824-0.067=0.4154>0.415,
\]

contradicting the assumed cap. We have proved

\[
 \boxed{\displaystyle
 q\ge12/5\quad\Longrightarrow\quad
                  \Gamma(C)>283/200>\sqrt2.}             \tag{19}
\]

All decimal-looking quantities in this section denote exact
terminating rationals. Equations (15), (16), and (18) are new finite
rational comparisons on top of the already verified Gaussian
enclosures; no floating-point evaluation is required for the proof.

## 5. Why this is stronger than an effective-rank observation

The baseline completion estimate alone, under a small Gamma cap,
can force s to be bounded below. That only bounds `1/s`, an effective
frame rank. It does not bound the actual rank of U or V and does not
authorize insertion of a finite-rank Grothendieck theorem.

Instead, (6) uses the actual common frame matrix to force a positive
cubic correction whenever s is nonzero. The same s that enlarges
the completion fluctuation penalty also enlarges this correction.
The optimization (17) controls the two terms jointly, without
asserting any rank truncation or rank-sensitive external result.

## 6. Crossing the previously isolated weak-Dirac barrier

For the restricted scalar active-cross diagnostic, STIPULATE leading
Boolean energy f=sqrt(2), and set

\[
 u={\sqrt2\over q},\qquad m=q^{-2}=u^2/2,
 \qquad \nu=\delta_m,
\]

and take its weak-feedback limit chi=0. At this stipulated active
normalization q>=sqrt(2), so 0<u<=1 and 0<m<=1/2. In particular u is
now variable; it is not frozen at the older diagnostic's Krivine
endpoint. The assumption Gamma(C)<=sqrt(2) used below is an upper
certificate condition. It does NOT imply saturation of an actual
Boolean norm or supply the separately stipulated value f=sqrt(2).

For completeness, substitute chi=0 and nu=delta_m into the GENERAL
two-trace formula (4) of *Evaluated scalar-diagonal diagnostic: strong
and weak feedback*, source `original_mo_scalar_moment_feedback_diagnostic.md`,
SHA256 `cc3869aa35b88ae50425c29cb78e3d4ced9b73e24731f54556fbd0b39fab1e9c`.
With kappa=2/pi, this gives for 0<=eta<1

\[
 A_\eta(m)={1+(\eta^2-2\eta)m\over(1-\eta^2m)^2},
 \qquad
 B_\eta(m)={1+\eta^2m\over(1-\eta^2m)^2},                 \tag{20a}
\]
\[
 U_\eta=\eta\sqrt{(1-u)A_\eta(m)}
             +(1-\eta)\sqrt\kappa\sqrt{B_\eta(m)}.
                                                                    \tag{20b}
\]

Because m<=1/2, the denominator stays separated from zero even at
eta=1. The squared first term in (20b) tends to
`(1-u)(1-m)/(1-m)^2`; the second term tends to zero. Therefore the
eta-up-to-one limit of the two-trace functional squared is

\[
                         {1-u\over1-u^2/2}.              \tag{20}
\]

This limiting substitution also respects the actual PSD covariance
repair. Equations (4a)--(4b) of the cited diagnostic do not require
its old numerical value of u: for bounded q the normalized repaired
resolvent traces obey

\[
 0\le\widetilde T_\eta-T_\eta
       \le{2+q^2\over n(1-\eta)},\qquad
 0\le\widetilde R_\eta-R_\eta
       \le{2+q^2\over n(1-\eta)^2}.                       \tag{20c}
\]

Evaluate the actual repaired functional at each FIXED eta<1 first.
Both trace errors vanish as n tends to infinity, and the bounded
continuous spectral integrands pass to nu=delta_m. Only afterwards
take eta to one in (20b). The fixed-eta inequality
`limsup_n inf_eta F_n(eta)<=limsup_n F_n(eta)` justifies the resulting
upper comparison. No uniform finite-n endpoint resolvent bound and
no PSD assumption on the unrepaired covariance are being used.

If its proposed finite-template certificate has Gamma(C)<=sqrt(2),
(19) requires q<12/5. Therefore

\[
 u>{5\sqrt2\over12}>2-\sqrt2.
\]

The last strict inequality is equivalent to `17sqrt(2)>24`, whose
square is the exact comparison 578>576. Equivalently, the q bound
has crossed `12/5<1+sqrt(2)`. Equation (20) is less than 1/2 precisely
when u>2-sqrt(2) in the indicated range. Thus this particular weak
Dirac moment diagnostic passes the target in the remaining
Gamma-certifiable range.

This is a conclusion about the stated diagnostic and the stated
template certificate. It does not prove that arbitrary actual B
have a Dirac bulk, a scalar dual, or Boolean norm bounded below by
Gamma. The general source/joint-shell comparison and the original
MO limit remain unresolved.
