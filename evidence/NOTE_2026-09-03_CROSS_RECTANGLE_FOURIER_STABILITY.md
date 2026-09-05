# Cross-rectangle Fourier stability and the diagonal-payment gate

**Status:** proved analytic rigidity and proved obstructions to norm-only and
Gram-only cross-block closures.  This does **not** prove the multiplier-two
bound.  It isolates an exact statewise diagonal-payment inequality that a
successful use of global minimality would still have to establish.

This note addresses the coupled cross-rectangle family left by the exact
bisection recurrence in
`NOTE_2026-09-03_TWO_HALF_SELF_GLUING_OBSTRUCTION.md`.  It uses
\[
 \mu _4=\{1,i,-1,-i\}.
\]
The fourth-phase cube used elsewhere in the repository is
\(e^{i\pi/4}\mu _4\); the common rotation does not change any of the
Hermitian quadratic or cross-bilinear quantities below.

## 1. Exact phase and four-linear forms

Let \(C,R\in\{\pm1\}^{\ell\times r}\), and put
\[
 G=C+iR=(1+i)Q,
 \qquad Q\in\mu _4^{\ell\times r}.
\]
For diagonal sign matrices \(D,E\), set
\[
 X_{D,E}={1\over2}(C+DCE+DR-RE).                    \tag{1}
\]
Entrywise, \(X_{D,E}\) is \(C\) when the two endpoint labels agree,
\(R\) on the \((+,-)\) rectangle, and \(-R\) on the \((-,+)\)
rectangle.  Thus these are exactly the cross blocks in the bisection
half-cut family.

Write
\[
 \beta(X)=\max_{x\in\{\pm1\}^{\ell},\,y\in\{\pm1\}^{r}}
             |x^TXy|
\]
and
\[
 \Gamma(C,R)=\max_{D,E}\beta(X_{D,E}).              \tag{2}
\]
If
\(\rho(s)=\max(|\mathop{\rm Re}s|,|\mathop{\rm Im}s|)\), then
\[
 \Gamma(C,R)
 =\max_{z\in\mu _4^\ell,w\in\mu _4^r}\rho(z^*Gw). \tag{3}
\]
Indeed, take the endpoint phase to be \(1\) for label \(+1\) and
\(-i\) for label \(-1\).  The real part of the resulting entry is,
in the four label cells, respectively \(C,R,-R,C\).  Conversely every
fourth-phase vector is a real sign vector times such a label phase.  A
global multiplication by \(i\) interchanges the real and imaginary parts,
which accounts for \(\rho\).

There is also an exact four-linear form.  Put \(x'=Dx\) and \(y'=Ey\)
in (1).  As \(D,E,x,y\) vary, the four sign vectors are independent, so
\[
 \boxed{
 \Gamma(C,R)={1\over2}\max_{x,x',y,y'}
 \left|x^TCy+x'^TCy'+x'^TRy-x^TRy'\right|.}         \tag{4}
\]
It follows immediately that
\[
 \max(\beta(C),\beta(R))
 \le \Gamma(C,R)
 \le \beta(C)+\beta(R)
 \le \sqrt2\sqrt{\beta(C)^2+\beta(R)^2}.            \tag{5}
\]

Both losses in the last two bounds are real.  Let
\[
 K_2=\begin{pmatrix}1&1\\-1&1\end{pmatrix}.
\]
For \(C=R=K_2\otimes J_k\), one has
\(\beta(C)=\beta(R)=2k^2\).  Taking
\(D=E=\mathop{\rm diag}(I_k,-I_k)\) in (1) gives
\(X_{D,E}=J_{2k}\), and hence
\[
 \Gamma(C,R)=4k^2
 =\beta(C)+\beta(R)
 =\sqrt2\sqrt{\beta(C)^2+\beta(R)^2}.               \tag{6}
\]
Thus no strict improvement over the \(\ell_1\) merge can follow from the
two separate cut norms alone.

The tempting Pythagorean candidate
\[
 \Gamma(C,R)\le\sqrt{\beta(C)^2+\beta(R)^2}          \tag{7}
\]
also fails on a low-cut Hadamard example.  Let
\[
 H_4=\begin{pmatrix}
 1&1&1&1\\
 1&-1&1&-1\\
 1&1&-1&-1\\
 1&-1&-1&1
 \end{pmatrix},
 \qquad C=R=K_2\otimes H_4.
\]
Then \(C\) is Hadamard of order eight and \(\beta(C)=20\).  Indeed,
for a sign vector \(y\), write \(a_j=|(Cy)_j|/2\).  The \(a_j\) are
nonnegative integers with \(\sum a_j^2=16\).  Since
\((a_j-1)(a_j-2)\ge0\), summing gives \(\sum a_j\le10\).  Equality is
attained, for example, by
\[
 y=(-1,-1,-1,-1,-1,-1,-1,1),
\]
whose image has absolute values \(6,2,2,2,2,2,2,2\).
Also \(\beta(H_4)=8\), by Cauchy--Schwarz and the sign vector
\((-1,-1,-1,1)\).  The same \(D,E\) as above give
\(X_{D,E}=J_2\otimes H_4\), so
\[
 \Gamma(C,R)\ge32>20\sqrt2
 =\sqrt{\beta(C)^2+\beta(R)^2}.                      \tag{8}
\]

## 2. The universal row-averaging floor

For fixed \(w\), define
\[
 F_Q(w)=\sum_{a=1}^{\ell}
 \rho((Gw)_a).                                      \tag{9}
\]
The \(\ell_\infty\)-triangle inequality gives
\(\rho(z^*Gw)\le F_Q(w)\).  Rotating each coordinate of \(z\) by a
multiple of \(i\), one can put the larger Cartesian component of every
summand on the positive real axis, proving equality.  Hence
\[
 \boxed{\Gamma(C,R)=\max_{w\in\mu _4^r}F_Q(w).}      \tag{10}
\]

Let \(S_r,S'_r\) be independent sums of \(r\) independent Rademacher
signs and put
\[
 \eta_r=\mathbb E\max(|S_r|,|S'_r|).                \tag{11}
\]
For uniform \(w\), every row of \(Gw\) has the law
\(S_r+iS'_r\).  Therefore
\[
 \mathbb EF_Q=\ell\eta_r,
 \qquad
 \Gamma(C,R)\ge\max(\ell\eta_r,r\eta_\ell).        \tag{12}
\]
Here the second floor follows by applying the same argument to \(G^*\).
The planar central limit theorem and uniform integrability give
\[
 \eta_r=\left({2\over\sqrt\pi}+o(1)\right)\sqrt r. \tag{13}
\]
Thus a balanced square cross completion has the orientation-independent
floor
\((2/\sqrt\pi+o(1))r^{3/2}\), already at the critical square-root-two
scale relative to the ordinary Rademacher cut floor.

## 3. An exact nonzero Fourier coefficient

For every \(r\ge2\), define
\[
 a_r={\binom{2r-2}{r-1}\over 2^{2r-2}}.            \tag{14}
\]
Then for every pair of distinct columns \(j,k\), normalized expectation
on \(\mu _4^r\) gives
\[
 \boxed{
 \mathbb E_w\bigl[F_Q(w)\,\overline{w_j}w_k\bigr]
   =a_r\sum_{a=1}^{\ell}q_{aj}\overline{q_{ak}}.}   \tag{15}
\]

To prove (15), first consider one row and absorb its phases by setting
\(u_j=q_{aj}w_j\).  Write
\[
 (1+i)u_j=\alpha_j+i\beta_j,
\]
where all \(\alpha_j,\beta_j\) are independent Rademacher signs.  The
row function becomes
\[
 h_r(u)=\max\left(\left|\sum_j\alpha_j\right|,
                    \left|\sum_j\beta_j\right|\right),
\]
while
\[
 \overline{u_1}u_2
 ={\alpha_1\alpha_2+\beta_1\beta_2
   +i(\alpha_1\beta_2-\beta_1\alpha_2)\over2}.
\]
The imaginary terms average to zero, and exchangeability of the two
Rademacher families gives
\[
 \mathbb E[h_r(u)\overline{u_1}u_2]
 =\mathbb E[h_r(u)\alpha_1\alpha_2].                \tag{16}
\]
Condition on
\(A_0=\sum_{j=3}^r\alpha_j\) and
\(B=\sum_{j=1}^r\beta_j\), and let
\(\phi(t)=\mathbb E_B\max(|t|,|B|)\).  The remaining average is
\[
 {1\over4}\mathbb E_{A_0}
 \bigl(\phi(A_0+2)+\phi(A_0-2)-2\phi(A_0)\bigr).    \tag{17}
\]
On the common parity lattice \(t\equiv b\pmod 2\), which is the only lattice
sampled here, the step-two second difference of \(\max(|t|,b)\) is \(4\) at
\(t=b=0\), is \(2\) when \(|t|=b>0\), and is zero otherwise.  By symmetry,
(17) is therefore
\[
 \mathbb P(S_{r-2}=S_r)
 =\mathbb P(S_{2r-2}=0)
 ={\binom{2r-2}{r-1}\over2^{2r-2}}=a_r,
\]
where the two sums in the first probability are independent.  Restoring
the row phase contributes \(q_{aj}\overline{q_{ak}}\), and summing the
rows proves (15).  Stirling's formula gives
\[
 a_r={1\over\sqrt{\pi(r-1)}}(1+O(r^{-1})).          \tag{18}
\]

## 4. Quantitative stability and equality

Put
\[
 \Delta=\Gamma(C,R)-\ell\eta_r\ge0.
\]
Parseval and (15) imply
\[
 \mathop{\rm Var}F_Q
 \ge a_r^2\|Q^*Q-\ell I_r\|_F^2.                  \tag{19}
\]
On the other hand, if \(X=F_Q-\mathbb EF_Q\), then
\(-\ell\eta_r\le X\le\Delta\) and
\(\mathbb EX_+=\mathbb EX_-\le\Delta\).  Consequently
\[
 \mathop{\rm Var}F_Q
 \le \Delta(\ell\eta_r+\Delta).                   \tag{20}
\]
Combining (19)--(20) gives the stability estimate
\[
 \boxed{
 \|Q^*Q-\ell I_r\|_F
 \le {\sqrt{\Delta(\ell\eta_r+\Delta)}\over a_r}.} \tag{21}
\]
Thus equality in the row-averaging floor forces
\[
 Q^*Q=\ell I_r.                                    \tag{22}
\]
Equivalently,
\[
 C^TC+R^TR=2\ell I_r,
 \qquad C^TR=R^TC.                                 \tag{23}
\]
This is the rectangular complex-Hadamard analogue of the
equal-square/commutator rigidity in the directed-half analysis.

When \(\ell<r\), the Welch trace bound gives
\[
 \|Q^*Q-\ell I_r\|_F^2\ge r\ell(r-\ell).
\]
It follows that
\[
 \Delta\ge
 {\sqrt{(\ell\eta_r)^2+4a_r^2r\ell(r-\ell)}-\ell\eta_r\over2}. \tag{24}
\]

There is also a complete harmonic characterization.  If
\(\chi_\gamma(w)=\prod_jw_j^{\gamma_j}\) and
\(\widehat h_r(\gamma)=\mathbb E[h_r(w)\overline{\chi_\gamma(w)}]\),
then
\[
 \widehat F_Q(\gamma)
 =\widehat h_r(\gamma)\sum_{a=1}^{\ell}\chi_\gamma(q_a). \tag{25}
\]
Hence \(F_Q\) is constant if and only if the right side vanishes for
every nontrivial \(\gamma\).  Equation (22) is exactly the family of
constraints in (25) at characters \(w_j\overline{w_k}\), whose base
coefficient is the nonzero number \(a_r\).

In the balanced square case there is no nontrivial equality case at all.
Indeed, if \(\ell=r\) and (22) holds, then also \(QQ^*=rI_r\).  For any
row \(q_s\), take \(w=\overline{q_s}\).  Then
\[
 Qw=re_s,
 \qquad F_Q(w)=r.                                  \tag{26}
\]
If the averaging floor were attained, \(F_Q\) would be the constant
\(r\eta_r\), contradicting \(\eta_r>1\) for every \(r\ge2\).  Thus
\[
 \boxed{\Gamma(C,R)>r\eta_r\quad(r\ge2)}            \tag{27}
\]
for every balanced sign pair, although (21) alone does not give a useful
order-\(r^{3/2}\) gap.

## 5. Complex-Hadamard rigidity is not sufficient

There is an infinite analytic family showing that (22) may coexist with
the largest value allowed by the spectral estimate.  Let
\(k=4^d\), let \(W_k\) be the Sylvester Walsh matrix, and put
\[
 Q=C=R=K_2\otimes W_k,
 \qquad r=2k.
\]
Then \(Q^*Q=rI_r\).  Choose
\(D=E=\mathop{\rm diag}(I_k,-I_k)\).  Direct block multiplication gives
\[
 X_{D,E}=J_2\otimes W_k.
\]
The Walsh matrix has a quadratic bent sign vector, so
\(\beta(W_k)=k^{3/2}\); Cauchy--Schwarz gives the matching upper bound.
Therefore
\[
 \Gamma(C,R)\ge\beta(J_2\otimes W_k)=4k^{3/2}.
\]
But \(G^*G=2rI_r\), and hence
\[
 F_Q(w)\le\sum_a|(Gw)_a|
 \le\sqrt r\,\|Gw\|_2
 =\sqrt2\,r^{3/2}=4k^{3/2}.
\]
Thus
\[
 \boxed{\Gamma(C,R)=\sqrt2\,r^{3/2}}               \tag{28}
\]
on this infinite family.  Since
\(r\eta_r=(2/\sqrt\pi+o(1))r^{3/2}\), its ratio to the averaging floor
tends to \(\sqrt{\pi/2}\).  Therefore complex-Hadamard geometry is a
necessary near-floor condition, but it is very far from sufficient.  No
square Butson matrix attains the floor exactly for \(r\ge2\), and the
stability argument by itself does not exclude an asymptotically attaining
Butson sequence; that is a separate higher-Fourier discrepancy problem.

## 6. Exact diagonal-payment variable

Now let the full Hermitian completion be
\[
 H=\begin{pmatrix}H_L&G\\G^*&H_R\end{pmatrix}.
\]
For fourth-phase vectors \(z,w\), put
\[
 a(z)=z^*H_Lz,\qquad b(w)=w^*H_Rw,
 \qquad \rho(z,w)=\rho(z^*Gw).
\]
The two diagonal quantities are real.  Rotating \(w\) by
\(1,i,-1,-i\) leaves \(b(w)\) fixed and rotates the cross scalar.
Consequently
\[
 \nu_4(H)=\max_{z,w}
 \bigl(|a(z)+b(w)|+2\rho(z,w)\bigr).                \tag{29}
\]
With \(F_Q(w)=\max_z\rho(z,w)\), define the exact excess
\[
 P_H(w)=
 \max_z\bigl(|a(z)+b(w)|+2\rho(z,w)\bigr)-2F_Q(w). \tag{30}
\]
Then
\[
 0\le P_H(w)\le\max_z|a(z)+b(w)|
\]
and, if \(Z(w)=\mathop{\rm argmax}_z\rho(z,w)\),
\[
 P_H(w)\ge\max_{z\in Z(w)}|a(z)+b(w)|.             \tag{31}
\]
Most importantly,
\[
 \boxed{\nu_4(H)=\max_w\bigl(2F_Q(w)+P_H(w)\bigr).} \tag{32}
\]

In the intended completion \(H=A+iS\), let its real part \(A\) be a globally
minimal signing and put \(M=\Phi(A)\).  The required multiplier-two estimate
for this orientation is therefore exactly the statewise anti-correlation
condition
\[
 \boxed{
 P_H(w)\le2\bigl(\sqrt2M-F_Q(w)\bigr)
             +o_{\rm Dini}(n^{3/2})
 \quad\hbox{for every }w.}                          \tag{33}
\]
In particular \(\Gamma(C,R)>\sqrt2M+o(n^{3/2})\) is fatal regardless of
the diagonal blocks, since \(P_H\ge0\).

Global minimality alone supplies \(\Phi(A^F)\ge M\) for every flip; this
is the reverse direction from the upper control of \(P_H\) required in
(33).  The Fourier stability theorem controls only \(F_Q\), while
\(P_H\) contains the diagonal blocks and the choice of a cross-maximizing
row phase.  The maximally bad family (28) also shows that even exact
second-order Hadamard rigidity does not control the cross maximum.

Accordingly, (21) is a positive rigidity theorem, but not yet a
multiplier-two construction.  Equation (33) is the precise live conditional
sublemma.  A successful use of global minimality must add genuinely
higher-order information--for example an active-state or endpoint-slack
theorem that forces \(P_H(w)\) to decrease when \(F_Q(w)\) increases.
Any route using only the row-average floor and the Gram relations (23) is
ruled out by (28).

## Scope

The proved content is:

1. the exact cross-rectangle phase and four-linear formulas (3)--(4);
2. sharp obstruction examples for separate-norm and Pythagorean merges;
3. the exact Fourier coefficient (15) and stability estimate (21);
4. impossibility of exact balanced row-floor equality for \(r\ge2\);
5. an infinite Gram-perfect family attaining the spectral maximum; and
6. the exact diagonal-payment decomposition (32).

No finite-order census, small-prime computation, storage scan, or
self-gluing lift is used.  The still-open step is to construct, for a
globally minimal real signing, compatible skew blocks satisfying (33), or
to derive (33) from a new global-minimality coupling theorem.

The 2026-09-05 follow-up retains the actual diagonal blocks:
`NOTE_2026-09-05_DIAGONAL_PAYMENT_COMPATIBILITY.md` proves the exact
conjugate-pair transfer identity and all-alternating-cycle criterion for
arbitrary real payment functions, with a separate quadratic/sign
realization requirement. It also proves that diagonal skew cannot lower a
coherent `R=+C` or `R=-C` baseline. The exact mixed moments in
`NOTE_2026-09-05_ACTUAL_DIAGONAL_MIXED_MOMENTS.md` include the surviving
left/right interaction even at perfect Gram geometry. Neither follow-up
establishes (33).
