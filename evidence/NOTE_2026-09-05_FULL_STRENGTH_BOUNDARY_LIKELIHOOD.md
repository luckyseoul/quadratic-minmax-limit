# Full-strength sign rounding: the exact boundary likelihood identity

2026-09-05. **Exact identity and a narrow actual-host counterexample;
the full-strength averaged comparison remains OPEN.**

Fix a complete signing (A). Retain the actual paired Gibbs prior

\[
\nu_\eta=\tfrac12(\mu_{A,\eta}\otimes\mu_{A,-\eta}
                       +\mu_{A,-\eta}\otimes\mu_{A,\eta}),
\qquad v=y\otimes x,
\]

and the canonical matrices (H,\mu,\Sigma_\rho=I+\rho H/\mu)
from the reviewed Gaussian-sign information note. Write (T=H/\mu),

\[
\ell(B)=\log\frac{Q_{\nu_\eta}(B)}{P_0(B)}
 =\log\mathbb E_{\nu_\eta}e^{\gamma v^T b}
                         -n^2\log\cosh\gamma,
\qquad J(\rho)=\mathbb E_{Q_{A,\rho}}\ell(B).
\]

Here (b=\operatorname{vec}B), and the host and all temperatures are
fixed while varying \(\rho\). The exact planted-channel identity is

\[
\mathbb E_{Q_{A,\rho}}F_{A,B}(t)
 =2a_A(\eta)+n^2\log\cosh\gamma+J(\rho),
\quad
J(\rho)=D(Q_{A,\rho}\Vert P_0)
          -D(Q_{A,\rho}\Vert Q_{\nu_\eta}).                 \tag{1}
\]

Only discrete relative entropies occur; both reference sign laws have
full support, even when \(\Sigma_1\) is singular.

## 1. Exact derivative, including the endpoint limit

For distinct cross coordinates \(e,f\), and an assignment \(b_{-ef}\)
of all other cross signs, define

\[
\Delta_{ef}\ell(b_{-ef})
 =\ell(++\,;b_{-ef})+\ell(--\,;b_{-ef})
   -\ell(+-\,;b_{-ef})-\ell(-+\,;b_{-ef}).
\]

If \(G\sim N(0,\Sigma_\rho)\) and \(0\le\rho<1\), then

\[
\boxed{\quad
J'(\rho)=\frac1{2\pi}\sum_{e<f}
 \frac{T_{ef}}{\sqrt{1-\rho^2T_{ef}^2}}
 \mathbb E\!\left[
   \Delta_{ef}\ell(\operatorname{sign}G_{-ef})
                 \mid G_e=G_f=0\right].\quad}             \tag{2}
\]

The conditional law in (2) is the Gaussian conditional law, with the
other coordinates subsequently thresholded. It depends on the pair
\((e,f)\), on \(\rho\), and on the fixed host. It is not the
unconditioned sign law or the original paired Gibbs law.

For a direct proof, differentiate the positive-definite Gaussian density:
\(\partial_\rho\varphi_{\Sigma_\rho}
=\tfrac12\sum_{ij}T_{ij}\partial_{ij}\varphi_{\Sigma_\rho}\).
Integrate by parts against the bounded piecewise-constant function
\(\ell(\operatorname{sign}g)\). Its distributional mixed derivative
is \(\Delta_{ef}\ell(\operatorname{sign}g_{-ef})
\delta_0(g_e)\delta_0(g_f)\). The two symmetric off-diagonal terms
cancel the factor \(1/2\); the density of \((G_e,G_f)\) at zero is
\([2\pi\sqrt{1-\rho^2T_{ef}^2}]^{-1}\). This proves (2).
Equivalently the same computation follows by smoothing the finitely
many orthant indicators and passing to the limit.

The mixed difference retains the complete Gibbs posterior. More precisely,
let \(\langle\cdot\rangle_{h}\) be the law proportional to
\(e^{h^Tv}\nu_\eta\). Holding \(h_j=\gamma b_j\) for \(j\ne e,f\),

\[
\Delta_{ef}\ell(b_{-ef})
 =\int_{-\gamma}^{\gamma}\!\int_{-\gamma}^{\gamma}
     \operatorname{Cov}_{h_e=s,h_f=u}(v_e,v_f)\,ds\,du.
                                                               \tag{3}
\]

In particular, \(|\Delta_{ef}\ell|\le4\gamma^2\). Since
\(|T_{ef}|\le1\), the resulting bound on (2) is integrable over
\(0\le\rho<1\), using
\(\int_0^1 |T_{ef}|/\sqrt{1-\rho^2T_{ef}^2}\,d\rho
=\arcsin|T_{ef}|\). Also couple
\(G_\rho=\sqrt\rho G_1+\sqrt{1-\rho}W\), with independent
\(W\sim N(0,I)\). Every coordinate of \(G_1\) has variance one,
so its sign is almost surely the limiting sign. Bounded convergence gives

\[
 J(1)=J(0)+\int_0^1 J'(\rho)\,d\rho.                       \tag{4}
\]

Thus (2)--(4) are an exact full-strength comparison, not a determinant
bound at a singular covariance or an annealed replacement.

## 2. A coordinatewise sign premise fails for an actual minimizer

A sufficient premise for \(J'\le0\) would be

\[
 T_{ef}\Delta_{ef}\ell(b_{-ef})\le0
 \quad\hbox{for every pair and every boundary sign context}.     \tag{5}
\]

That premise is false even under the actual global-minimizer and
balanced-step hypotheses.

Take \(n=3\), \(A=J_3-I_3\), and \(t=1\), so
\(\eta=\gamma=\beta/\sqrt2\). Every order-three signing is a
switching and possibly a global negative of this (A). Switching
preserves both partition functions; global negation swaps them.
Consequently this (A) is a global half-product minimizer at every
\(\beta>0\).

At any positive temperature \(s\) used to construct the covariance,
put \(z=e^{4s}\). Its actual opposite-phase off-diagonal correlations are

\[
 u=\frac{z-1}{z+3},\qquad
 v=\frac{1-z}{1+3z},\qquad
 \alpha=u+v=\frac{2(z-1)^2}{(z+3)(1+3z)}>0.
\]

Its extremal eigenvalues are \(2,-1\), hence \(\mu=2+\alpha\).
Choose \(e=(1,1)\), \(f=(2,1)\), and fix the other signs by

\[
 B(e,f)=
 \begin{pmatrix}
 e&1&1\\ f&-1&-1\\ -1&1&1
 \end{pmatrix}.
\]

The shared-column entry is \(H_{ef}=-\alpha<0\). Define the exact
zero-temperature width

\[
 W(B)=\max_{x,y}|Q_A(x)-Q_A(y)+x^TBy|
     =\max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TBy|\bigr).
\]

The second equality uses the independent global reversal of (x).
For this (A), \(Q_A(x)=3\) for a uniform spin vector and \(-1\)
otherwise. This gives the following short exact certificate:

| \((e,f)\) | one-sided uniform (x), nonuniform (y) | nonuniform (x), uniform (y) | (W(B)) |
|---|---:|---:|---:|
| \((+,+)\) | 5 | 9 | 9 |
| \((-,-)\) | 9 | 9 | 9 |
| \((+,-)\) | 7 | 11 | 11 |
| \((-,+)\) | 7 | 7 | 9 |

For verification without a signing census, write the total entry sum
as \(S=e+f+1\), row sums as \((e+2,f-2,1)\), and column sums
as \((e+f-1,1,1)\). A nonuniform three-spin vector, up to its
global negative, is obtained by reversing exactly one coordinate.
Thus the middle two table columns are respectively
\(4+\max_j|S-2c_j|\) and \(4+\max_i|S-2r_i|\).
When both spin vectors have the same uniform/nonuniform type, their
internal difference is zero and the cross magnitude is at most nine.
In the last row the bound nine is attained because
\(B=(-1,1,-1)^T(1,-1,-1)\), with both factors nonuniform.
This proves every entry of the width column.

There are 64 spin pairs. The maximum-term bound and \(\cosh z\le e^{|z|}\)
therefore give

\[
 \eta W(B)-\log128\le F_{A,B}(1)\le\eta W(B).
\]

The likelihood differs from this pressure by a constant independent
of (B). Hence the table implies

\[
 \Delta_{ef}\ell
 \le-2\eta+2\log128<0\qquad(\eta>\log128).
\]

It follows that \(T_{ef}\Delta_{ef}\ell>0\), disproving (5).
For every \(\rho<1\), the Gaussian covariance and its conditional
covariance after setting \(G_e=G_f=0\) are positive definite. This
specified boundary sign context therefore has strictly positive
conditional probability at every such \(\rho\).

This is **not** a counterexample to the averaged derivative or the
desired asymptotic finite-step bound. In fact, the particular context
above disappears from the singular endpoint boundary support: for this
order-three host, every matrix in the range of \(\Sigma_1\) has all
row and column sums equal. At \(G_{11}=G_{21}=0\), the first row of
the chosen context has positive sum and the second has negative sum,
which is impossible. This explicitly shows why the changing boundary
support and the weighted average in (2) cannot be discarded.

Indeed, on the eigenspaces of (A) with eigenvalues (2,-1\), the
mixed tensor sectors have covariance eigenvalue zero, while the two
same-sector eigenvalues are \(3(2-\alpha)/(2+\alpha)>0\) and
\(3(1+\alpha)/(2+\alpha)>0\), since \(0<\alpha<2/3\).
Thus the covariance range is exactly
\(\operatorname{span}(\mathbf1\mathbf1^T)
 \oplus(\mathbf1^\perp\otimes\mathbf1^\perp)\), the claimed
equal-row-and-column-sum space.

## 3. The precise remaining comparison

The actual finite-step target at a half-product minimizer is

\[
J(1)\le2[a_A(\beta)-a_A(\eta)]-n^2\log\cosh\gamma+o(n).
\]

Equations (2)--(4) reduce this to a weighted, pair-dependent conditional
posterior-covariance integral. The original negative prior contraction
\(\operatorname{tr}(\mathcal K H)<0\) does not evaluate that integral.
The counterexample above only rules out the coordinatewise premise (5).
It does not rule out cancellations in (2), control using its actual
boundary support, a direct likelihood comparison, or a selected outcome.
No computation or signing census was run for this note.
