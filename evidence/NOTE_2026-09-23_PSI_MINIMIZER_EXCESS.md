# Coherent excess of a minimizer is pinned at the doubling constant

2026-09-23. Original limit OPEN. This does not claim existence, and it
does not claim the factor-2 half of CORE §7.

## 1. A deterministic comparison

Let \(A\in\mathcal S_n\) and
\[
\Psi(A)=\Phi\begin{pmatrix} A & A \\ A & -A \end{pmatrix}.
\]
The matrix
\[
K=\begin{pmatrix} A & A+I \\ A+I & -A \end{pmatrix}
\]
is Seidel of order \(2n\): its diagonal is zero, and every off-diagonal
entry is \(\pm1\). For \(x,y\in\{\pm1\}^n\),
\[
Q_K(x,y)=Q_A(x)-Q_A(y)+x^TAy+x\cdot y,
\]
so \(\lvert Q_K\rvert\le\Psi(A)+n\). Therefore
\[
m_{2n}\le\Psi(A)+n \tag{*}
\]
for every signing, and in particular for every exact minimizer.
On the order-8 minimizer below, \(\Phi(K)=\Psi=40\), so the \(+n\) room
is not always taken.

The same identity says \(\Psi(A)=2\max_T\Phi(A^{(T)})\), where
\(A^{(T)}\) is \(A\) with the principal block on \(T\) negated. The empty
block gives \(\Psi(A)\ge2\Phi(A)\).

Display (21) of the shifted-threshold note, with threshold
\(h=\sqrt{6\ln n}\), also yields
\(m_{2n}\le\Psi(A_n)+O(n^{16/11})\). That route is not needed for \((*)\).

## 2. An \(O(n)\) excess is incompatible with the lower bound

CORE §4 gives \(\liminf\alpha_n>0\). Suppose \(\gamma<2\sqrt2\) and
\(\Psi(A_n)\le\gamma\, m_n\) for every exact minimizer of every order
\(n\ge N\). Display \((*)\) gives
\[
\alpha_{2n}
\le\frac{\gamma}{2\sqrt2}\,\alpha_n+O(n^{-1/2}).
\]
The coefficient \(\lambda=\gamma/(2\sqrt2)\) is strictly less than 1.
Along \(n_j=2^j N\) the orbit tends to 0, including the boundary case
where the error ratio equals \(\lambda\). That contradicts the positive
liminf. Only positivity of the liminf is used.

Hence
\[
\limsup_{n\to\infty}\frac{\Psi(A_n)}{m_n}\ge 2\sqrt2.
\]
For every fixed \(C\), \(\Psi(A_n)\le 2m_n+Cn\) fails for infinitely many
\(n\). An \(O(n)\) excess is a stricter ceiling than \(2\sqrt2\, m_n\).

Because \(m_n=\Theta(n^{3/2})\),
\(\limsup\Psi(A_n)/m_n\le 2\sqrt2\) is the same statement as
\(\Psi(A_n)\le 2\sqrt2\, m_n+o(n^{3/2})\). By \((*)\) that upper bound is
sufficient for \(m_{2n}\le 2\sqrt2\, m_n+o(n^{3/2})\), not equivalent to
it. Summability of the resulting error, which CORE §7 needs, is stronger.
Neither statement is proved, and neither one closes existence.
On C5 and the Paley matrix of order 6 the sampled Boolean thresholds
decrease toward \(\Psi\) and do not undercut it.

## 3. The finite form \(\Psi\le 2m_n+2n\) is already false

The matrix
`evidence/psi_minimizer_excess_20260923/A8_psi40.txt`
is an exact minimizer of order 8: \(\Phi=10=m_8\) and \(\Psi=40\).
Thus
\[
\Psi=40>36=2m_8+2\cdot 8,
\qquad
\frac{\Psi}{m_8}=4>2\sqrt2.
\]
A census of all \(2^{21}\) switching representatives with first row
positive found 4200 classes with \(\Phi=10\). The largest value of
\(\lvert Q-2Q_T\rvert\) on those classes is 20, so the largest \(\Psi\)
is 40. The additive gap over \(2m_n+2n\) is exactly 4 at this order, and
it is attained. That census is finite evidence about order 8. It is not
the asymptotic theorem in §2.

Previously recorded minimizers remain above the same constant:
\(\Psi/m_n\) equals \(3, 3.6, 3.2, 10/3, 46/13, 3\) at orders
\(5,6,8,9,10,13\) for the saved matrices, and equals 4 for the matrix above.

## 4. Consult

Claude referee `deep_review` on 2026-09-23, model `claude-opus-5-5`,
effort `max`: **PASS-WITH-NOTE**, `do_not_branch`. The limsup in §2
stands. The note supplied the deterministic comparison \((*)\), which
replaces the Gaussian saturation argument. No new attack was opened.

## 5. What this does not do

The limit stays open. Factor 3 in CORE §7 is untouched. Section 2 does
not prove \(\Psi\le 2\sqrt2\, m_n+o(n^{3/2})\).
