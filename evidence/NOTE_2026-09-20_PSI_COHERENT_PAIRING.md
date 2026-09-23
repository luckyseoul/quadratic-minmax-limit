# Coherent pairing: identity, a sharp bound, conditional doubling

2026-09-20. Original limit OPEN. This does not claim existence.
Doubling alone does not (CORE §7 also needs the factor 3; a pure
doubling recurrence admits log-periodic sequences).

## 1. Identity

Let \(A\in\mathcal S_n\) and
\[
\Psi(A)=\Phi\begin{pmatrix} A & A \\ A & -A \end{pmatrix}.
\]
For \(x\in\{\pm1\}^n\) and \(T\subset[n]\) write
\(Q_T(x)=\sum_{i<j,\,i,j\in T}A_{ij}x_ix_j\). Then
\[
\Psi(A)=\max_{x,T}\bigl|2Q_A(x)-4Q_T(x)\bigr|.
\]

*Proof.* For \(y_i=s_ix_i\) and \(T=\{i:s_i=-1\}\), the coefficient of
\(A_{ij}x_ix_j\) in \(Q(x)-Q(y)+x^TAy\) is
\(1-s_is_j+s_i+s_j\). That equals \(-2\) when both indices lie in \(T\)
and \(+2\) otherwise. Summing gives \(2Q_A(x)-4Q_T(x)\). Every pair
\((x,y)\) arises this way. \(\square\)

In particular \(\Psi(A)\ge2\Phi(A)\), by \(T=\emptyset\).

## 2. The bound \(|Q-2Q_T|\le\Phi+n\) is false

It holds for every signing of order \(n\le7\) (exhaustive; minimum slack
\(\Phi+n-|Q-2Q_T|\) is \(3,2,3,0,1\)). It is sharp at \(n=6\): a plus
triangle and every other edge minus, \(Q(1)=-9=-\Phi\), \(Q_T=3\),
\(|Q-2Q_T|=15=\Phi+n\).

It fails at \(n=8\). Let \(T\) be a 4-set, put \(+1\) on every edge inside
\(T\) and \(-1\) on every other edge. Then \(Q(1)=-16\), \(Q_T=6\),
\(|Q-2Q_T|=28\), while a direct Boolean sweep gives \(\Phi=16\). Thus
\[
28=\lvert Q-2Q_T\rvert > \Phi+n=24.
\]
The same block pattern with \(t=\lfloor n/2\rfloor\) gives
\(|Q-2Q_T|=\binom{n}{2}\) and \(\Phi+n<\binom{n}{2}\) for every even
\(n\ge8\). So \(\Psi\) can exceed \(2\Phi\) by \(\Theta(n^2)\). The
\(|h|\to\infty\) coherent pairing is then the trivial \(O(n^2)\) regime,
not a doubling comparison.

Recorded minimizers still satisfy the false general bound
(\(\Psi/2\le\Phi+n\)), with equality at the order-10 minimizer
(\(\Psi=46\), \(\Phi=13\)). That is a finite observation, not a theorem.

## 3. Finite threshold, not the coherent limit

On the \(n=8\) counterexample, Monte Carlo of
\(B_h=\mathrm{sign}(G+hA)\) (40 draws, spectral-midpoint \(G\)) has its
smallest sample mean at \(h=0\): about \(42.2\), below
\(2\sqrt2\cdot16\approx45.25\). At \(h\ge3\) the mean has already
saturated at \(\Psi=56\), above the target. So for this signing the
useful threshold is the centered one, and sending \(h\to\infty\) makes
the bound worse.

On C5 the opposite happens: the sample means decrease toward
\(\Psi=12>2\sqrt2\cdot4\), and no tested \(h\) goes under the target.
Which \(h\) wins depends on the signing. A conditional theorem
“slack \(\Rightarrow m_{2n}\le2\sqrt2 m_n\) for large \(n\)” was checked
by `gpt-6-astra` (`PASS-WITH-NOTE`) and is not used: the hypothesis is
false.

## 4. What remains

The identity of §1 stands. Do not prove \(|Q-2Q_T|\le\Phi+n\).

An \(O(n)\) bound on \(\Psi-2\Phi\) for minimizers is too strong: together
with display (21) of the shifted-threshold note it would force
\(\alpha_n\to0\), against CORE §4. The coherent endpoint is pinned from
below,
\[
\limsup\Psi(A_n)/m_n\ge 2\sqrt2,
\]
Since \(m_n=\Theta(n^{3/2})\), \(\limsup\Psi/m_n\le 2\sqrt2\) is the same
statement as \(\Psi\le 2\sqrt2\, m_n+o(n^{3/2})\). Summability is stronger.
See
`NOTE_2026-09-23_PSI_MINIMIZER_EXCESS.md`. The finite inequality
\(\Psi\le 2m_n+2n\) is already false: an order-8 minimizer has
\(\Psi=40\). Limit OPEN.
