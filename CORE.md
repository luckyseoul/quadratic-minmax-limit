# Original MO problem: verified core

The question is whether the sequence below converges. Convergence and the
value of a possible limit are not proved here. These arguments require no
optimality theorem for a special matrix family, no completion assumption,
and no archived proof-status ledger.

## 1. Definition

For an integer \(n\ge2\), let \(\mathcal S_n\) be the symmetric real
\(n\times n\) matrices with zero diagonal and off-diagonal entries in
\(\{-1,1\}\). Define
\[
Q_A(x)=\sum_{i<j}A_{ij}x_ix_j=\tfrac12x^TAx,\qquad
\Phi(A)=\max_{x\in\{-1,1\}^n}|Q_A(x)|,
\]
\[
m_n=\min_{A\in\mathcal S_n}\Phi(A),\qquad
\alpha_n=\frac{m_n}{n^{3/2}}.
\]
All minima and maxima are over finite sets. The original problem asks
whether \(\lim_{n\to\infty}\alpha_n\) exists; it does not prescribe its value.

## 2. Restriction and padding

For \(2\le n\le N\),
\[
m_n\le m_N,\qquad 0\le m_{n+1}-m_n\le n.
\]
Indeed, take any order-\(N\) matrix \(A\), fix signs on a set \(S\) of
\(n\) coordinates, and average over independent uniform signs outside
\(S\). Every term meeting the complement has mean zero, so
\(Q_{A[S]}(x)=\mathbb E[Q_A(X)\mid X_S=x]\). Consequently
\(\Phi(A[S])\le\Phi(A)\); minimize using an optimal \(A\).
Conversely, adjoining one vertex with arbitrary incident signs adds a
term of absolute value at most \(n\). Iteration gives, for integers \(d\ge0\),
\[
m_{n+d}\le m_n+dn+\frac{d(d-1)}2.
\]
It is \(m_n\), not necessarily \(\alpha_n\), that is nondecreasing.

## 3. An elementary all-orders upper bound

For every \(n\ge2\),
\[
m_n\le\sqrt{\log2}\,n^{3/2}.
\]
Choose the \(K=\binom n2\) upper-triangular entries independently and
uniformly in \(\{-1,1\}\). For each fixed \(x\), \(Q_A(x)\) is a sum of
\(K\) independent uniform signs. Since
\(\mathbb E e^{sQ_A(x)}=(\cosh s)^K\le e^{Ks^2/2}\), the exponential
Markov inequality, optimized at \(s=t/K\), gives
\[
\Pr\{|Q_A(x)|\ge t\}\le2e^{-t^2/(2K)}.
\]
A union bound over all \(2^n\) states yields
\[
\Pr\{\Phi(A)\ge t\}\le2^{n+1}e^{-t^2/(n(n-1))}.
\]
At \(t=\sqrt{\log2}\,n^{3/2}\), the right side is
\(2^{-1/(n-1)}<1\), so a signing below this threshold exists.

Together with Section 2 this also gives
\[
|\alpha_{n+1}-\alpha_n|
\le\frac{n}{(n+1)^{3/2}}
 +m_n\bigl(n^{-3/2}-(n+1)^{-3/2}\bigr)
=O(n^{-1/2}).
\]
Boundedness and vanishing consecutive gaps do not themselves imply
convergence.

They do show that the set of subsequential limits is the closed interval
between the lower and upper limits. Every value strictly between those
ends is crossed infinitely often; at each crossing a sequence term lies
within one consecutive gap of that value. The gaps vanish, giving the
required convergent subsequence. Boundedness supplies both endpoints.

## 4. A uniform Gaussian lower bound

For every \(A\in\mathcal S_n\),
\[
\Phi(A)\ge\frac{n\sqrt{n-1}}\pi.
\]
For \(n=2\), \(\Phi(A)=1\ge2/\pi\), so this case is immediate.
Assume henceforth \(n\ge3\). Let \(g\) be a vector of independent standard
normal variables, and put
\[
z^\pm=(I\pm A/\sqrt{n-1})g,\qquad x^\pm=\operatorname{sgn}(z^\pm).
\]
Every coordinate of \(z^\pm\) has variance two. For \(i\ne j\), its pair
correlations satisfy
\[
A_{ij}r_{ij}^\pm=u_{ij}\pm v,\qquad
u_{ij}=\frac{A_{ij}(A^2)_{ij}}{2(n-1)},\qquad
v=\frac1{\sqrt{n-1}}.
\]
Since \(|(A^2)_{ij}|\le n-2\),
\[
|u_{ij}\pm v|\le\frac{n-2}{2(n-1)}+\frac1{\sqrt{n-1}}
=1-\frac12\left(1-\frac1{\sqrt{n-1}}\right)^2<1.
\]
The Gaussian sign identity is
\(\mathbb E[\operatorname{sgn}Z_1\operatorname{sgn}Z_2]
=(2/\pi)\arcsin r\) for correlation \(r\). To see it, realize the pair
as projections of an isotropic planar Gaussian onto two unit directions
at angle \(\arccos r\); their signs disagree with probability
\(\arccos r/\pi\).

Oddness of arcsine and \((\arcsin)'\ge1\) on \((-1,1)\) now give
\[
\begin{aligned}
\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
&=\frac2\pi\sum_{i<j}
 [\arcsin(u_{ij}+v)-\arcsin(u_{ij}-v)]\\
&\ge\frac2\pi\binom n2\,2v
=\frac{2n\sqrt{n-1}}\pi.
\end{aligned}
\]
Both expectations belong to \([-\Phi(A),\Phi(A)]\). Dividing their
difference bound by two proves the claim. Minimizing over \(A\) yields
\[
m_n\ge\frac{n\sqrt{n-1}}\pi,\qquad
\liminf_n\alpha_n\ge\frac1\pi.
\]
In particular \(m_n=\Theta(n^{3/2})\), without assuming that any particular
family is optimal.

**Later strict asymptotic corollary (2026-09-06).** The reviewed
[second-moment-tail source theorem](evidence/NOTE_2026-09-06_ALL_LAW_SECOND_MOMENT_TAIL_GAIN.md),
Section 8, strengthens the asymptotic lower to
\[
\liminf_{n\to\infty}\alpha_n>\frac1\pi.
\]
Equivalently, some unspecified \(\varepsilon_0>0\) and \(n_0\) give
\(\Phi(A)/n^{3/2}\ge1/\pi+\varepsilon_0\) for every complete signing
of order \(n\ge n_0\). The proof excludes hypothetical saturation of
the exact nuclear lower via full-source spectral-tail control and gain.
It supplies neither a numerical \(\varepsilon_0\), convergence, nor a
limit value; its \(F_1(1)\) is not an unconditional lower bound.

**Explicit paired-field improvement (2026-09-17).** The author-reviewed
[paired polynomial-phase theorem](evidence/NOTE_2026-09-17_PAIRED_POLYNOMIAL_FIELD_LOWER.md)
now gives, with \(\kappa=2/\pi\),
\[
\liminf_n\alpha_n\ge
B:=\frac{81\kappa/2+9\sqrt{\kappa(2/3+\kappa/2)}}{101}
=0.3258407554\ldots>\frac{13}{40}.
\]
The paired local-field bound is independent of a fixed operator cap;
same-order regularization then removes that cap. This supplies neither
an explicit finite-order cutoff nor convergence. Its exact scalar
certificate and author-review scope are recorded with the proof.

## 5. Transfer along ratio-dense subsequences

Let \(2\le n_1<n_2<\cdots\) be an unbounded sequence of integers with
\(n_{k+1}/n_k\to1\). Then
\[
\liminf_n\alpha_n=\liminf_k\alpha_{n_k},\qquad
\limsup_n\alpha_n=\limsup_k\alpha_{n_k}.
\]
For \(n_k\le N\le n_{k+1}\), restriction gives
\[
\alpha_{n_k}\left(\frac{n_k}{N}\right)^{3/2}
\le\alpha_N\le
\alpha_{n_{k+1}}\left(\frac{n_{k+1}}N\right)^{3/2}.
\]
Both ratios tend to one uniformly in this interval, and \(\alpha\) is
bounded. Taking lower and upper limits proves the assertion. Thus
convergence along such a subsequence is equivalent to convergence of the
whole sequence. This says nothing about which signing attains \(m_{n_k}\).

## 6. Conference constructions give \(\limsup\alpha_n\le1/2\)

Here the number-theoretic input is the classical prime-number theorem in
the fixed progression \(1\pmod4\):
\[
\#\{q\le X:q\text{ prime},\ q\equiv1\pmod4\}
\sim\frac{X}{2\log X}.
\]
In particular, the successive such primes satisfy
\(q_k\sim2k\log k\), hence \((q_{k+1}+1)/(q_k+1)\to1\).

For completeness the needed matrix construction is explicit. For such a
prime \(q\), let \(\chi\) be the quadratic character on \(\mathbf F_q\),
with \(\chi(0)=0\), and index \(C\) by \(\{\infty\}\cup\mathbf F_q\):
\[
C_{\infty\infty}=0,\quad C_{\infty u}=C_{u\infty}=1,\quad
C_{uv}=\chi(u-v).
\]
It is symmetric because \(\chi(-1)=1\), and its off-diagonal entries are
signs. The character identities
\[
\sum_t\chi(t)=0,\qquad
\sum_t\chi(t)\chi(t-a)=-1\quad(a\ne0)
\]
give \(C^2=qI\). To verify the second identity, scale \(t=as\) and
count \(y^2=s(s-1)\): the factorization
\((2s-1-2y)(2s-1+2y)=1\) gives exactly \(q-1\) pairs \((s,y)\),
whereas their number is \(q+\sum_s\chi(s(s-1))\).

Writing \(r=q+1\), the spectral bound on the Boolean cube gives
\[
m_r\le\Phi(C)\le\tfrac12 r\|C\|_{\rm op}
=\tfrac12 r\sqrt{r-1}.
\]
Apply Section 5 to the increasing orders \(r_k=q_k+1\). It follows that
\[
\frac{13}{40}<B\le\liminf_n\alpha_n
\le\limsup_n\alpha_n\le\frac12.
\]
The lower bound is the September 17 result in Section 4; the conference
construction and ratio-dense transfer supply only the upper bound here.
This uses conference matrices only as admissible constructions for an
upper bound. It does not assert that they minimize \(\Phi\), or that the
two ends of the sandwich coincide.

## 7. Optional conditional two-ray criterion

This is a sufficient criterion, not a required architecture for solving
the original problem. No necessity assertion is made, and its
amplification hypotheses are not proved here.

Put \(H(n)=m_n^{2/3}\) and \(h(n)=H(n)/n=\alpha_n^{2/3}\). For a
nonnegative function \(\eta\) on the integers \(n\ge2\), define
\[
\eta^*(N)=\sup_{u\in\mathbf Z,\ u\ge N}\eta(u),\qquad
E(N)=\sum_{j\ge0}\eta^*(2^jN).
\]
Suppose \(E(N)\to0\) and, for all sufficiently large integers \(n\),
\[
H(2n)\le2H(n)+2n\eta(n),\qquad
H(3n)\le3H(n)+3n\eta(n).
\]
Then \(\alpha_n\) converges.

Proof: dividing by the new argument gives
\(h(qn)\le h(n)+\eta(n)\), for \(q=2,3\). Along a word in these
multipliers, the argument before step \(j\) is at least \(2^jn\).
Consequently, uniformly for integers \(a,b\ge0\),
\[
h(2^a3^bn)\le h(n)+E(n).
\]
The sorted semigroup \(\{2^a3^b:a,b\ge0\}\) has consecutive ratios
tending to one. Indeed, irrationality of \(\log2/\log3\) makes the
residues of \(a\log2\) dense modulo \(\log3\). For any \(\varepsilon>0\),
a finite prefix has all circular gaps below \(\varepsilon\). Adding
nonnegative multiples of \(\log3\) therefore puts a semigroup logarithm
in every sufficiently large interval of length \(\varepsilon\).

Fix an eligible \(n\). For large \(N\), choose \(s=2^a3^b\) with
\(N\le sn\) and \(sn/N\to1\). Monotonicity of \(H\) gives
\[
h(N)\le\frac{sn}{N}h(sn)
\le\frac{sn}{N}\bigl(h(n)+E(n)\bigr).
\]
Hence \(\limsup_Nh(N)\le h(n)+E(n)\). Let \(n\) tend to infinity
along a liminf subsequence; since \(E(n)\to0\), the lower and upper
limits agree. Finally \(\alpha_n=h(n)^{3/2}\) converges.

For example, \(\eta(n)=O(n^{-\delta})\), \(\delta>0\), or
\(\eta(n)=O((\log n)^{-1-\varepsilon})\), \(\varepsilon>0\), satisfies
the stated tail condition. Merely knowing \(\eta(n)\to0\) does not
establish that condition.

## 8. The Paley sequence \(n=p^2+1\) determines the cluster set

Let \(p_k\) be the odd primes and \(n_k=p_k^2+1\). The prime-number
theorem gives \(p_{k+1}/p_k\to1\), hence \(n_{k+1}/n_k\to1\). Section 5
therefore yields
\[
\liminf_n\alpha_n=\liminf_k\alpha_{n_k},\qquad
\limsup_n\alpha_n=\limsup_k\alpha_{n_k}.
\]
In particular, \(\alpha_n\) converges if and only if \(\alpha_{n_k}\)
does, and the possible limits coincide. Multipliers two and three are
not required for this equivalence.

For each odd prime \(p\), the Paley conference matrix \(C\) of order
\(n=p^2+1\) over \(\mathbf F_{p^2}\) admits an explicit Boolean
eigenvector with \(Cx=px\) (halfspace of an \(\mathbf F_p\)-linear form;
see the archived \(\rho=1\) note). Consequently
\[
\Phi(C)=\frac12 np=\frac12 n\sqrt{n-1},
\]
and \(m_n\le\Phi(C)\). Write
\[
\gamma_p=\frac{\Phi(C)-m_n}{n^{3/2}}
=\frac12\sqrt{1-\frac1n}-\alpha_n\ge0.
\]
Then \(\alpha_{n_k}=\frac12\sqrt{1-1/n_k}-\gamma_{p_k}\). As
\(\sqrt{1-1/n_k}\to1\),
\[
\liminf_n\alpha_n=\frac12-\limsup_k\gamma_{p_k},\qquad
\limsup_n\alpha_n=\frac12-\liminf_k\gamma_{p_k}.
\]
Existence of \(\lim\alpha_n\) is therefore equivalent to existence of
\(\lim\gamma_{p_k}\). The value, if the limit exists, is
\(\frac12-\lim\gamma_p\).

Two exclusive alternatives, still open:

- \(\gamma_{p_k}\to0\). Then \(\alpha_{n_k}\to\frac12\), hence
  \(\lim\alpha_n=\frac12\).
- \(\limsup_k\gamma_{p_k}=c>0\). Then
  \(\liminf_n\alpha_n\le\frac12-c\). If in fact \(\gamma_{p_k}\to c\),
  the limit exists and equals \(\frac12-c\). If \(\gamma_{p_k}\) has two
  distinct cluster points, \(\alpha_n\) diverges.

A matching-size or \(O(n)\) undercut of \(C\) forces \(\gamma_p=O(n^{-1/2})\to0\)
and would settle the first alternative. That implication is not available
as a proof: the only certified Paley-order undercut at \(n=26\) has
switching distance \(122=\Theta(n^{3/2})\) from \(C\) (see
`evidence/NOTE_2026-09-19_BEATER_SHELL.md`).

The CORE Gaussian bound is saturated at conference: \(A^2=(n-1)I\)
forces the pair quantities \(u_{ij}=0\), so the arcsine difference is
exactly \(2\arcsin(1/\sqrt{n-1})\) and the extra factor
\(\arcsin(v)/v\to1\). That estimate cannot be improved to \(\tfrac12\)
at Paley orders by the same method. Lipschitz plus small switching
distance is likewise unavailable: the certified \(n=26\) undercutter
has switching Hamming distance \(122=\Theta(n^{3/2})\) from Paley, with
\(A^2\) off-diagonal in \(\{0,\pm4\}\) and six nonzeros per row. It is
the plus-I lift \(K(B)=\bigl(\begin{smallmatrix}B&B+I\\B+I&-B\end{smallmatrix}\bigr)\)
of an exact \(m_{13}=20\) block (see
`evidence/NOTE_2026-09-19_PLUS_I_COHERENT_LIFT.md`). The same lift of the
5-cycle is the exact minimizer at \(n=10\).
No reviewed construction gives a uniform \(c>0\), and no reviewed
argument forces \(\gamma_p\to0\). The uniform lower bound of Section 4
only yields \(\gamma_p\le\tfrac12-B+o(1)\). This section does not reopen
Paley gap-2 covers or Hadamard families.
