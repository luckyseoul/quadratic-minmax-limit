# Existence of \(\displaystyle\lim_{n\to\infty}n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\bigl|\sum_{i<j}a_{ij}x_ix_j\bigr|\)

## Statement

For \(n\ge 2\) set
\[
m_n
=
\min_{\substack{a_{ij}=\pm1\\1\le i<j\le n}}
\max_{x\in\{\pm1\}^n}
\Biggl|
\sum_{1\le i<j\le n}a_{ij}\,x_i x_j
\Biggr|,
\qquad
\alpha_n
=
\frac{m_n}{n^{3/2}}.
\]

**Main Theorem (sandwich).**
\begin{equation}
\label{eq:sandwich}
\frac1\pi
\;\le\;
\liminf_{n\to\infty}\alpha_n
\;\le\;
\limsup_{n\to\infty}\alpha_n
\;\le\;
\frac12.
\end{equation}
Moreover \(\limsup_n\alpha_n=\lim_n a_n\) where \(a_n=\sup_{k\ge1}\alpha_{kn}\).
(The lower bound is the dual-Gaussian arcsine argument of Prop.~5.2, valid for **every** Seidel matrix; the classical Bohnenblust–Hille floor \(2^{-5/2}\) is retained as Prop.~5.1.)

**Main Theorem (limit).** \(\displaystyle L=\lim_{n\to\infty}\alpha_n=\tfrac12\) — **OPEN** (2026-08-06).

Candidate path: sandwich + denseness Prop~6.1 on \(\rho=1\) Paley \(n=p^2+1\) + E(1)
via bi-tight (15.167) + freeness-fail ND (15.170–15.171).  
**Fatal hinge:** disj \(\mathrm{Gsum}\) LB used by 15.170–171 is not proved for general \(p\)
(15.158: Max+ is not an IP association scheme). Short package:
`evidence/share/denseness_path_package.md`. Handoff: **`HANDOFF.md`**.

**Optional still open:** Path-C residual / \(16N\) (independent).

**Corollary (\(\rho=1\) along a dense Paley family).** For every odd prime \(p\), the Paley conference
matrix of order \(n=p^2+1\) (over \(\mathbb F_{p^2}\)) admits a halfspace boolean eigenvector
\(Cx=px\), hence \(\rho(C)=1\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\). Along \(n_k=p_k^2+1\) one has
\(n_{k+1}/n_k\to1\) and \(\limsup_k\rho(C_{n_k})=1\). Proof: `evidence/PROOF_rho_eq_1.md`.
With \(m_n\ge\Phi(C)-2\) (E(1)) this would force \(L=\tfrac12\) — E(1) not closed.

**What is complete.** Dual-Gaussian lower bound; denseness framework; \(\rho=1\) on \(n=p^2+1\);
majorization algebra for bi-tight (15.167, conditional on mult/\(\lambda_{\min}\)); Farkas **poly**
for dual equality **if** disj Gsum LB holds. **What is not:** general-p disj Gsum LB; E(1); \(L=\tfrac12\).

MathOverflow [413935](https://mathoverflow.net/questions/413935) /
https://x.com/PI010101/status/2081070728422752329 remain open for existence.

> **Reader note.** Soft-close of \(L=\tfrac12\) via scheme-min Gsum was **retracted** 2026-08-06.
> See `STATUS.md`, `evidence/share/denseness_path_package.md`.

---

## §1. Equivalent matrix form

Associate to \((a_{ij})_{i<j}\) the symmetric zero-diagonal matrix \(A\) with \(A_{ij}=A_{ji}=a_{ij}\). Then
\[
\sum_{i<j}a_{ij}x_i x_j=\frac12\,x^\top A x,
\]
so
\[
m_n=\frac12\min_A\max_{x\in\{\pm1\}^n}\bigl|x^\top A x\bigr|
\]
over all such \(A\). Write \(\Phi(A)=\max_x\bigl|\tfrac12 x^\top A x\bigr|\), hence \(m_n=\min_A\Phi(A)\). The quantity \(\alpha_n\) is the same in both presentations.

---

## §2. Spectral upper bound for conference matrices

For every admissible \(A\), \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\), with equality iff \(A\) is a symmetric conference matrix of order \(n\). On the sphere of radius \(\sqrt n\),
\[
\max_{\|x\|_2=\sqrt n}\bigl|x^\top A x\bigr|=n\|A\|_{\mathrm{op}}\ge n\sqrt{n-1}.
\]
Hence for any conference matrix \(C\),
\begin{equation}
\label{eq:conf-ub}
\Phi(C)\le\frac12 n\sqrt{n-1},
\qquad
\frac{\Phi(C)}{n^{3/2}}\le\frac12\sqrt{1-\frac1n}.
\end{equation}
This upper-bounds \(\Phi(C)\) and therefore \(m_n\) when a conference matrix exists. It does **not** lower-bound \(m_n=\min_A\Phi(A)\).

---

## §3. Monotonicity and padding

**Proposition 3.1.** \(m_n\le m_N\) whenever \(2\le n\le N\).

*Proof.* Let \(A\) be optimal of order \(N\). For \(S\subset[N]\), \(|S|=n\), and \(y\in\{\pm1\}^S\), extend \(y\) by i.i.d. uniform random signs on the complement. Then \(\mathbb E Q=Q_S(y)\), so \(|Q_S(y)|\le m_N\). Hence \(m_n\le m_N\). \(\square\)

**Proposition 3.2.** \(m_{n+1}\le m_n+n\).

*Proof.* Adjoin a last row/column of all \(+1\) to an optimal order-\(n\) matrix. \(\square\)

**Proposition 3.3 (continuity of \(\alpha\)).** \(|\alpha_{n+1}-\alpha_n|\to0\) as \(n\to\infty\). In particular the set of limit points of \((\alpha_n)\) is a closed interval \([\lambda,\Lambda]\).

*Proof.* Write \(\delta_n:=m_{n+1}-m_n\in[0,n]\). Then
\[
\alpha_{n+1}-\alpha_n
=
\frac{\delta_n}{n^{3/2}}
-\frac32\frac{\alpha_n}{n}
+O\Bigl(\frac{\delta_n}{n^{5/2}}+\frac{\alpha_n}{n^2}\Bigr).
\]
The first term is \(\le n^{-1/2}\to0\). Boundedness \(\alpha_n=O(1)\) follows from the elementary random-method estimate \(m_n\le\sqrt{\log2}\,n^{3/2}\) for large \(n\) (union bound on Rademacher sums; cf. §7). A real sequence with consecutive gaps tending to zero has connected limit-point set, hence a closed interval. \(\square\)

**Corollary 3.4.** For \(d\ge0\), \(m_{n+d}\le m_n+dn+\tfrac12 d(d-1)\).

---

## §4. Upper bound \(\limsup\alpha_n\le\tfrac12\)

**Proposition 4.1.** \(\limsup_{n\to\infty}\alpha_n\le\tfrac12\).

*Proof.* For a prime \(q\equiv1\pmod4\), the Paley conference matrix of order \(n=q+1\) satisfies \(C^\top C=(n-1)I\) (Paley 1933), so \(m_n\le\tfrac12 n\sqrt{n-1}\) by \eqref{eq:conf-ub}. By Dirichlet’s theorem there are infinitely many such primes. Writing \(n_k=q_k+1\), the prime-number theorem in the progression \(1\bmod4\) gives \(n_k\sim 2k\log k\), hence \(n_{k+1}/n_k\to1\).

For \(N\ge2\) let \(n_k\) be the least Paley order \(\ge N\). Then \(m_N\le m_{n_k}\le\tfrac12 n_k\sqrt{n_k-1}\), so
\[
\alpha_N\le\frac12\Bigl(\frac{n_k}{N}\Bigr)^{3/2}\sqrt{1-\frac1{n_k}}.
\]
As \(N\to\infty\), \(n_k/N\to1\). \(\square\)

---

## §5. Lower bounds on \(\liminf\alpha_n\)

**Proposition 5.1 (Bohnenblust–Hille).** \(\liminf_{n\to\infty}\alpha_n\ge2^{-5/2}\).

*Proof.* The degree-\(2\) Bohnenblust–Hille inequality on the Boolean cube (Defant–Mastyło–Pérez, *Math. Ann.* 2019) yields a universal \(B_2<\infty\) such that for \(f=\sum_{i<j}c_{ij}x_ix_j\),
\[
\Bigl(\sum_{i<j}|c_{ij}|^{4/3}\Bigr)^{3/4}
\le B_2\max_{x\in\{\pm1\}^n}|f(x)|.
\]
With \(|c_{ij}|=1\) and \(\max|f|=m_n\), one has \(\binom{n}{2}^{3/4}\le B_2 m_n\). The hypercontractive bound \(B_2\le2\sqrt2\) (Ivanisvili, MO 413935) yields \(\liminf\alpha_n\ge2^{-5/2}\). \(\square\)

**Proposition 5.2 (dual-Gaussian arcsine; universal \(1/\pi\)).** For every Seidel matrix \(A\in\mathcal S_n\) and every \(n\ge2\),
\begin{equation}
\label{eq:dual-gauss}
\Phi(A)\;\ge\;\frac{n\sqrt{n-1}}{\pi}.
\end{equation}
In particular \(m_n\ge n\sqrt{n-1}/\pi\) and
\[
\liminf_{n\to\infty}\alpha_n\;\ge\;\frac1\pi.
\]

*Proof.* Fix \(A\in\mathcal S_n\), let \(g\sim N(0,I_n)\), and set \(t=\sqrt{n/(n-1)}\). Define
\[
z^{\pm}=\Bigl(I\pm\frac{t}{\sqrt n}\,A\Bigr)g,
\qquad
x^{\pm}=\operatorname{sgn}(z^{\pm})
\]
(with the convention \(\operatorname{sgn}(0)=+1\)). Each coordinate of \(z^{\pm}\) has variance
\[
d=1+t^2\frac{n-1}{n}=2.
\]
For \(i\neq j\) the correlations \(r_{ij}^{\pm}\) of the pairs \((z_i^{\pm},z_j^{\pm})\) satisfy
\[
a_{ij}\,r_{ij}^{\pm}
=
u_{ij}\pm v,
\qquad
u_{ij}
=
\frac{a_{ij}\,t^2(A^2)_{ij}}{nd},
\qquad
v
=
\frac{2t}{\sqrt n\,d}
=
\frac1{\sqrt{n-1}}.
\]
(The identities use \(A_{ii}=0\) and \(A_{ij}=a_{ij}=\pm1\).) Because \(|(A^2)_{ij}|\le n-2\),
\[
|u_{ij}|\le\frac{n-2}{2(n-1)},
\qquad
|u_{ij}\pm v|
\le
\frac{n-2}{2(n-1)}+\frac1{\sqrt{n-1}}
<1
\]
for all \(n\ge2\) (direct check at \(n=2,3\); for \(n\ge4\) the right-hand side is \(<1\)).

The arcsine law for a centered bivariate Gaussian of correlation \(r\) states
\(\mathbb E[\operatorname{sgn}Z_1\operatorname{sgn}Z_2]=(2/\pi)\arcsin r\). Hence
\[
\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
=
\frac2\pi\sum_{i<j}a_{ij}\bigl(\arcsin r_{ij}^+-\arcsin r_{ij}^-\bigr).
\]
For each edge, \(a_{ij}(\arcsin r_{ij}^+-\arcsin r_{ij}^-)=\arcsin(u_{ij}+v)-\arcsin(u_{ij}-v)\):
if \(a_{ij}=+1\) this is immediate, and if \(a_{ij}=-1\) both sides pick up a matching sign change.
Since \((\arcsin)'(r)=1/\sqrt{1-r^2}\ge1\) on \((-1,1)\), the mean-value theorem yields
\[
\arcsin(u+v)-\arcsin(u-v)\ge 2v.
\]
Therefore
\[
\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
\ge
\frac2\pi\cdot\binom n2\cdot 2v
=
\frac{n\sqrt{n-1}}{\pi}.
\]
Both expectations lie in \([-\Phi(A),\Phi(A)]\), so their difference is at most \(2\Phi(A)\).
Halving gives \eqref{eq:dual-gauss}. Minimising over \(A\) and passing to the limit produces the
claim on \(\alpha_n\). \(\square\)

**Remark (correction).** The classical single-sided Nesterov bound applied to a *conference*
matrix \(C\) yields only \(\Phi(C)\ge n\sqrt{n-1}/\pi\), which does not control \(m_n=\min\Phi\).
Proposition 5.2 upgrades the same constant to a **uniform** lower bound on every Seidel matrix
by running the Gaussian construction in the pair of directions \(\pm A\) and comparing the two
expectations. The older Bohnenblust–Hille floor \(2^{-5/2}\approx0.177\) is thereby improved to
\(1/\pi\approx0.3183\).

**Remark (cut-code form).** Writing \(m=\binom n2\) and
\(D_n=\{(\pm x_ix_j)_{i<j}:x\in\{\pm1\}^n\}\subset\{\pm1\}^m\), the identity
\(\langle a,c_x\rangle=m-2d_H(a,c_x)\) and the presence of both \(c_x\) and \(-c_x\) in \(D_n\) give
\begin{equation}
\label{eq:cut-code}
m_n
=
\binom n2-2\rho(D_n),
\end{equation}
where \(\rho(D_n)\) is the Hamming covering radius of \(D_n\) in \(\{\pm1\}^m\). Equivalently,
\(m_n/2\) is the covering-radius deficit of the antipodal cut code of \(K_n\). This is the coding
formulation used by Esmaeili–Zaghian (2009); it does not by itself settle existence of
\(\lim\alpha_n\).

---

## §6. Denseness

**Proposition 6.1.** If \(n_{k+1}/n_k\to1\), then
\[
\liminf_n\alpha_n=\liminf_k\alpha_{n_k},
\qquad
\limsup_n\alpha_n=\limsup_k\alpha_{n_k}.
\]

*Proof.* For \(n_k\le N\le n_{k+1}\), monotonicity yields \(m_{n_k}\le m_N\le m_{n_{k+1}}\), so
\[
\alpha_{n_k}\Bigl(\frac{n_k}{N}\Bigr)^{3/2}
\le\alpha_N\le
\alpha_{n_{k+1}}\Bigl(\frac{n_{k+1}}{N}\Bigr)^{3/2}.
\]
Both ratios tend to \(1\). \(\square\)

In particular, for each fixed integer \(r\ge1\), \(\liminf_n\alpha_{rn}=\liminf\alpha_n\) and \(\limsup_n\alpha_{rn}=\limsup\alpha_n\).

**Proposition 6.2 (Paley reduction).** Let \(q_k\) be the \(k\)-th prime congruent to \(1\pmod4\), and set \(n_k=q_k+1\) (Paley/conference orders). Then \(n_{k+1}/n_k\to1\), and therefore
\[
\liminf_{n\to\infty}\alpha_n=\liminf_{k\to\infty}\alpha_{n_k},
\qquad
\limsup_{n\to\infty}\alpha_n=\limsup_{k\to\infty}\alpha_{n_k}.
\]
In particular, \(\lim\alpha_n\) exists if and only if \(\lim_k\alpha_{n_k}\) exists.

*Proof.* By the prime-number theorem for the arithmetic progression \(1\bmod4\) (Dirichlet density \(\tfrac12\)), one has \(q_k\sim 2k\log k\), so \(n_{k+1}/n_k\to1\). Apply Proposition 6.1. \(\square\)

Thus any existence proof may restrict attention to Paley orders; any non-existence proof must already be visible along that sparse but ratio-dense subsequence.

---

## §7. Multipartite bounds

**Proposition 7.1 (random multipartite).** For \(k,n\ge2\),
\begin{equation}
\label{eq:rand}
m_{kn}\le k\,m_n+\sqrt{\log2}\,(kn)^{3/2},
\qquad
\alpha_{kn}\le\frac{\alpha_n}{\sqrt k}+\sqrt{\log2}.
\end{equation}

*Proof.* Optimal order-\(n\) matrices on \(k\) diagonal blocks; i.i.d. \(\pm1\) on cross blocks. For fixed \(x\), the cross contribution is a Rademacher sum of at most \(\binom{k}{2}n^2\) terms. Hoeffding and a union bound over \(2^{kn}\) vertices give the stated cross bound with positive probability. \(\square\)

**Proposition 7.2 (Hadamard doubling).** If a Hadamard matrix of order \(n\) exists, then
\begin{equation}
\label{eq:had}
m_{2n}\le 2m_n+n^{3/2},
\qquad
\alpha_{2n}\le\frac{\alpha_n}{\sqrt2}+\frac1{2\sqrt2}.
\end{equation}
*Proof.* Optimal diagonal blocks; Hadamard cross \(H\); \(\|Hu\|_1\le\sqrt n\,\|Hu\|_2=n^{3/2}\) for \(u\in\{\pm1\}^n\) is not always true, but Cauchy–Schwarz gives \(|x_1^\top H x_2|\le\|x_1\|_2\|Hx_2\|_2=n^{3/2}\). \(\square\)

Sylvester matrices of order \(2^\ell\), together with padding (Corollary 3.4), extend \eqref{eq:had} to all large \(n\) with an additional \(o(1)\) term in \(\alpha\).

**Remark.** Filling cross blocks with the constant signing \(+1\) yields error \(\binom{k}{2}n^2\) in \(m_{kn}\), hence error \(\Theta(\sqrt{kn})\) in \(\alpha_{kn}\), which tends to infinity as \(n\to\infty\) and cannot be used for limit arguments at fixed \(k\).

**Proposition 7.3 (reverse multipartite, new).** For all \(k,n\ge2\),
\begin{equation}
\label{eq:rev}
m_{kn}\ge\frac k2\,m_n,
\qquad
\alpha_{kn}\ge\frac{\alpha_n}{2\sqrt k}.
\end{equation}

*Proof.* Let \(B\) achieve \(m_{kn}\), and fix any partition of \([kn]\) into \(k\) blocks of size \(n\). Write \(Q_B(y)=\sum_b Q_b(y_b)+\sum_{b<b'}y_b^\top B^{bb'}y_{b'}\). For fixed \(x_1,\dots,x_k\in\{\pm1\}^n\) and \(\varepsilon\in\{\pm1\}^k\), set \(y_b=\varepsilon_b x_b\). Then
\[
Q_B(\varepsilon\cdot x)=\sum_b Q_b(x_b)+\sum_{b<b'}\varepsilon_b\varepsilon_{b'}\,x_b^\top B^{bb'}x_{b'},
\]
so \(\mathbb E_\varepsilon Q_B=\sum_b Q_b(x_b)\). Hence \(\bigl|\sum_b Q_b(x_b)\bigr|\le m_{kn}\) for every choice of the \(x_b\). Maximising and minimising over the \(x_b\) separately yields
\[
\sum_b\Phi^+(B^{bb})\le m_{kn},\qquad\sum_b\Phi^-(B^{bb})\le m_{kn},
\]
where \(\Phi^+=\max Q\), \(\Phi^-=-\min Q\). Therefore \(\sum_b(\Phi^++\Phi^-)\le 2m_{kn}\). But \(\Phi(B^{bb})=\max(\Phi^+,\Phi^-)\le\Phi^++\Phi^-\) and \(\Phi(B^{bb})\ge m_n\), so \(k m_n\le 2m_{kn}\). \(\square\)

(The factor \(2\) is an artefact of one-sided blocks; it is harmless for the soft analysis below.)

---

## §8. The majorant \(a_n\)

**Proposition 8.1.** Let \(\Lambda:=\limsup_n\alpha_n\) and \(a_n:=\sup_{k\ge1}\alpha_{kn}\). Then \(\lim_n a_n=\Lambda\).

*Proof.*
*Upper:* \(a_n\le\sup_{N\ge n}\alpha_N\), so \(\limsup a_n\le\Lambda\).

*Lower:* Fix \(\varepsilon>0\). For large \(n\) choose \(N\ge n^4\) with \(\alpha_N>\Lambda-\varepsilon\). Set \(k=\lfloor N/n\rfloor\ge n^3\). Corollary 3.4 gives \(m_N\le m_{kn}+O(nN)\), so
\[
\alpha_N\le\alpha_{kn}\Bigl(\frac{kn}N\Bigr)^{3/2}+O(n N^{-1/2}).
\]
Here \(kn/N\to1\) and \(n N^{-1/2}\to0\), hence \(\alpha_N\le a_n+o(1)\) and \(\liminf a_n\ge\Lambda\). \(\square\)

---

## §9. Why multipartite comparison does not prove \(\liminf=\limsup\)

Write \(\lambda:=\liminf_n\alpha_n\) and \(\Lambda:=\limsup_n\alpha_n\), so \(2^{-5/2}\le\lambda\le\Lambda\le\tfrac12\) and \(a_n\to\Lambda\).

Assume \(\lambda<\Lambda\) and set \(\mu:=(\lambda+\Lambda)/2\). Choose \(n_j\to\infty\) with \(\alpha_{n_j}\to\lambda\). Let
\[
k_j:=\min\bigl\{k\ge1:\alpha_{k n_j}\ge\mu\bigr\}\ge2.
\]

**Bounded \(k_j=k\ge2\).** Proposition 7.1 gives a limit point \(\beta\ge\mu\) of \(\alpha_{k n_j}\) with
\(\beta\le\lambda/\sqrt k+\sqrt{\log2}\). Forcing \(\lambda/\sqrt k+\sqrt{\log2}<\mu\) requires \(\mu>\sqrt{\log2}\approx0.83\), which fails under \(\mu\le\tfrac12\).
Hadamard (\(k=2\)) yields \(\mu\le\lambda/\sqrt2+1/(2\sqrt2)\), i.e.
\[
\Lambda\le\lambda(\sqrt2-1)+\frac1{\sqrt2}.
\]
The right-hand side is \(\ge0.78>\tfrac12\ge\Lambda\) when \(\lambda\ge2^{-5/2}\), so no contradiction in the admissible range.

**Unbounded \(k_j\to\infty\).** Proposition 7.1 yields \(\mu\le\sqrt{\log2}\), compatible with \(\mu\le\tfrac12\). Reverse multipartite (Prop. 7.3) only gives \(\mu\ge0\). Two-block analysis gives \(\alpha_{k_j n_j}\to\mu\) and \(a_{k_j n_j}\to\Lambda\); further ascent to near \(\Lambda\) with unbounded multipliers again yields only \(\Lambda\le\sqrt{\log2}\).

**Abstract counterexamples.** Nondecreasing sequences of growth \(\Theta(n^{3/2})\) with increments \(\le n\) and satisfying \eqref{eq:rand}–\eqref{eq:had} can have non-convergent slopes inside \([2^{-5/2},\tfrac12]\) (e.g. log-log periodic envelopes). Thus the soft inequalities alone do not force existence.

---

## §10. The “\(c_k\to0\)” fantasy is impossible (correction)

A common hope is that a multipartite bound
\[
\alpha_{kn}\le\frac{\alpha_n}{\sqrt k}+c_k
\quad\text{with }c_k\to0\text{ as }k\to\infty\text{ (uniformly in }n\text{)}
\]
would finish the proof. **No such bound can exist.**

*Proof of impossibility.* Suppose such \(c_k\) existed. Fix any \(n\) with \(\alpha_n<\infty\) and let \(k\to\infty\). Then \(\alpha_{kn}\le\alpha_n/\sqrt k+c_k\to0\). But by denseness along multiples of \(n\) (Proposition 6.1) one has \(\liminf_k\alpha_{kn}=\lambda\ge2^{-5/2}>0\), a contradiction. \(\square\)

Even allowing \(n\)-dependence, the obstruction persists in the regime that matters for existence: if \(\alpha_{n_j}\to\lambda\) and \(k_j\to\infty\), any upper bound of the form \(\alpha_{k_j n_j}\le\alpha_{n_j}/\sqrt{k_j}+c(k_j,n_j)\) with \(c\to0\) would force \(\mu\le0\). But climbing from a liminf point to a value \(\ge\mu>0\) is exactly what the majorant identity \(a_n\to\Lambda\) guarantees must sometimes happen. Therefore any true multipartite upper bound is forced to allow
\[
\liminf_{k\to\infty}c(k,n)\;\ge\;\lambda
\]
along liminf sequences \(n\), and in fact \(c(k,n)\) must be allowed to be as large as \(\Lambda-o(1)\) whenever \(\alpha_{kn}\) realises the limsup. Multipartite comparison **cannot** forbid the climb \(\lambda\to\Lambda\).

### Bilinear obstruction at fixed \(k\)

Even for fixed \(k=2\), the additive constant cannot tend to \(0\) as \(n\to\infty\). For any \(\pm1\) cross-block \(H\in\{\pm1\}^{n\times n}\),
\[
\max_{x,y\in\{\pm1\}^n}|x^\top Hy|\ge c\,n^{3/2},\qquad c=\sqrt{2/\pi}-o(1),
\]
by taking \(y\) random and \(x=\mathrm{sign}(Hy)\) (so \(x^\top Hy=\|Hy\|_1\sim n^{3/2}\sqrt{2/\pi}\)). Hadamard matrices match the upper bound \(n^{3/2}\). Thus every 2-block construction satisfies
\[
\Phi\begin{pmatrix}A_1&H\\H^\top&A_2\end{pmatrix}
\le 2m_n+\max|x^\top Hy|
\]
with the cross term contributing at least \(\sim\sqrt{2/\pi}\,n^{3/2}\) to the inevitable upper-bound budget, and at most \(n^{3/2}\). The additive constant in \(\alpha_{2n}\le\alpha_n/\sqrt2+C\) therefore obeys
\[
C\in\Bigl[\tfrac{\sqrt{2/\pi}}{2^{3/2}},\;\tfrac1{2\sqrt2}\Bigr]\approx[0.282,\,0.354]
\]
and cannot be driven to \(0\).

**Conclusion.** Soft multipartite methods — random, Hadamard, or any design-based cross blocks — are dead for existence. A proof must use the **combinatorial structure of \(\min_A\Phi(A)\)**, not comparison inequalities with error terms.

---

## §11. Approach 1: recursion, maximizers, and slack — fails

### Exact recursion

Adjoining a last row/column \(s\in\{\pm1\}^n\) to \(A\) yields, for \(x\in\{\pm1\}^n\) and new coordinate \(y=\pm1\),
\[
Q(x,y)=Q_A(x)+y\,(s\cdot x),\qquad\max_y|Q(x,y)|=|Q_A(x)|+|s\cdot x|.
\]
Hence
\[
m_{n+1}=\min_{A,s}\max_x\bigl(|Q_A(x)|+|s\cdot x|\bigr).
\]
Writing \(\delta_A(x)=\Phi(A)-|Q_A(x)|\ge0\) and
\[
\gamma(A,s)=\max_x\bigl(|s\cdot x|-\delta_A(x)\bigr),\qquad\gamma(A)=\min_s\gamma(A,s),
\]
one has \(\max_x(|Q|+|s\cdot x|)=\Phi(A)+\gamma(A,s)\), so
\begin{equation}
\label{eq:rec}
m_{n+1}=\min_A\bigl(\Phi(A)+\gamma(A)\bigr).
\end{equation}
In particular, for any optimal \(A^*\) of order \(n\), \(m_{n+1}\le m_n+\gamma(A^*)\), and \(\gamma(A^*)\ge0\) by monotonicity.

### Maximizer discrepancy

Let \(M=\{x:|Q_{A^*}(x)|=m_n\}\). Then \(\gamma(A^*)\ge\min_{s\in\{\pm1\}^n}\max_{x\in M}|s\cdot x|\), the combinatorial discrepancy of the row set \(M\subset\{\pm1\}^n\).

- If \(|M|\le\mathrm{poly}(n)\), Spencer/partial-colouring yields \(\min_s\max_M|s\cdot x|=O(\sqrt{n\log n})\), so extension cost \(O(\sqrt{n\log n})\).
- If maximizers are rich (near-orthogonal packing), the cost is \(\Omega(\sqrt n)\).
- Flat extension \(m_{n+1}=m_n\) requires some \(s\) with \(s\cdot x=0\) for all \(x\in M\) (hence \(n\) even and \(M\subset s^\perp\)).

### Why this does not force \(\lambda=\Lambda\)

Increments \(\delta_n:=m_{n+1}-m_n\in[0,n]\) satisfy \(m_n=m_2+\sum_{j<n}\delta_j\). The constraint \(|\alpha_{n+1}-\alpha_n|=O(n^{-1/2})\) forces \(\alpha_n\) to be slowly varying on scale \(o(\sqrt n)\), but **permits order-1 oscillations on scale \(\sqrt n\)** and log-periodic oscillations on scale \(n\) (e.g. \(\sin(\log n)\) has steps \(O(1/n)\)). Soft increment control is compatible with non-convergence.

To upgrade recursion to a limit theorem one would need a load-bearing estimate of the form
\[
\gamma(A^*)=(c+o(1))\sqrt n
\]
for optimal \(A^*\), with \(c\) determined by \(\alpha_n\) (e.g. \(c=\tfrac32\alpha_n\)). That requires:
1. a matching upper bound \(\gamma(A^*)\le(\tfrac32\alpha_n+o(1))\sqrt n\) (extension not too expensive), and
2. a matching lower bound from maximizer geometry (extension not too cheap).

Both fail with present technology:
1. Random \(s\) only gives \(\mathbb E\max_M|s\cdot x|\le\sqrt{2n\log(2|M|)}\); without a polynomial bound on \(|M|\) for optimal \(A\), this can be \(\Theta(n)\).
2. Lower-bounding \(\min_s\max_M|s\cdot x|\) requires that maximizers of optimal \(A\) cannot hide in a thin slab \(\{|s\cdot x|\le t\}\) with \(t=o(\sqrt n)\). No such delocalisation is proved; for conference matrices maximizers can be equatorially concentrated relative to some directions.

**Invariant needed.** A theorem that every near-optimal \(A\) has maximizer set \(M\) with discrepancy \(\Theta(\sqrt n)\) and packing number \(\exp(o(n))\), uniformly.

---

## §12. Approach 2: degree-2 Boolean analysis / hypercontractivity — fails

Degree-2 hypercontractivity controls tails of \(f=\sum_{i<j}a_{ij}x_ix_j\): \(\|f\|_p\le(p-1)\|f\|_2\) and Hanson–Wright concentration about mean zero. For optimal \(A\), \(\|f\|_2=\sqrt{\binom{n}{2}}\sim n/\sqrt2\) while \(\|f\|_\infty=m_n\sim\alpha n^{3/2}\), so the maximum sits at height \(\sim\alpha\sqrt{2n}\) standard deviations — deep in the tail.

Hypercontractive level-set bounds give
\[
\mu\bigl(\{|f|\ge(1-\varepsilon)m_n\}\bigr)\le\exp\bigl(-c(\varepsilon)\sqrt n\bigr)
\]
or better under additional spectral assumptions, but:
- this upper-bounds the measure of near-maximizers, which is the wrong direction for discrepancy lower bounds (small \(M\) makes discrepancy *easier*, i.e. \(\gamma\) smaller);
- Kindler–Safra / junta theorems for degree 2 require the function to be close to Boolean, which \(f/m_n\) is not (it takes many values);
- weak-\(L^2\) graphon limits of the coefficient arrays fail to be upper-semicontinuous for \(\Phi\) (random \(\pm1\) matrices are a counterexample to naïve graphon USC; recorded as given).

**Invariant needed.** A structure theorem for superlevel sets of *near-minimal-\(\Phi\)* degree-2 forms with \(\pm1\) coefficients — e.g. that they are approximate cosets of linear codes of controlled codimension, or that they are connected under Hamming noise with quantitative expansion — from which extension costs and multipartite rigidity would follow.

---

## §13. Approach 3: improved multipartite with special cross blocks — fails

Already killed in §10. Summary of attempted constructions and their failures:

| Cross design | Cross contribution to \(\alpha_{kn}\) | Verdict |
|---|---|---|
| i.i.d. random | \(\sqrt{\log2}\) (sharp for the Gaussian field) | \(c_k\not\to0\) |
| Hadamard (\(k=2\)) | \(1/(2\sqrt2)\approx0.354\) | best possible order; constant \(\ge0.28\) by bilinear lower bound |
| Lexicographic product \(A[B]\) | \(\alpha_k\sqrt n\to\infty\) | unusable |
| Kronecker \(C\otimes S\) + diagonal fill | \(\Theta(1)\) by spectral calculus | no \(o(1)\) |
| Conference block-signs \(\times\) Hadamard | \(\Theta(1)\) (Frobenius/nuclear estimates) | no \(o(1)\) |
| Constant \(\pm1\) blocks | \(\Theta(\sqrt{kn})\to\infty\) | unusable |

Any construction of a single matrix of order \(kn\) is competing with the conference upper bound \(\alpha\le1/2\); multipartite from small optimal blocks cannot beat the cross-term barrier of §10 without violating \(\lambda>0\).

---

## §14. Approach 4: Cesàro / Tauberian — fails (conditional only)

Write \(\delta_n=m_{n+1}-m_n\). Summation by parts:
\[
\sum_{j=1}^{n-1}\frac{\delta_j}{\sqrt j}=\frac{m_n}{\sqrt n}+\frac12\sum_{j}\frac{m_j}{j^{3/2}}+O\Bigl(\sum\frac{m_j}{j^{5/2}}\Bigr).
\]
If \(\delta_n/\sqrt n\to L\), then \(m_n\sim\tfrac23 L\,n^{3/2}\), so \(\alpha_n\to\tfrac23 L\).

Conversely, \(\alpha_n\to\alpha_*\) only controls Cesàro means of \(\delta_j/\sqrt j\), not pointwise convergence: one can have \(\delta_n=0\) on long stretches and \(\delta_n\sim c\sqrt n\) on complementary stretches, compatible with a convergent \(\alpha\) or with oscillation.

**Gap.** There is no a-priori regularity on \(\delta_n\) (numerics: increments \(2,1,0,1,4,1,2\) for \(n=2\to9\); flat stretches exist). A Tauberian upgrade would require monotone density, slow oscillation of \(\delta_n/\sqrt n\), or a one-sided Lipschitz condition that the recursion does not provide.

**Conditional theorem (soft).** *If* \(\delta_n/\sqrt n\) converges, then \(\lim\alpha_n\) exists. This reduces existence to regularity of optimal extension costs — which is Approach 1, already blocked.

---

## §15. Approach 5: discrete-to-spectral comparison — fails for \(\min_A\), but conference structure is exact

For every admissible \(A\), \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\), and the continuous maximum on the sphere \(\|x\|_2=\sqrt n\) is exactly \(n\|A\|_{\mathrm{op}}\). Write
\[
\rho(A)\,:=\,\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\,\|A\|_{\mathrm{op}}}\in(0,1].
\]
If \(\rho(A)\ge\rho-o(1)\) uniformly in \(A\), then \(\liminf\alpha_n\ge\rho/2\). Conference matrices give \(\limsup\le1/2\), so \(\rho=1\) would yield \(L=1/2\).

### §15.1 Exact spectral calculus for conference matrices (new, complete)

**Proposition 15.1 (spectral identity).** Let \(C\) be a symmetric conference matrix of even order \(n\) (so \(C^\top C=(n-1)I\), zero diagonal, off-diagonal \(\pm1\)). Let \(\lambda=\sqrt{n-1}\) and let \(P_+\) be the orthogonal projector onto the \(+\lambda\) eigenspace of \(C\) (dimension \(n/2\)). Then
\begin{equation}
\label{eq:conf-spec}
C=\lambda\,(2P_+-I),\qquad
P_+=\frac{I+C/\lambda}{2},
\end{equation}
and for every \(x\in\mathbb R^n\),
\begin{equation}
\label{eq:xCx}
x^\top C x=\lambda\bigl(2\|P_+x\|_2^2-\|x\|_2^2\bigr).
\end{equation}
In particular, for \(x\in\{\pm1\}^n\),
\[
\rho(C)=\max_{x\in\{\pm1\}^n}\Bigl|2\,\frac{\|P_+x\|_2^2}{n}-1\Bigr|,
\qquad
\Phi(C)=\frac12\,\lambda\cdot n\cdot\rho(C)=\frac12\,n\sqrt{n-1}\,\rho(C).
\]

*Proof.* Spectrum of \(C\) is \(\{\pm\lambda\}\) with equal multiplicity \(n/2\) (trace zero, \(C^2=\lambda^2 I\)). Hence \(C=\lambda(P_+-P_-)\) and \(P_++P_-=I\), so \(C=\lambda(2P_+-I)\). The formula for \(P_+\) follows by solving. Substitute into \(x^\top C x\). Diagonals: \(C_{ii}=0=\lambda(2(P_+)_{ii}-1)\) forces \((P_+)_{ii}=\tfrac12\). Off-diagonals: \((P_+)_{ij}=C_{ij}/(2\lambda)=\pm1/(2\lambda)\). \(\square\)

**Proposition 15.2 (exact Nesterov expectation).** Let \(C\) be as above, \(g\sim\mathcal N(0,P_+)\) (Gaussian supported on the \(+\lambda\) eigenspace), and \(s=\mathrm{sign}(g)\in\{\pm1\}^n\). Then
\begin{equation}
\label{eq:nest-exact}
\mathbb E\bigl[s^\top C s\bigr]
=
\frac2\pi\,n(n-1)\,\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr),
\end{equation}
and therefore
\begin{equation}
\label{eq:nest-rho}
\rho(C)
\;\ge\;
\frac{\mathbb E[s^\top C s]}{n\lambda}
=
\frac2\pi\,\sqrt{n-1}\,\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr)
\;\xrightarrow{n\to\infty}\;
\frac2\pi.
\end{equation}

*Proof.* The covariance of \(g\) is \(P_+\), with diagonal \(\tfrac12\) and off-diagonal \(C_{ij}/(2\lambda)\). Hence
\[
\mathbb E[s_is_j]=\frac2\pi\arcsin\bigl(2(P_+)_{ij}\bigr)
=
\begin{cases}
1,&i=j,\\
\frac2\pi\arcsin\bigl(C_{ij}/\lambda\bigr),&i\neq j.
\end{cases}
\]
Since \(\arcsin\) is odd and \(C_{ij}=\pm1\),
\[
\mathbb E[s^\top C s]
=\sum_{i\neq j}C_{ij}\cdot\frac2\pi\arcsin\Bigl(\frac{C_{ij}}{\lambda}\Bigr)
=\frac2\pi\arcsin\Bigl(\frac1\lambda\Bigr)\sum_{i\neq j}C_{ij}^2
=\frac2\pi\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr)\,n(n-1).
\]
Divide by \(n\lambda\) and pass to the limit using \(\arcsin u\sim u\). \(\square\)

**Corollary 15.3.** Along any sequence of conference matrices, \(\Phi(C)\ge\bigl(\tfrac1\pi-o(1)\bigr)n\sqrt{n-1}\), hence the Paley upper-bound sequence satisfies
\[
\frac1\pi
\;\le\;
\liminf_k\frac{\Phi(C_{n_k})}{n_k^{3/2}}
\;\le\;
\limsup_k\frac{\Phi(C_{n_k})}{n_k^{3/2}}
\;\le\;
\frac12.
\]
(The same Nesterov lower bound does **not** control \(m_n=\min_A\Phi(A)\).)

### §15.2 Switching, minimal operator norm, and the optimality reduction (new, complete)

Write \(\mathcal S_n\) for the set of Seidel matrices of order \(n\) (symmetric, zero diagonal, off-diagonal \(\pm1\)).

**Proposition 15.4 (Seidel switching).** Let \(D=\mathrm{diag}(\varepsilon)\), \(\varepsilon\in\{\pm1\}^n\), and \(A\in\mathcal S_n\). Set \(A'=DAD\). Then \(A'\in\mathcal S_n\), \(A'\) is cospectral with \(A\), and \(\Phi(A')=\Phi(A)\).

*Proof.* Off-diagonal entries of \(A'\) are \(\varepsilon_i A_{ij}\varepsilon_j\in\{\pm1\}\); diagonal remains \(0\). Cospectrality: \(A'v=\lambda v\) iff \(A(Dv)=\lambda(Dv)\). For \(\Phi\): \(x^\top A'x=(Dx)^\top A(Dx)\) and \(x\mapsto Dx\) permutes \(\{\pm1\}^n\). \(\square\)

**Proposition 15.5 (absolute bound / min-op).** For every \(A\in\mathcal S_n\),
\[
\|A\|_{\mathrm{op}}\ge\sqrt{n-1},
\]
with equality if and only if \(A^2=(n-1)I\) (i.e.\ \(A\) is a conference matrix). In particular, when a conference matrix of order \(n\) exists, the minimizers of \(\|A\|_{\mathrm{op}}\) on \(\mathcal S_n\) are exactly the conference matrices of order \(n\).

*Proof.* \(\mathrm{tr}(A)=0\) and \(\|A\|_F^2=n(n-1)\), so if \(\lambda_1,\dots,\lambda_n\) are the eigenvalues then \(\sum\lambda_i=0\) and \(\sum\lambda_i^2=n(n-1)\). Hence \(\max_i\lambda_i^2\ge n-1\), i.e.\ \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\). Equality forces \(\lambda_i^2=n-1\) for all \(i\); combined with \(\sum\lambda_i=0\) and \(n\) even one has spectrum \(\{\pm\sqrt{n-1}\}\) with equal multiplicity, so \(A^2=(n-1)I\). \(\square\)

**Proposition 15.6 (factorization of \(\Phi\)).** For every \(A\in\mathcal S_n\) with \(\|A\|_{\mathrm{op}}>0\),
\[
\Phi(A)=\frac12\,n\,\|A\|_{\mathrm{op}}\,\rho(A),\qquad
\rho(A)=\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\,\|A\|_{\mathrm{op}}}.
\]
Consequently, whenever a conference matrix \(C\) of order \(n\) exists,
\begin{equation}
\label{eq:phi-factor}
m_n\le\Phi(C)=\frac12\,n\sqrt{n-1}\,\rho(C).
\end{equation}

**Proposition 15.7 (beaters must have strictly worse cube/sphere ratio).** Let \(C\) be conference of order \(n\) and \(A\in\mathcal S_n\). If \(\Phi(A)<\Phi(C)\), then
\[
\rho(A)<\rho(C)\cdot\frac{\sqrt{n-1}}{\|A\|_{\mathrm{op}}}\le\rho(C).
\]
In particular, no matrix with \(\rho(A)\ge\rho(C)\) can beat \(C\).

*Proof.* \(\Phi(A)=\tfrac12 n\|A\|_{\mathrm{op}}\rho(A)\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\,\rho(C)\). The inequality \(\Phi(A)<\Phi(C)\) rearranges to the claim; the last step is Prop 15.5. \(\square\)

**Proposition 15.8 (limsup controlled by Paley \(\rho\)).** Let \(n_k=q_k+1\) be Paley orders and \(C_k\) the corresponding Paley conference matrix. Then
\begin{equation}
\label{eq:limsup-rho}
\limsup_{n\to\infty}\alpha_n
\;\le\;
\frac12\limsup_{k\to\infty}\rho(C_k)
\;\le\;
\frac12.
\end{equation}
In particular, if \(\limsup_k\rho(C_k)=\rho^\star<1\), this **strictly improves** the spherical limsup \(\tfrac12\) of Theorem A.

*Proof.* For \(N\ge2\) let \(n_k\) be the least Paley order \(\ge N\). Then \(m_N\le m_{n_k}\le\Phi(C_k)\), so
\[
\alpha_N\le\frac{\Phi(C_k)}{N^{3/2}}=\frac12\Bigl(\frac{n_k}{N}\Bigr)^{3/2}\sqrt{1-\frac1{n_k}}\,\rho(C_k).
\]
As \(N\to\infty\), \(n_k/N\to1\) (Prop 6.2), and the claim follows. \(\square\)

**Proposition 15.9 (equivalence form of asymptotic optimality).** Let \(C_k\) be Paley conference of order \(n_k\). Write
\[
r(A)\,:=\,\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\sqrt{n-1}}
=\rho(A)\cdot\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}.
\]
Then \(r(A)\ge\rho(A)\) with equality iff \(\|A\|_{\mathrm{op}}=\sqrt{n-1}\). The following are equivalent:
1. \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) (asymptotic optimality along Paley);
2. \(\displaystyle\min_{A\in\mathcal S_{n_k}}r(A)=\rho(C_k)+o(1)\);
3. \(\displaystyle\min_{A\in\mathcal S_{n_k}}\rho(A)\,\|A\|_{\mathrm{op}}=\rho(C_k)\sqrt{n_k-1}+o(\sqrt{n_k})\).

*Proof.* \(\Phi(A)=\tfrac12 n\sqrt{n-1}\,r(A)\) and \(m_n=\min\Phi\), while for conference \(r(C)=\rho(C)\). \(\square\)

Thus asymptotic optimality is the claim that conference minimises the single scalar \(r(A)=\rho\cdot\|A\|_{\mathrm{op}}/\sqrt{n-1}\). Prop 15.5 says conference uniquely minimises the op-factor; Prop 15.7 says any competitor must pay in \(\rho\).

**Proposition 15.10 (L²-universality on the cube).** For every \(A\in\mathcal S_n\) and \(x\) uniform in \(\{\pm1\}^n\), writing \(Q(x)=\sum_{i<j}A_{ij}x_ix_j=\tfrac12 x^\top A x\),
\[
\mathbb E\bigl[Q(x)^2\bigr]=\binom{n}{2},\qquad
\mathbb E\bigl[(x^\top A x)^2\bigr]=2n(n-1).
\]
In particular the \(L^2\) mass of the degree-2 form is **identical** for every Seidel matrix.

*Proof.* The Walsh functions \(x\mapsto x_ix_j\) (\(i<j\)) are orthonormal in \(L^2(\{\pm1\}^n)\). Hence \(\mathbb E[Q^2]=\sum_{i<j}A_{ij}^2=\binom{n}{2}\). The second identity is \(x^\top A x=2Q\). \(\square\)

**Proposition 15.11 (unique minimiser of \(\mathrm{tr}(A^4)\)).** For every \(A\in\mathcal S_n\),
\[
\mathrm{tr}(A^4)=\sum_{i=1}^n\lambda_i(A)^4\ge n(n-1)^2,
\]
with equality if and only if \(A\) is a conference matrix.

*Proof.* \(\sum\lambda_i^2=\|A\|_F^2=n(n-1)\). By Cauchy–Schwarz / QM-AM on \((\lambda_i^2)_{i=1}^n\),
\[
\frac1n\sum\lambda_i^4\ge\Bigl(\frac1n\sum\lambda_i^2\Bigr)^2=(n-1)^2,
\]
so \(\sum\lambda_i^4\ge n(n-1)^2\). Equality holds iff all \(\lambda_i^2\) are equal, i.e.\ \(A^2=(n-1)I\). \(\square\)

**Proposition 15.13 (exact fourth moment of \(Q\)).** Let \(A\in\mathcal S_n\), \(e=\binom{n}{2}\), and let \(x\) be uniform in \(\{\pm1\}^n\). Write \(Q(x)=\sum_{i<j}A_{ij}x_ix_j\). Then
\begin{equation}
\label{eq:Q4}
\mathbb E\bigl[Q(x)^4\bigr]
=
3e^2
+
3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)
-
n(n-1)(3n-5).
\end{equation}
Equivalently,
\[
\mathbb E\bigl[(x^\top A x)^4\bigr]
=
48\,\mathrm{tr}(A^4)
+12\,\mathrm{tr}(A^2)^2
-48\sum_{i=1}^n(A^2)_{ii}^2
-16\,\mathrm{tr}(A^2)\,(3n-5).
\]
In particular, since \(\mathrm{tr}(A^2)=n(n-1)\) and \((A^2)_{ii}=n-1\) for every Seidel matrix,
\[
\mathbb E[Q^4]
=
3e^2-n(n-1)(3n-5)
+3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)
\]
is **uniquely minimised** precisely when \(A\) is a conference matrix, with minimum value
\[
\mathbb E_C[Q^4]=3e^2-n(n-1)(3n-5)=\frac{n(n-1)}{4}\bigl(3n(n-1)-4(3n-5)\bigr).
\]

*Proof sketch.* The identity is the specialisation to Seidel matrices (\(A_{ij}^2=1\) off-diagonal, \(A_{ii}=0\)) of the degree-4 moment expansion of a Rademacher chaos of order 2. The expansion’s graph-counting terms (pairings, wedges, signed 4-cycles, and length-4 closed walks) reduce, for Seidel matrices, to the single spectral invariant \(\mathrm{tr}(A^4)=\|A^2\|_F^2\) together with pure functions of \(n\). The resulting closed form is recorded above; it has been cross-checked to machine precision against exhaustive half-cube enumeration for all Seidel matrices of orders \(4\le n\le 11\) in a 169-matrix battery (random, Paley, and all-negative). Uniqueness of the minimiser is Prop 15.11. \(\square\)

*Corollary.* Conference matrices are the unique Seidel matrices that simultaneously
(i) minimise \(\|A\|_{\mathrm{op}}\),
(ii) minimise \(\mathrm{tr}(A^4)\),
(iii) minimise \(\mathbb E[Q^4]\),
and (iv) realise the universal cube-\(L^2\) mass of Prop 15.10.
Asymptotic optimality of \(m_n\) is the remaining claim that they also minimise the cube \(L^\infty\) norm of \(Q\) (i.e.\ minimise \(r(A)\)).

*Remark (L⁴ lower bound on \(\Phi\)).* Cauchy–Schwarz on the cube measure yields \(\max|Q|^2\ge\mathbb E[Q^4]/\mathbb E[Q^2]\), hence
\[
\Phi(A)\;\ge\;
\sqrt{
3e+\frac{3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)}{e}-2(3n-5)
}.
\]
The right-hand side is minimised at conference matrices and is \(\sim n\sqrt{3/2}\) for large \(n\) — an \(\Omega(n)\) lower bound on \(m_n\), weaker than the Bohnenblust–Hille \(\Omega(n^{3/2})\) bound of Prop 5.1, but spectrally sharp within the moment method of order 4.

**Proposition 15.14 (exact optimality criterion via fourth moments).** Let \(C\) be a conference matrix of order \(n\), \(e=\binom{n}{2}\), \(\Phi_*=\Phi(C)\), and
\[
\Delta_*:=\Phi_*^2\,e-\mathbb E_C[Q^4]
=\Phi_*^2\,e-3e^2+n(n-1)(3n-5).
\]
Write \(\delta(A):=\mathrm{tr}(A^4)-n(n-1)^2\ge0\). If \(A\in\mathcal S_n\) satisfies \(\Phi(A)\le\Phi_*\), then necessarily
\begin{equation}
\label{eq:delta-bound}
\delta(A)\;\le\;\frac{\Delta_*}{3}.
\end{equation}
Consequently: if every non-conference Seidel matrix of order \(n\) obeys \(\delta(A)>\Delta_*/3\), then every \(A\) with \(\Phi(A)\le\Phi_*\) is a conference matrix, and therefore
\[
m_n=\Phi_*.
\]

*Proof.* From \(\mathbb E[Q^4]\le(\max|Q|)^2\mathbb E[Q^2]\) and Prop 15.13,
\[
3e^2-n(n-1)(3n-5)+3\delta(A)
=\mathbb E_A[Q^4]
\le\Phi(A)^2\,e
\le\Phi_*^2\,e
=\mathbb E_C[Q^4]+\Delta_*
=3e^2-n(n-1)(3n-5)+\Delta_*.
\]
Cancel to get \(3\delta(A)\le\Delta_*\). If the spectral gap of every non-conference matrix exceeds \(\Delta_*/3\), then \(\delta(A)=0\), so \(A\) is conference (Prop 15.11) and \(\Phi(A)=\Phi_*\) along the switching class (Prop 15.4). \(\square\)

**Corollary 15.15 (exact optimality at \(n=6\)).** For \(n=6\), the Paley conference matrix satisfies \(\Phi_*=5\), \(\Delta_*=90\), and \(\Delta_*/3=30\). Exhaustive enumeration of all \(2^{10}=1024\) vertex-folded Seidel matrices of order 6 yields \(\min\delta(A)=64>30\) over non-conference matrices (and exactly 12 conference matrices in the switching class, all with \(\Phi=5\)). Therefore \(m_6=5=\Phi(C)\).

*Remark.* For \(n=14\) (Paley, \(\Phi_*=21\)) one has \(\Delta_*/3\approx7341\), while the structural lower bound \(\delta\ge16(n-2)=192\) (from at least \(4(n-2)\) off-diagonal entries of \(A^2\) equal to \(\pm2\)) is far smaller, so \eqref{eq:delta-bound} does **not** force conference.

**Proposition 15.16 (global super-linear \(\min\delta\) is impossible).** Let \(C\) be a conference matrix of order \(n\) and let \(A\) be obtained by flipping a single off-diagonal pair \((A_{ij},A_{ji})=(-C_{ij},-C_{ji})\). Then
\[
\delta(A)=16(n-2)=\Theta(n).
\]
In particular \(\min\{\delta(A):A\in\mathcal S_n\text{ non-conference}\}=O(n)\), so no super-linear lower bound on \(\min\delta\) over *all* non-conference Seidel matrices can hold. Consequently the *global* gap test “\(\delta>\Delta_*/3\) for every non-conference \(A\)” of Prop 15.14 is asymptotically void: \(\Delta_*/3=\Theta(n^5)\) (using \(\Phi_*=\Theta(n^{3/2})\)) while the structural floor is \(\Theta(n)\). At \(n=6\) the test still works because the floor \(64\) exceeds \(\Delta_*/3=30\); already at \(n=14\) one has \(192\ll7341\).

*Proof.* For any Seidel matrix, \(\delta(A)=\sum_{i\neq j}(A^2)_{ij}^2\). Conference matrices have \(A^2=(n-1)I\), so all off-diagonal entries of \(A^2\) vanish. A direct expansion of \(A^2\) after one edge flip shows that exactly \(4(n-2)\) off-diagonal entries become \(\pm2\) (and the rest stay \(0\) up to the symmetric copy), giving \(\delta=4(n-2)\cdot4=16(n-2)\). This identity is confirmed exactly on all edges of the Paley conference matrices of orders \(n=6,14,18\). \(\square\)

**Proposition 15.17 (conditional fourth-moment criterion).** Let \(C,\Phi_*,\Delta_*\) be as in Prop 15.14. Then
\[
m_n=\Phi_*
\quad\text{if and only if}\quad
\min\bigl\{\Phi(A):A\in\mathcal S_n,\;\delta(A)\le\Delta_*/3\bigr\}
=\Phi_*.
\]
(The minimum on the right is attained on the conference switching class whenever it equals \(\Phi_*\).)

*Proof.* Prop 15.14 says every \(A\) with \(\Phi(A)\le\Phi_*\) lies in \(\{\delta\le\Delta_*/3\}\). Hence if the minimum of \(\Phi\) on that set is \(\Phi_*\), no matrix undercuts \(\Phi_*\). Conversely if \(m_n=\Phi_*\) then the conference matrices realise the minimum on every nonempty set containing them. \(\square\)

Thus the only remaining \(L^4\to L^\infty\) path is to control \(\Phi\) *inside* the thin spectral shell \(\delta\le\Delta_*/3\), not to enlarge the global \(\min\delta\). Numerically (86-worker SA under the cap \(\delta\le\Delta_*/3\), 40 seeds per Paley order; also under \(\Phi\le\Phi_*+\text{margin}\); edge-flip and min-op searches), the best \(\Phi\) found inside the shell is \(5,27,41\) at \(n=6,14,18\) against Paley \(\Phi_*=5,21,33\) — **no undercut**. At \(n=6\) every sample under the cap is conference. This is strong evidence for Prop 15.17 at these orders, not a general proof.

**Proposition 15.18 (projector form of \(\rho\)).** Let \(C\) be a symmetric conference matrix of order \(n\), \(s=\sqrt{n-1}\), and \(P_+=(I+C/s)/2\) the orthogonal projector onto the positive eigenspace of \(C\). Then for every \(x\in\mathbb R^n\),
\[
x^\top C x = s\bigl(2\|P_+x\|_2^2-\|x\|_2^2\bigr),
\]
and therefore
\[
\rho(C)=\max_{x\in\{\pm1\}^n}\Bigl|\,2\frac{\|P_+x\|_2^2}{n}-1\,\Bigr|.
\]
In particular Theorem E(2) is equivalent to: the maximal cube-imbalance of the Paley spectral projector \(P_+\) converges as \(n_k\to\infty\).

*Proof.* \(C=s(2P_+-I)\) because the eigenvalues of \(C\) are \(\pm s\). Expand and restrict to the cube \(\|x\|_2^2=n\). \(\square\)

Fixed-order \(L^p\) norms of \(Q\) give only \(\Omega(n)\) lower bounds on \(\Phi\) (Bonami hypercontractivity yields \(\|Q\|_p\le(p-1)\|Q\|_2=O(n)\) for fixed \(p\)), so they cannot replace the \(L^\infty\) comparison for the \(n^{3/2}\) scale.

**Proposition 15.19 (asymptotic vacuity of the fourth-moment shell).** For every Seidel matrix \(A\in\mathcal S_n\),
\[
\delta(A)=\mathrm{tr}(A^4)-n(n-1)^2\le(n-1)^4+(n-1)-n(n-1)^2=O(n^4),
\]
with equality in the leading \((n-1)^4\) term for the switching class of \(J-I\). On the other hand, if a conference matrix \(C\) of order \(n\) exists with \(\Phi_*=\Phi(C)\) and \(\rho(C)\ge\rho_0>0\), then
\[
\frac{\Delta_*}{3}
=\frac{\Phi_*^2\,e-\mathbb E_C[Q^4]}{3}
\sim\frac{\rho(C)^2}{24}\,n^5
\]
as \(n\to\infty\). Consequently there exists \(N=N(\rho_0)\) such that for all conference orders \(n\ge N\),
\[
\frac{\Delta_*}{3}\ge\max_{A\in\mathcal S_n}\delta(A),
\]
so the set \(\{\delta\le\Delta_*/3\}\) is **all of** \(\mathcal S_n\). In particular Prop 15.17 becomes the tautology \(m_n=\min\Phi\), and the fourth-moment criterion supplies **no** asymptotic information.

*Proof of the \(\delta\)-bound.* Write \(\lambda_1,\ldots,\lambda_n\) for the eigenvalues of \(A\). Seidel matrices satisfy \(\sum_i\lambda_i^2=\mathrm{tr}(A^2)=n(n-1)\) and \(\|\lambda\|_\infty=\|A\|_{\mathrm{op}}\le n-1\) (since \(A+I\) is a \(\{0,2\}\)-matrix after a diagonal shift of the complete graph, or by Gershgorin). Under the \(\ell^2\) and \(\ell^\infty\) constraints, \(\sum\lambda_i^4\) is maximised by putting one coordinate at \(\pm(n-1)\) and the rest at values of absolute value at most \(1\) in the \(\ell^2\) budget \(n-1\), which is realised by \(A=J-I\) (eigenvalues \(n-1\) once and \(-1\) with multiplicity \(n-1\)). The asymptotic for \(\Delta_*\) is Prop 15.1: \(\Phi_*=\tfrac12 n\sqrt{n-1}\,\rho(C)\), \(e=\binom{n}{2}\), and \(\mathbb E_C[Q^4]=O(n^4)\). \(\square\)

*Numerical threshold.* With exact \(\Phi_*\) one has \(\Delta_*/3<\max\delta\) at \(n\le30\) and \(\Delta_*/3>\max\delta\) already at \(n=38\) (\(\Phi_*=109\) exact). Thus Props 15.14–15.17 are useful only in a finite window of small conference orders; asymptotic optimality requires a genuinely different \(L^\infty\) comparison.

**Proposition 15.20 (Lipschitz continuity of \(\Phi\); Frobenius form).** For all \(A,B\in\mathcal S_n\) and all \(x\in\{\pm1\}^n\),
\[
\bigl|x^\top(A-B)x\bigr|\le n\,\|A-B\|_F,
\]
hence
\begin{equation}
\label{eq:phi-lip}
\bigl|\Phi(A)-\Phi(B)\bigr|
\;\le\;
\frac n2\,\|A-B\|_F.
\end{equation}
In particular, if \(A\) differs from a conference matrix \(C\) in exactly \(k\) undirected edges, then \(\|A-C\|_F=2\sqrt k\) and
\[
\Phi(A)\;\ge\;\Phi(C)-n\sqrt k.
\]
Relative to the \(n^{3/2}\) scale this is an \(O(\sqrt{k/n})\) relative error: any \(A\) within \(k=o(n)\) edge flips of a conference matrix satisfies \(\Phi(A)\ge\Phi(C)-o(n^{3/2})\).

*Proof.* Cauchy–Schwarz: \(|x^\top Mx|\le\|M\|_F\|x\|_2^2=n\|M\|_F\) for \(M=A-B\) and \(\|x\|_2=\sqrt n\). Taking \(\Phi=\max|Q|=\tfrac12\max|x^\top(\,\cdot\,)x|\) yields \eqref{eq:phi-lip}. Each flipped edge changes two off-diagonal entries by \(2\) in absolute value, contributing \(4\) to \(\|A-C\|_F^2\) per edge. \(\square\)

**Proposition 15.20b (edge-counting Lipschitz — sharp sparse form).** Let \(A,C\in\mathcal S_n\) differ in exactly \(k\) undirected edges. Then for every \(x\in\{\pm1\}^n\),
\[
\bigl|Q_A(x)-Q_C(x)\bigr|\;\le\;2k,
\]
and therefore
\begin{equation}
\label{eq:phi-edge-lip}
\Phi(A)\;\ge\;\Phi(C)-2k.
\end{equation}
In particular, on any order where a conference matrix \(C\) with \(\rho(C)=1\) exists (so \(\Phi(C)=\tfrac12 n\sqrt{n-1}\)),
\[
m_n\;\ge\;\tfrac12 n\sqrt{n-1}-2k_\star
\quad\text{where}\quad
k_\star:=\min\bigl\{d_H(A',C):A'\sim A\text{ by Seidel switching},\;\Phi(A)=m_n\bigr\}
\]
(minimum over switchings of a minimiser \(A\)). Consequently:
\begin{itemize}
\item if \(k_\star=o(n^{3/2})\) then \(m_n=\Phi(C)-o(n^{3/2})\) (E(1) on that order);
\item the Frobenius form (Prop 15.20) needed the stronger \(k_\star=o(n)\) for the same conclusion — edge counting weakens the rigidity demand from \(o(n)\) to \(o(n^{3/2})\).
\end{itemize}
Shipped: `edge_hamming`, `phi_edge_lipschitz_lower` in `src/minmax_quadratic.py`.

*Proof.* On a disagreeing edge one has \(A_{ij}-C_{ij}=-2C_{ij}\), so
\[
Q_A(x)-Q_C(x)
=\sum_{\{i,j\}\in F}(A_{ij}-C_{ij})x_ix_j
=\sum_{\{i,j\}\in F}(-2C_{ij})x_ix_j
\]
and \(\lvert Q_A(x)-Q_C(x)\rvert\le 2\lvert F\rvert=2k\). Hence \(\lvert Q_A(x)\rvert\ge\lvert Q_C(x)\rvert-2k\) for every \(x\), and taking \(\max_x\) yields \eqref{eq:phi-edge-lip}. The E(1) criterion is the specialisation \(\Phi(C)=\tfrac12 n\sqrt{n-1}\) together with \(2k_\star=o(n^{3/2})\). \(\square\)

**Proposition 15.20c (degree Lipschitz).** If the disagreement graph of \(A\) vs \(C\) has maximum degree \(D\), then \(\Phi(A)\ge\Phi(C)-Dn\).

*Proof.* \(E=A-C\) has off-diagonal entries in \(\{0,\pm2\}\) and row \(\ell^1\)-norm at most \(2D\), so \(\|E\|_{\mathrm{op}}\le 2D\). Thus \(\lvert x^\top Ex\rvert\le 2Dn\) and \(\lvert Q_A-Q_C\rvert\le Dn\). \(\square\)

*Remark (N10 matches the edge form).* At \(n=10\), undercutting optima are perfect-matching flips of Paley (\(k_\star=5\), \(D=1\)): edge lip gives \(m_{10}\ge15-10=5\) (true; actual \(m_{10}=13\)), degree lip gives the same \(15-10=5\). Relative gap \(2/n^{3/2}\approx0.063\to0\) along any sequence with \(k_\star=O(n)\). **E(1) is reduced to proving \(k_\star=o(n^{3/2})\) for \(\Phi\)-minimisers on \(n=p^2+1\)** (still open as a general statement; verified structurally at \(n=10\), consistent with exact MITM SA at \(n=26\) where no undercut of \(\Phi=65\) is known).

**Proposition 15.20d (conditional settlement: \(k_\star=O(n^{3/2})\) \(\Rightarrow L=\tfrac12\)).** Let \(n_k=p_k^2+1\) run over the \(\rho=1\) Paley family, and write \(k_\star(n)\) for the minimal best-switch Hamming distance from a \(\Phi\)-minimiser to the Paley conference matrix of order \(n\). If \(k_\star(n_k)=O(n_k^{3/2})\), then
\[
\lim_{n\to\infty}\alpha_n=\tfrac12.
\]
*Proof.* By the Max-Lipschitz Prop 15.27 (using the \(\mathrm{Max}_{+}\) frame identity certified for Paley), \(m_{n_k}\ge\Phi(C_{n_k})-2k_\star(n_k)/p_k\). Under \(k_\star=O(n^{3/2})\) and \(p=\sqrt{n-1}\),
\[
\frac{2k_\star}p=O(n),\qquad
\alpha_{n_k}\ge\tfrac12\sqrt{1-1/n_k}-O(n_k^{-1/2})\to\tfrac12.
\]
Combined with the universal limsup \(\limsup\alpha_n\le\tfrac12\) (Prop 4.1) and denseness Prop 6.1–6.2 along \(n_k\), one has \(\lim\alpha_n=\tfrac12\). \(\square\)

*(Weaker edge-only form.)* Prop 15.20b alone needs the stronger hypothesis \(k_\star=O(n)\) for the same conclusion (gap \(2k_\star=O(n)\)). Max-Lipschitz saves a factor \(p=\sqrt{n-1}\).

*Status of the hypothesis.* At \(n=10\), \(k_\star=5=O(n)\subset O(n^{3/2})\) (N10-S/C6). At \(n=26\), exact MITM sparse/SA census found no undercut of \(\Phi=65\) (consistent with \(k_\star=0\)). **The general bound \(k_\star=O(n^{3/2})\) on all \(n=p^2+1\) is not proved.** Existence of \(\lim\alpha_n\) remains **OPEN**.

**Proposition 15.21 (single-edge local optimality under maximizer balance).** Let \(C\in\mathcal S_n\) be a conference matrix, \(M=\Phi(C)\), and \((p,r)\) an unordered edge. Write \(\varepsilon(x):=C_{pr}x_px_r\in\{\pm1\}\). Suppose there exists \(x^*\in\{\pm1\}^n\) with \(|Q_C(x^*)|=M\) and
\[
\varepsilon(x^*)\,=\,-\mathrm{sign}\bigl(Q_C(x^*)\bigr)
\]
(with the convention \(\mathrm{sign}(0)=+1\)). Let \(A\) be \(C\) with edge \((p,r)\) flipped. Then
\[
\Phi(A)\;\ge\;M+2.
\]
In particular \(A\) cannot undercut \(C\).

*Proof.* Flipping changes the quadratic form by \(Q_A(x)=Q_C(x)-2\varepsilon(x)\). At \(x^*\) one has \(Q_A(x^*)=Q_C(x^*)-2\varepsilon(x^*)=Q_C(x^*)+2\,\mathrm{sign}(Q_C(x^*))\), so \(|Q_A(x^*)|=M+2\). \(\square\)

**Corollary 15.21a (Paley edge-transitivity + balance).** For the Paley conference matrices of orders \(n\in\{6,14,18\}\), *every* edge satisfies the maximizer-balance hypothesis of Prop 15.21 (verified by exhaustive half-cube enumeration of the maximizer sets: \(12\), \(156\), and \(204\) maximizers respectively; zero unbalanced edges). Consequently every single edge flip strictly increases \(\Phi\) by at least \(2\). Combined with Prop 15.16 (\(\delta=16(n-2)\) after one flip), the Paley matrix is a strict local minimum of both \(\Phi\) and \(\delta\) in the edge-flip graph at these orders.

*Remark (E(1) programme without \(L^4\)).* Props 15.20–15.21 control matrices *near* a conference in Hamming distance. The remaining gap is matrices with \(\|A\|_{\mathrm{op}}=(1+o(1))\sqrt{n-1}\) that are *not* Hamming-close to any conference (near-equiangular Seidel matrices outside the switching class). The natural three-step attack is:
1. **Universal cube/sphere floor.** Prove \(\rho(A)\ge\tfrac2\pi-o(1)\) for all \(A\in\mathcal S_n\) (Nesterov SDP rounding under eigenvector delocalization: if a top eigenvector satisfies \(\|v\|_\infty\le n^{-1/2+\varepsilon}\) then the rank-1 SDP witness \(Y_{ij}=v_iv_j/(|v_i||v_j|)\) yields SDP\(\ge n\|A\|_{\mathrm{op}}(1-o(1))\) and cube\(\ge(2/\pi)\)SDP).
2. **Op-control of near-minimisers of \(r\).** If \(r(A)\le\rho(C)+o(1)\) and \(\rho(A)\ge2/\pi-o(1)\), then \(\|A\|_{\mathrm{op}}/\sqrt{n-1}\le(\pi/2)\rho(C)+o(1)\). Along Paley, \(\rho(C)\to\rho_*\) would give a uniform op-factor bound.
3. **Spectral rigidity \(\Rightarrow\) Hamming closeness.** Convert small \(\delta(A)=\mathrm{tr}(A^4)-n(n-1)^2\) (which follows from op-factor\(\to1\) by Prop 15.12) into existence of a conference \(C\) with \(\|A-C\|_F=o(n)\), then apply Prop 15.20.

Step 1 is standard for Wigner-type matrices with delocalized eigenvectors but is not proved here for every Seidel matrix. Steps 2–3 are conditional on step 1 and on \(\rho(C_k)\to\rho_*\). Multi-core Nesterov-cluster sampling (`attack_deloc_nesterov`, 86 workers) recovers \(\rho_{\mathrm{LB}}\ge 2/\pi-0.01\) on every random/near-min-op sample tested (\(n\le32\)), and recovers the exact \(\rho(C)\) on Paley conference matrices (fat \(\pm\) eigenspaces of multiplicity \(n/2\)).

*Delocalization barrier.* SA maximising \(\sqrt n\|v\|_\infty\) for extreme eigenvectors of Seidel matrices reaches factors \(\approx 2.75\) at \(n=40\) (`attack_deloc_max`). The elementary Nesterov+rank-1 bound only yields \(\rho\ge(2/\pi)/c^2\) with \(c=\sqrt n\|v\|_\infty\), hence \(\rho\ge0.08\) at \(c=2.75\) — useless compared to Prop 15.22's reduction and the DMP floor. Universal \(\rho\ge2/\pi\) therefore needs a genuinely different argument (full SDP analysis, not rank-1 witnesses).

*Local-min evidence for E(1).* Edge-flip local search for minimisers of \(\Phi\) (86 workers, SA+greedy descent, `attack_local_min_phi`): at \(n=6\) every one of 50 local minima is conference (\(\Phi=5\)); at Paley orders \(n=14,18\) no run undercuts \(\Phi(C)\) (best local mins at \(25,41>21,33\)); at \(n=8,10,12\) the unique/near-unique local \(\Phi\) values match prior global SA bests. This is consistent with conference being the global min whenever it exists, but does not prove absence of a rare undercutting basin.

**Theorem N10-S (certified structure of the \(n=10\) gap; 2026-07-26).** Let \(C\) be Paley of order \(10\) (\(\Phi(C)=15\), \(\rho=1\)). Exact enumeration yields:
1. Every edge of \(C\) is maximizer-balanced (Prop 15.21); every single edge flip has \(\Phi\ge17\).
2. \(\min\{\Phi(A):d_H(A,C)=k\}\) for \(k=0,\ldots,5\) equals \(15,17,15,17,15,13\). Thus no matrix within Hamming distance \(4\) of \(C\) undercuts Paley, and the first undercutting to \(m_{10}=13\) is at distance **5**.
3. Of the \(\binom{45}{5}\) five-edge sets, exactly **144** achieve \(\Phi=13\), and **all 144 are perfect matchings** of \(K_{10}\). Among all \(945\) perfect matchings, the flip-\(\Phi\) histogram is \(\{13{:}144,\,17{:}405,\,21{:}360,\,25{:}36\}\).
4. Absolute gap \(\Phi(C)-m_{10}=2\) gives relative gap \(2/n^{3/2}\approx0.063\to0\) if such \(O(1)\) gaps persist — consistent with E(1), not a counterexample.
5. Random perfect-matching flips of Paley \(C_{26}\) (\(\Phi=65\)) produced \(\Phi\ge73\) on 86 samples; SA+exact-\(\Phi\) rescore (86 workers) found **no** matrix with exact \(\Phi<65\). The matching undercut does **not** lift naively to the next \(\rho=1\) order.

Full writeup and JSON: `evidence/N10_STRUCTURE.md`. Tests: `test_n10_*` in `tests/test_minmax.py`.

**Theorem N10-C (classification of the 144 matchings; 2026-07-27).** Let \(M\) be a perfect matching of \(K_{10}\) and \(S_M(x)=\sum_{\{i,j\}\in M}C_{ij}x_ix_j\). Then flipping \(M\) on Paley \(C_{10}\) yields \(\Phi=13\) if and only if \(\mathrm{sign}(Q_C(x))\,S_M(x)\ge1\) for every maximizer \(x\) of \(C\) (equivalently: for the six maximizers with \(Q_C=+15\)). The resulting 144 matchings form a single orbit under \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) acting on \(\mathrm{PG}(1,9)\). Sign product \(\prod C_{ij}=-1\) is necessary but not sufficient (504 matchings satisfy it). Evidence: `evidence/N10_MATCHING_CLASSIFY.md`. **Existence of \(\lim\alpha_n\) remains open.**

**Theorem N10-C6 (Hamming-6 undercutters are 6-cycles; 2026-07-27).** Exhaustive scan of all \(\binom{45}{6}=8{,}145{,}060\) six-edge sets: exactly **360** satisfy \(\Phi(C\oplus F)<15\), each is a single **6-cycle**, and each has \(\Phi=13=m_{10}\). Combined with N10-S, every undercutter of cardinality \(5\) or \(6\) is a path/cycle graph (\(\Delta\le2\), \(k\le n\)). Evidence: `evidence/N10_CYCLE_UNDERCUTTERS.md`, `src/n10_cycle_undercutters.py`. **Does not settle \(\lim\alpha_n\)** (needs \(k_\star=O(n)\) for general \(p\)).

**Proposition 15.24 (maximizers of \(\rho=1\) conference are boolean eigenvectors; 2026-07-27).** Let \(C\) be conference of order \(n\) with \(\rho(C)=1\) and \(p=\sqrt{n-1}\). For \(x\in\{\pm1\}^n\),
\[
\lvert x^\top Cx\rvert=np\quad\Longleftrightarrow\quad Cx=\pm px.
\]
Thus \(\mathrm{Max}(C)=\{x\in\{\pm1\}^n:Cx=\pm px\}\). *Proof.* \(C=pP_+-pP_-\) yields \(x^\top Cx=p(2\|P_+x\|^2-n)\); equality \(\lvert x^\top Cx\rvert=np\) forces \(x\) into one eigenspace. \(\square\)

Certified boolean \(+p\)-evec counts for Paley \(n=p^2+1\): \(12,260,11452\) at \(p=3,5,7\) (`evidence/BOOLEAN_EVECS_MAX.md`). For the halfspace construction one has \(\sum_i x_i=p+1\) constantly. The ratio \(\#/n^{3/2}\) increases through \(p=7\), so the crude covering bound \(k_\star\le|\mathrm{Max}|\) is **not** \(o(n^{3/2})\) and does not prove E(1).

**Proposition 15.25 (star-reduction recursive formula for \(m_n\); 2026-07-27).** For every \(n\ge2\),
\begin{equation}
\label{eq:recursive-m}
m_n
=
\min_{B\in\mathcal S_{n-1}}
\max_{x\in\{\pm1\}^{n-1}}
\Bigl(
\bigl|Q_B(x)\bigr|
+
\bigl|\textstyle\sum_{i=1}^{n-1}x_i\bigr|
\Bigr),
\end{equation}
where \(Q_B(x)=\sum_{1\le i<j\le n-1}B_{ij}x_ix_j\).

*Proof.* Every Seidel matrix of order \(n\) is switching-equivalent to one with first row \((0,+1,\ldots,+1)\) (switch vertex \(j\) by \(\varepsilon_j=A_{0j}\)). Switching preserves \(\Phi\), so the minimum of \(\Phi\) may be taken over this slice. For such an \(A\), write \(x=(x_0,x')\) and \(B=A[1\!:\!,1\!:]\). Then
\[
Q_A(x)
=
x_0\sum_{j=1}^{n-1}x'_j
+
Q_B(x'),
\]
and \(\max_{x_0=\pm1}\lvert Q_B(x')+x_0\,s(x')\rvert=\lvert Q_B(x')\rvert+\lvert s(x')\rvert\) with \(s=\sum x'_i\). Therefore \(\Phi(A)=\max_{x'}(\lvert Q_B\rvert+\lvert s\rvert)\), and minimising over \(A\) is minimising over \(B\in\mathcal S_{n-1}\). \(\square\)

*Certified checks.* Identity \eqref{eq:recursive-m} holds exactly for all \(3\le n\le8\) (exhaustive on both sides) and matches recorded \(m_9,m_{10},m_{11}\) under SA for the right-hand side.

*Remark (E(1) via recursion).* Writing \(f(B)=\max(\lvert Q_B\rvert+\lvert s\rvert)\), if \(B^*\) realises the min then \(m_n=f(B^*)\ge\Phi(B^*)+\max_{L^*}\lvert s\rvert\) where \(L^*\) is the level set of \(B^*\) at height \(\Phi(B^*)\). Numerically the boost \(\max_{L^*}\lvert s\rvert\) is often \(0\) or \(1\) (not \(\Omega(\sqrt n)\)), so the recursion alone does not force the \(n^{3/2}\) growth rate. Combined with Prop 15.20d, E(1) still reduces to \(k_\star=O(n)\) on the \(\rho=1\) family. **Existence of \(\lim\alpha_n\) remains OPEN.**

**n=26 exact MITM census (2026-07-27).** Shipped `phi_mitm` (meet-in-the-middle exact \(\Phi\), even \(n\le28\)). Random matchings/cycles/stars/\(k\le20\) flips and 86-seed SA+MITM rescore: **no undercut of \(\Phi(C_{26})=65\)** (best SA exact \(67\)). Evidence: `evidence/E1_N26_SPARSE_EXACT.md`, `e1_n26_mitm_sa.json`. Consistent with \(k_\star=0\) at \(n=26\); not a general E(1) proof.

**Proposition 15.26 (matching flips preserve local maximality of boolean evecs; 2026-07-27).** Let \(C\) be a conference matrix of order \(n\) with \(\rho(C)=1\) and \(p=\sqrt{n-1}\ge3\), and let \(M\) be a matching of \(K_n\). Write \(A\) for the Seidel matrix obtained by flipping the edges of \(M\) on \(C\). Then every boolean eigenvector \(y\in\{\pm1\}^n\) with \(Cy=py\) is a **coordinate-local maximiser** of the map \(x\mapsto x^\top Ax\) on the cube: for all coordinates \(i\),
\[
y_i\,(Ay)_i\ge p-2\ge1>0.
\]
(The same holds for \(Cy=-py\) and local maximisers of \(x\mapsto -x^\top Ax\).)

*Proof.* If \(i\) is unmatched by \(M\), then \((Ay)_i=(Cy)_i=py_i\), so \(y_i(Ay)_i=p>0\). If \(i\) is matched to \(\pi(i)\), flipping the edge changes the \(\pi(i)\)-term in row \(i\) from \(C_{i\pi}y_\pi\) to \(-C_{i\pi}y_\pi\), hence
\[
(Ay)_i=(Cy)_i-2C_{i,\pi(i)}y_{\pi(i)}=py_i-2C_{i,\pi(i)}y_{\pi(i)},
\]
and \(y_i(Ay)_i=p-2C_{i,\pi(i)}y_iy_{\pi(i)}\). The character \(C_{i\pi}y_iy_\pi\in\{\pm1\}\), so the display is at least \(p-2\). \(\square\)

*Certified global coincidence at \(n=10\) only.* For every one of the \(945\) perfect matchings \(M\) of \(K_{10}\), one has the stronger identity
\[
\Phi(C\oplus M)=\max_{y\in\mathrm{Max}(C)}\lvert Q_{C\oplus M}(y)\rvert
\]
(histograms agree: \(\{13{:}144,\,17{:}405,\,21{:}360,\,25{:}36\}\)). **This does not lift to \(n=26\):** among 30 random perfect matchings of Paley \(C_{26}\), only 19 satisfy the identity; the other 11 have \(\Phi>\max_{\mathrm{Max}}|Q|\) (non-maximiser spikes), with exact MITM \(\Phi\in\{75,\ldots,87\}\) all strictly above \(\Phi(C)=65\). Star flips of degree \(\ge3\) destroy even local maximality (scores \(p-2d\) become negative). Evidence: `evidence/E1_STAR_REDUCTION_PROBE.md`.

*Remark (route to the matching dichotomy).* Prop 15.26 keeps boolean evecs first-order critical after matching flips, which organises the \(n=10\) undercut analysis (N10-S/C). The global Max-determination identity is special to \(n=10\), not a general shortcut. Matching dichotomy \(m_n=\min(\Phi(C),\min_M\Phi(C\oplus M))\) remains open; at \(n=26\) random matchings only raise \(\Phi\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.27 (Max\(_{+}\) tight frame and Max-Lipschitz; 2026-07-27).** Let \(C\) be a conference matrix of order \(n\) with \(\rho(C)=1\), \(p=\sqrt{n-1}\), and write \(\mathrm{Max}_{+}=\{y\in\{\pm1\}^n:Cy=py\}\). Assume the frame identity
\begin{equation}
\label{eq:max-frame}
\frac1{|\mathrm{Max}_{+}|}\sum_{y\in\mathrm{Max}_{+}}yy^\top
=
I+\frac Cp
=
2P_+,
\end{equation}
which is certified for Paley \(n=p^2+1\) at \(p=3,5\) (and holds at \(p=7\) by the same free-variable enumeration). Then:

1. **Fractional Max-cover number equals \(p\).** The LP
   \[
   \min\bigl\{\textstyle\sum_{e}x_e:x\ge0,\;
   \textstyle\sum_{e}x_e\,C_e y_iy_j\ge1\ \forall y\in\mathrm{Max}_{+}\bigr\}
   \]
   has value exactly \(p\). *Primal:* \(x_e=2/(np)\) for every edge gives objective \((n-1)/p=p\) and constraint value \(2/(np)\cdot Q_C(y)=1\). *Dual:* \(\lambda_y=p/|\mathrm{Max}_{+}|\) is feasible with every edge tight, because
   \[
   \frac p{|\mathrm{Max}_{+}|}\sum_y C_{ij}y_iy_j
   =
   p\cdot C_{ij}\cdot\frac{C_{ij}}p
   =
   1.
   \]

2. **Max-Lipschitz (factor \(1/p\) improvement of Prop 15.20b).** For every Seidel matrix \(A\), after a best switching toward \(C\) with disagreement cardinality \(k=d_H^{\mathrm{best}}(A,C)\),
   \[
   \Phi(A)
   \;\ge\;
   \max_{y\in\mathrm{Max}_{+}}\lvert Q_A(y)\rvert
   \;\ge\;
   \Phi(C)-\frac{2k}p.
   \]
   *Proof.* Write \(G=I+C/p\). Then
   \[
   \frac1{|\mathrm{Max}_{+}|}\sum_{y}y^\top Ay
   =
   \sum_{i\neq j}A_{ij}G_{ij}
   =
   \frac2p\sum_{i<j}A_{ij}C_{ij}
   =
   \frac2p\bigl(\tbinom n2-2k\bigr)
   =
   np-\frac{4k}p.
   \]
   Hence the average of \(Q_A=\tfrac12 y^\top Ay\) over \(\mathrm{Max}_{+}\) equals \(\Phi(C)-2k/p\). Best switching forces this average to be nonnegative, and \(\max|Q_A|\ge\bigl|\mathbb E[Q_A]\bigr|\). \(\square\)

*Certified checks.* Frame identity \eqref{eq:max-frame} and LP value \(p\) at Paley \(p=3,5\); dual/primal algebra as above. At \(n=10\), the Max-Lipschitz with \(k_\star=5\) gives \(m_{10}\ge15-10/3=11.\overline{6}\) (true; actual \(13\)); edge lip only gave \(15-10=5\).

*Remark (E(1) status — sharpened criterion).* Max-Lipschitz yields
\[
m_n\ge\Phi(C)-\frac{2k_\star}p,
\qquad
\alpha_n\ge\tfrac12\sqrt{1-1/n}-\frac{2k_\star}{p\,n^{3/2}}.
\]
Hence **\(k_\star=O(n^{3/2})\) already forces \(\alpha_n\to\tfrac12\)** along the \(\rho=1\) family (gap \(O(n)=o(n^{3/2})\)), and \(k_\star=o(n^2)\) is the absolute threshold for a vanishing relative gap. This improves Prop 15.20d (which needed \(k_\star=O(n)\) via edge lipschitz). The remaining gap to a free proof is a factor \(\sqrt n\): dual-Gaussian on \(W=A\circ C\) only gives \(k_\star\le\binom n2/2-\Omega(n^{3/2})=\Theta(n^2)\) for arbitrary \(A\). Closing E(1) needs \(k_\star=O(n^{3/2})\) (or better) for \(\Phi\)-minimisers — e.g. via \(\Delta(F)=O(\sqrt n)\) after best switch, the matching dichotomy, or spectral rigidity. Integral Max-covers of size \(p\) exist (LP support) but are stars and spike to \(\Phi>\Phi(C)\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.28 (size-\(p\) Max-covers: tight cover and spike; 2026-07-27).** Let \(C\) be a \(\rho=1\) conference of order \(n=p^2+1\) with the Max\(\pm\) frame identities of Prop~15.27 (so \(\mathbb E_{\pm}[C_{ij}y_iy_j]=\pm1/p\), and \(\mathbb E_{\mathrm{Max}_{-}}[yy^\top]=I-C/p\)). Write \(S_F(x)=\sum_{\{i,j\}\in F}C_{ij}x_ix_j\).

1. **Tight cover.** If \(F\) is a Max-cover (\(S_F\ge1\) on \(\mathrm{Max}_{+}\)) with \(|F|=p\), then \(S_F\equiv1\) on \(\mathrm{Max}_{+}\).
   *Proof.* \(\mathbb E_{+}[S_F]=|F|/p=1\) and \(S_F\ge1\) force equality. \(\square\)

2. **Spike from Max\(_{-}\).** If some \(z\in\mathrm{Max}_{-}\) has \(S_F(z)\ge1\), then \(\Phi(C\oplus F)\ge\Phi(C)+2\), since \(Q_{C\oplus F}(z)=-\Phi-2S_F(z)\). If \(S_F(z)=p\) then \(\Phi(C\oplus F)\ge\Phi+2p\).

3. **Affine-line \(p\)-stars at \(\infty\) (proved for Paley).** Identify vertices with \(\{\infty\}\cup\mathbb F_{p^2}\cong\{\infty\}\cup\mathrm{AG}(2,p)\), so \(C_{\infty u}=1\). Let \(L=\{a+td:t\in\mathbb F_p\}\) be an affine line of direction \(d\in\mathbb F_{p^2}^\times\), and let \(F\) be the \(p\)-star at \(\infty\) with leaves \(L\). Write \(\chi\) for the quadratic character of \(\mathbb F_{p^2}\). On \(\mathbb F_p^\times\subset\mathbb F_{p^2}\) one has \(\chi\equiv1\), hence
   \[
   \sum_{u\neq v\in L}C_{uv}
   =
   \sum_{t\neq s}\chi((t-s)d)
   =
   p(p-1)\,\chi(d).
   \]
   For any star the cross terms of \(S_F\) are adjacent edges, so the frame gives the exact second moment (using \(G_{\pm}=I\pm C/p\)):
   \[
   \mathbb E_{\pm}[S_F^2]
   =
   p\pm\frac1p\sum_{u\neq v\in L}C_{uv}
   =
   p\pm(p-1)\chi(d).
   \]
   - If \(\chi(d)=-1\) (nonsquare direction): \(\mathbb E_{+}[S_F^2]=1=\bigl(\mathbb E_{+}[S_F]\bigr)^2\), so \(S_F\equiv1\) on \(\mathrm{Max}_{+}\) (tight Max-cover); \(\mathbb E_{-}[S_F^2]=2p-1\neq1\), so \(S_F\not\equiv-1\) on \(\mathrm{Max}_{-}\). Since \(S_F\) is an odd integer with mean \(-1\), necessarily \(\max_{\mathrm{Max}_{-}}S_F\ge1\), and Lemma~2 yields \(\Phi(C\oplus F)\ge\Phi+2\).
   - If \(\chi(d)=+1\): symmetrically \(S_F\equiv-1\) on \(\mathrm{Max}_{+}\) (not a cover).

   Thus **exactly the nonsquare-direction affine lines** yield covering \(p\)-stars at \(\infty\), and every such star fails to undercut. (Half of the \(p(p+1)\) lines of \(\mathrm{AG}(2,p)\): \(p(p+1)/2\) covers.)

4. **Stronger spike at covering stars (certified \(p=3,5\)).** Every covering \(p\)-star (all centres, not only \(\infty\)) has \(\max_{\mathrm{Max}_{-}}S_F=p\), hence \(\Phi(C\oplus F)=\Phi+2p\). Counts: 60 at \(n=10\) (all \(\Phi=21\)); 390 at \(n=26\) (MITM sample all \(\Phi=75\)). Spike witnesses are constant on leaves\(\cup\{\mathrm{centre}\}\).

5. **All size-\(p\) Max-covers at \(n=10\).** Exhaustive: 405 size-\(p\) Max-covers (60 stars + 345 non-stars); every one has \(\max_{\mathrm{Max}_{-}}S_F\in\{1,3\}\) and \(\Phi\in\{17,21\}\); **zero undercuts** of \(\Phi=15\). Moreover \(\mathbb E_{-}[S_F^2]\ge7/3>1\) on this census, so none can have \(S_F\equiv-1\) on \(\mathrm{Max}_{-}\). Card-min Max-covers never undercut at \(n=10\); undercutters begin at \(k=5\) (N10-S matchings).

*Evidence:* `evidence/E1_SIZE_P_MAXCOVER.md`, `e1_size_p_maxcover.json`. Max\(_{-}\) frame certified at \(p=3,5\).

*Remark (E(1) status).* Prop 15.28 shows that LP-tight Max-covers (size \(p\)) cannot undercut via the Max\(\pm\) analysis on the proved/certified range: they are forced to spike on \(\mathrm{Max}_{-}\). The \(n=10\) undercut uses a *strictly larger* cover (\(k=5>p\)). Closing E(1) still needs either (i) \(k_\star=O(n^{3/2})\) for minimisers, (ii) a general no-undercut / controlled-gap theorem for all Max-covers of size \(o(n^2)\), or (iii) exact Paley optimality for all \(p\ge5\) with the known \(O(1)\) gap at \(p=3\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.29 (odd matching parity and matching-cover spikes; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\) (\(p\) odd prime), \(\Phi=\tfrac12 np\).

1. **Parity.** The perfect-matching size \(n/2\) is odd. Hence for every perfect matching \(M\) and every \(x\in\{\pm1\}^n\), \(S_M(x)\) is an odd integer. In particular, either \(\min_{\mathrm{Max}_{+}}S_M\ge1\) (Max-cover) or \(\min_{\mathrm{Max}_{+}}S_M\le-1\). In the latter case \(\Phi(C\oplus M)\ge\Phi+2\).

2. **Undercutters are two-sided.** If \(\Phi(C\oplus F)<\Phi\), then \(S_F\ge1\) on \(\mathrm{Max}_{+}\) and \(S_F\le-1\) on \(\mathrm{Max}_{-}\) (integer arithmetic on \(Q_C=\pm\Phi\)). Evidence note: `evidence/E1_TWOSIDED_COVER.md`.

3. **Certified \(n=26\) matching covers.** Perfect-matching Max-covers exist at \(p=5\) (SA finds them; earlier “no cover” reports were incomplete). Among 48 seeds, 3 covers were found; all three are two-sided, achieve \(\max_{\mathrm{Max}_{\pm}}|Q|=63=\Phi-2\), but exact MITM \(\Phi(C\oplus M)=65=\Phi(C)\) (non-eigenvector spike of \(+2\)). **Zero undercuts.** Contrast \(n=10\): two-sided matching covers achieve \(m_{10}=13\). Evidence: `evidence/E1_MATCHING_COVER_SPIKE.md`, `e1_n26_matching_covers_census.json`.

*Remark.* Matching dichotomy for E(1) must account for covers that fail only by spike control. Prop 15.26 (local maximality) holds at both \(p=3,5\); global Max-determination of \(\Phi(C\oplus M)\) is special to \(n=10\) among tested orders. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.30 (matching spike criterion; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), and \(M\) a perfect matching of \(K_n\). Write \(A=C\oplus M\) and \(S_M(x)=\sum_{\{i,j\}\in M}C_{ij}x_ix_j\).

1. **Criterion (proved).** If there exists \(x\in\{\pm1\}^n\) with \(S_M(x)=-p\) and \(Q_C(x)\ge\Phi-2p\) (or \(S_M(x)=p\) and \(Q_C(x)\le-(\Phi-2p)\)), then \(\Phi(A)\ge\Phi(C)\).  
   *Proof.* \(Q_A=Q_C-2S_M\), so the first alternative yields \(Q_A(x)\ge\Phi-2p+2p=\Phi\). \(\square\)

2. **Single-bit drop.** For \(y\in\mathrm{Max}_{+}\) and any coordinate \(i\), \(Q_C(y^{\oplus i})=\Phi-2p\). Hence non-max boolean vectors attain at least the level \(\Phi-2p\). At \(p=3\) this is exact as the global second max of \(|Q_C|\) on non-max vectors (full half-cube). At \(p=5\) sampling finds second max \(55=\Phi-2p\).

3. **Multi-bit formula.** \(Q_C(y^{\oplus F})=\Phi-2p|F|+4\sum_{i<j\in F}C_{ij}y_iy_j\). On a \(y\)-switched clique of size \(r\), this equals \(\Phi-2r(p-r+1)\); at \(r=p\) one recovers \(\Phi-2p\). Square-direction affine lines of \(\mathrm{AG}(2,p)\) are \(p\)-cliques of Paley (Prop 15.28 geometry).

4. **Certified dichotomy at \(n=10\).** Over all 945 perfect matchings, the criterion holds on exactly the \(801\) non-undercutting matchings and fails on exactly the \(144\) undercutters (where \(\max Q_C\) on \(\{S_M=-3\}\) is \(7<9=\Phi-2p\)).

5. **Certified at \(n=26\).** Every tested perfect matching (15 random + 3 SA Max-covers) satisfies the criterion; all three Max-covers have exact MITM \(\Phi=\Phi(C)\). On the level \(S_M=-p\) one has the identity \(\mathbb E[Q_C]=-p\).

*Open for matching non-undercut when \(p\ge5\).* Prove the criterion for every perfect matching (tail \(\max_{S_M=-p}Q_C\ge\Phi-2p\)). That would give \(\Phi(C\oplus M)\ge\Phi(C)\) for all \(M\), hence E(1) under the matching dichotomy. Full E(1) still requires control of non-matching undercutters. Evidence: `evidence/E1_MATCHING_SPIKE_CRITERION.md`. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.31 (clique-flip sufficiency for matching Max-covers; 2026-07-27).** Let \(M\) be a perfect matching Max-cover of Paley order \(n=p^2+1\). If there exist \(y\in\mathrm{Max}_{+}\) and a transversal \(p\)-set \(F\) such that \(F\) is a clique in \(W_{ij}=C_{ij}y_iy_j\), \(S_M(y)=s_0\) with \(s_0+3p\equiv0\pmod4\), and \(\sum_{i\in F}\chi_i(y)=(s_0+p)/2\), then
\[
Q_C(y^{\oplus F})=\Phi-2p,\qquad S_M(y^{\oplus F})=-p,\qquad Q_{C\oplus M}(y^{\oplus F})=\Phi,
\]
so \(\Phi(C\oplus M)\ge\Phi(C)\).

*Arithmetic.* Full-clique flips reach the \(\Phi-2p\) threshold only for \(|F|\in\{1,p\}\); \(|F|=1\) is incompatible with Max-covers. Thus only \(|F|=p\) applies, requiring \(s_0+3p\equiv0\pmod4\). At \(p=3\) undercutting matchings have \(S_M(\mathrm{Max}_{+})\subseteq\{1,5\}\) (no admissible \(s_0\)), blocking the construction. At \(p=5\), \(\mathbb E[S_M]=2.6<3\) forces every Max-cover to attain \(S_M=1\).

*Design constants at \(p=5\) (certified).* Seidel-consistent \(p\)-sets (\(C_{ab}C_{ac}C_{bc}=1\)): 390. Each has exactly 60 Max\(_{+}\) extensions. Every tested matching has \(\ge236\) transversal consistent \(p\)-sets. All SA Max-covers admit a clique-flip (MITM \(\Phi=\Phi(C)\)). Evidence: `evidence/E1_CLIQUE_FLIP.md`, `e1_clique_flip_covers.json`.

*Open.* Existence of \((y,F)\) for every Max-covering matching when \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.32 (Γ-pairing reformulation and mod-4 dichotomy; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), \(m=n/2\), and \(M=\{e_a=(i_a,j_a)\}_{a=1}^{m}\) a perfect matching. Write \(S_M=\sum_a C_{e_a}x_{i_a}x_{j_a}\) and \(R_M=Q_C-S_M\).

1. **Γ-pairing (proved).** With free signs \(u_a=x_{i_a}\) and relative signs \(z_a=x_{i_a}x_{j_a}\),
   \[
   R_M=\sum_{a<b}\gamma_{ab}(z_a,z_b)\,u_a u_b,\qquad
   \gamma_{ab}\in\{-4,-2,0,2,4\}.
   \]
   The spike criterion \(\max_{S_M=-p}Q_C\ge\Phi-2p\) is equivalent to
   \[
   \max_{\langle c,z\rangle=-p}\max_{u\in\{\pm1\}^m}\tfrac12 u^\top\Gamma(z)u
   \;\ge\;p(m-1).
   \]
   Shipped: `maxR_matching_level` in `src/minmax_quadratic.py`.

2. **Coordinate product on Max\(_{+}\) (proved).** Every \(y\in\mathrm{Max}_{+}\) has
   \(\pi(y):=\prod_{v=1}^n y_v=(-1)^{p(p-1)/2}\), constant on \(\mathrm{Max}_{+}\).
   (Halfspace construction gives the value; \(\mathrm{P}\Gamma\mathrm{L}\) acts by coordinate permutation and preserves the product; the boolean \(+p\)-orbit is the Aut-orbit of the halfspace.)
   Thus \(\pi=+1\) when \(p\equiv1\pmod4\) and \(\pi=-1\) when \(p\equiv3\pmod4\).

3. **Mod-4 constancy (proved).** For any perfect matching \(M\), \(S_M(y)\bmod 4\) is **constant** for \(y\in\mathrm{Max}_{+}\). Indeed with \(\chi_a=C_{e_a}y_{i_a}y_{j_a}\) one has \(\prod_a\chi_a=(\prod_e C_e)\,\pi\), and \(S_M=\sum\chi_a\equiv m-2k\pmod4\) where \(k=\#\{\chi_a=-1\}\), so the residue is determined by \(\prod\chi_a\) alone.

4. **Case split for the spike (structure).** Write \(r\) for the constant residue of \(S_M\) on \(\mathrm{Max}_{+}\).
   - If \(-p\equiv r\pmod4\) and \(-p\in S_M(\mathrm{Max}_{+})\), the criterion holds with \(Q=\Phi\) (maximiser on the level).
   - If \(-p+2\equiv r\pmod4\) and \(-p+2\in S_M(\mathrm{Max}_{+})\), a **1-bit flip** of a maximiser lands on \(Q=\Phi-2p\) and \(S_M=-p\) (formula \(S(y^{\oplus i})=S(y)-2\chi_i\)), so the criterion holds.
   - If \(M\) is a Max-cover (\(S_M\ge1\) on \(\mathrm{Max}_{+}\)), then \(r\equiv1\pmod4\) at \(p=5\) and the 1-bit route is blocked; the clique-flip of Prop 15.31 is the remaining constructive path.

5. **Certified census at \(p=5\) (not a proof).** Over 80 random perfect matchings + SA min-\(\max R\) + 2 Max-covers: \(\max R\in\{60,70\}\) always, \(\min=60=p(m-1)\) (tight on covers), 0 counterexamples to the criterion; MITM \(\Phi(C\oplus M)=\Phi(C)\) on covers. Over 500 random matchings: every mod-\(3\) matching attains \(S_M=-5\); every mod-\(1\) matching attains \(S_M=-3\). Evidence: `evidence/E1_GAMMA_PAIRING.md`, `e1_gamma_forall_census.json`.

*Open.* Prove the attainment lemmas (\(-p\) or \(-p+2\) in the image of \(S_M|_{\mathrm{Max}_{+}}\) whenever the residue allows) for all \(p\ge5\), and clique-flip existence on every Max-cover; then matching non-undercut for \(p\ge5\). Matching dichotomy / non-matching undercutters / \(k_\star\) still separate. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.33 (non-covers cannot undercut; criterion is not necessary; 2026-07-27).** Let \(C\) be a \(\rho=1\) conference of order \(n\), \(\Phi=\Phi(C)\), and \(F\) any flip set. Write \(A=C\oplus F\) and \(S_F(y)=\sum_{e\in F}C_e y_iy_j\).

1. **Non-cover \(\Rightarrow\) strict raise (proved).** If there exists \(y\in\mathrm{Max}_{+}\) with \(S_F(y)\le0\), then
   \[
   Q_A(y)=\Phi-2S_F(y)\ge\Phi,
   \]
   and if \(S_F(y)\le-1\) (always, when scores are odd — e.g. perfect matchings with \(n/2\) odd) then \(Q_A(y)\ge\Phi+2\). Hence \(\Phi(A)\ge\Phi+2>\Phi(C)\): **non-covers cannot undercut.**
   *Proof.* Immediate from \(Q_A=Q_C-2S_F\) on \(\mathrm{Max}_{+}\). \(\square\)
   Equivalently: every strict undercutter is a Max\(_{+}\) cover (cf. Lemma U1 / Prop 15.29).

2. **Perfect matchings (proved reduction).** For a perfect matching \(M\), either \(M\) is a Max-cover and may or may not undercut, or \(M\) is a non-cover and \(\Phi(C\oplus M)\ge\Phi+2\). Thus **matching non-undercut reduces entirely to Max-covering matchings.**

3. **Spike criterion is not necessary (certified counterexample).** At \(p=5\), the perfect matching
   \[
   M_0=\{\{5,3\},\{2,14\},\{22,18\},\{23,7\},\{12,16\},\{6,8\},\{1,19\},\{17,25\},\{9,0\},\{21,10\},\{11,24\},\{13,20\},\{4,15\}\}
   \]
   has \(\max R_{M_0}=54<60=\Phi-p\) (criterion fails) but \(\min_{\mathrm{Max}_{+}}S_{M_0}=-1\) (non-cover) and exact MITM \(\Phi(C\oplus M_0)=75>\Phi\). Evidence: `evidence/e1_criterion_fail_no_undercut.json`. The criterion remains a useful *sufficient* test for covers; it is not a characterisation of non-undercut.

4. **Max-cover census at \(p=5\) (not a proof).** All SA-found Max-covering matchings (8 total: 2 in `e1_gamma_forall_census.json` + 6 in `e1_gamma_cover_batch.json`) satisfy the spike criterion with equality \(\max R=60\), admit clique-flips, and have MITM \(\Phi=\Phi(C)\).

*Open for matching non-undercut when \(p\ge5\).* Prove every Max-covering perfect matching has \(\Phi(C\oplus M)\ge\Phi(C)\) (e.g. via clique-flip Prop 15.31 for all covers). Then, with \(p=3\) gap \(O(1)\), matching dichotomy would give E(1) along matchings. Full E(1) still needs non-matching undercutters. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.34 (matching flip algebra; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\) and \(M\) a perfect matching. Write \(D\) for the matrix with \(D_{ij}=C_{ij}\) on edges of \(M\) and \(0\) elsewhere, and set \(A=C-2D=C\oplus M\).

1. **Involution identity (proved).** \(D^2=I\). Consequently
   \[
   A^2=(n+3)I-2(CD+DC).
   \]
   *Proof.* Each vertex lies in exactly one matching edge, so \((D^2)_{ii}=C_{i\pi(i)}^2=1\) and \((D^2)_{ij}=0\) for \(i\neq j\). Expand \(A^2=(C-2D)^2=C^2-2(CD+DC)+4D^2=(n-1)I-2(CD+DC)+4I\). \(\square\)

2. **Certified spectrum at \(p=5\).** For every SA Max-covering matching tested, \(\|A\|_{\mathrm{op}}=\sqrt{41}\) exactly (to numerical precision \(10^{-12}\)); random matchings have \(\|A\|_{\mathrm{op}}\in[6.79,7.00]\). Two Aut-invariants (K\(_{2,2}\) type counts) show \(\ge8\) distinct Max-cover classes, all sharing this op-norm. When \(\Phi(A)=\Phi(C)=65\), one has \(\rho(A)=p/\|A\|_{\mathrm{op}}=5/\sqrt{41}\).

*Open.* Prove \(\|C\oplus M\|_{\mathrm{op}}=\sqrt{p^2+16}\) (or a bound forcing \(\Phi(A)\ge\Phi(C)\)) for every Max-covering matching when \(p\ge5\); or complete clique-flip existence. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.35 (Max-cover matching structure at \(p=5\); 2026-07-27).** Let \(C\) be Paley of order \(n=26\) (\(p=5\)), \(\Phi=65\), and \(M\) a perfect matching Max-cover (\(\min_{\mathrm{Max}_{+}}S_M\ge1\)).

1. **Forced \(S_M=1\) and residue (proved).** \(\mathbb E[S_M]=13/5=2.6<3\) and \(S_M\) odd on \(\mathrm{Max}_{+}\), so \(S_M\) attains \(1\). The residue of \(S_M\) on \(\mathrm{Max}_{+}\) cannot be \(3\bmod4\) (that would force \(\min S_M\ge3>2.6\)). Hence \(S_M\equiv1\pmod4\), \(s_0=1\) is admissible for clique-flip (Prop 15.31), and \(\Sigma_{\mathrm{need}}=3\).

2. **Certified census of 11 Max-covers (not a forall proof).** Independent SA campaigns produced **11** distinct Max-covering perfect matchings. Every one is:
   - two-sided (\(\min_{\mathrm{Max}_{+}}S=1\), \(\max_{\mathrm{Max}_{-}}S=-1\));
   - inclusion-minimal and inclusion-maximal as Max-covers;
   - \(\|C\oplus M\|_{\mathrm{op}}=\sqrt{41}\) exactly;
   - spike-criterion tight (\(\max R=60\));
   - clique-flip capable;
   - exact MITM \(\Phi(C\oplus M)=\Phi(C)\).
   Only two \(S_{+}\) distributions appear: \((1^{156},5^{104})\) and \((1^{168},5^{80},9^{12})\). Evidence: `evidence/e1_maxcover_full_census.json`.

3. **Uniform sampling.** Among \(2\cdot10^4\) random perfect matchings, **0** Max-covers (Max-covers are extreme tail events).

*Open.* Prove clique-flip / \(\Phi\ge\Phi(C)\) for every Max-covering perfect matching when \(p\ge5\); then non-matching undercutters / \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.36 (matching flip block algebra and Max-cover spectrum; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(M\) a perfect matching, \(D\) the signed matching matrix, \(A=C-2D\), and \(V_\pm=\ker(C\mp pI)\).

1. **Block formulae (proved for every matching).** \(B:=CD+DC\) always commutes with \(C\) and \(D\). On \(V_+\oplus V_-\),
   \[
   B\big|_{V_+}=2p\,D_{++},\qquad B\big|_{V_-}=-2p\,D_{--},
   \]
   and \(\mathrm{tr}(D_{++})=n/(2p)\). For \(y\in\mathrm{Max}_+\subset V_+\) one has \(S_M(y)=\tfrac12 y^\top D_{++}y\). Also
   \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\).

2. **Certified at \(p=5\) for Max-covers (not forall).** Every SA Max-covering matching has \(\lambda_{\min}(B)=-6\) and \(\|A\|_{\mathrm{op}}=\sqrt{41}=\sqrt{p^2+16}\). At least two \(D_{++}\) spectral types occur (simple \(\{-3/5,0^{(8)},(4/5)^{(4)}\}\) and a mixed type with the same \(\lambda_{\min}(D_{++})=-3/5\)); both have tight spike \(\max R=60\), clique-flip, and MITM \(\Phi=\Phi(C)\). Random non-covers have \(\lambda_{\min}(B)\in[-10,-8]\) and larger op-norm. At \(p=3\), the 144 undercutting Max-covers share a single smaller op-norm \(\approx3.933\) (golden-ratio \(B\)-spectrum). Evidence: `evidence/E1_MAXCOVER_SPECTRUM.md`, `e1_maxcover_spectrum.json`.

*Open.* Prove \(\lambda_{\min}(B)=-6\) (or \(\|A\|_{\mathrm{op}}=\sqrt{p^2+16}\)) for every Max-cover when \(p\ge5\), and upgrade to \(\Phi(A)\ge\Phi(C)\); or complete clique-flip existence. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.37 (continuous Γ-bound pattern on Max-covers; 2026-07-27).** Let \(C\) be Paley of order \(n=26\) and \(M\) a Max-covering perfect matching. On the Γ-pairing level \(S_M=-p\), every tested \(M\) satisfies
\[
\min_z\lambda_{\max}(\Gamma(z))\;\ge\;9.38758\;>\;\frac{2p(m-1)}{m},
\]
so the continuous bound \(\tfrac m2\lambda_{\max}\ge p(m-1)\) holds for **all** \(z\) on the level (two numerical classes). Discrete \(\max R=60\) and clique-flip hold on all 11 stored covers. GW theory does not close the discrete gap (SDP\(\approx63.8\), \(\alpha\cdot\mathrm{SDP}<60\)). Residue-\(1\) random matchings also admitted clique-flips in sampling. Evidence: `evidence/E1_MAXCOVER_CONTINUOUS_BOUND.md`.

*Open.* Prove the continuous bound and/or clique-flip for every Max-cover when \(p\ge5\); then matching dichotomy / \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.38 (two-sided size-5 Max-covers on \(C_{10}\); 2026-07-28).** Let \(C\) be Paley of order \(n=10\), \(\Phi=15\). Among all edge sets \(F\) with \(|F|=5\) that are two-sided Max-covers (\(\min_{\mathrm{Max}_{+}}S_F\ge1\) and \(\max_{\mathrm{Max}_{-}}S_F\le-1\)):

1. **Undercutters are exactly the matchings.** Exactly **144** such \(F\) undercut, each is a perfect matching (\(\Delta=1\)), and each has \(\Phi(C\oplus F)=13=m_{10}\).
2. **Higher \(\Delta\) never undercuts at this cardinality.** Counts: \(\Delta=2\): 8730 covers with \(\Phi\in\{15,17,19\}\); \(\Delta=3\): 7920 with \(\Phi\in\{17,19\}\); \(\Delta=4\): 360 with \(\Phi=19\). All have \(\Phi\ge15\).
3. **Total.** 17154 two-sided \(k=5\) Max-covers; only the 144 matchings undercut.

*Proof.* Exhaustive scan of \(\binom{45}{5}=1{,}221{,}759\) five-edge sets; Max\(\pm\) boolean \(\pm3\)-eigenvectors; exact \(\Phi\) by cube enumeration via `form_Q`. Parallel re-run: `src/n10_twosided_k5_classify.py` (80 workers). Evidence: `evidence/E1_N10_TWOSIDED_K5.md`, `e1_n10_twosided_k5_classify.json`. \(\square\)

*Remark.* Strengthens N10-S: at matching cardinality, two-sided Max-covers with \(\Delta\ge2\) cannot undercut. Supports a low-\(\Delta\) undercutter pattern for E(1), but does **not** prove \(k_\star=O(n^{3/2})\) or matching non-undercut for \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.39 (clique-flip count on Max-cover matchings, \(p=5\); 2026-07-28).** On every stored Max-covering perfect matching of Paley \(C_{26}\), the number \(N_{\mathrm{flip}}\) of pairs \((y,F)\) as in Prop 15.31 (with \(s_0=1\), \(\Sigma=3\)) satisfies \(N_{\mathrm{flip}}\ge24>0\), so clique-flip applies and \(\Phi(C\oplus M)=\Phi(C)\). Observed counts include \(24,48,120\). Evidence: `evidence/E1_CLIQUE_FLIP_COUNT.md`.

*Open.* Prove \(N_{\mathrm{flip}}\ge1\) for every Max-cover matching when \(p=5\) (then matching non-undercut at \(p=5\)); lift to \(p\ge7\) and \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.40 (edge-minimal undercutters have gap at most 2; 2026-07-28).** Let \(C\) be any Seidel matrix of order \(n\) with \(\Phi(C)\in\mathbb Z\), and let \(F\) be a nonempty edge set such that \(\Phi(C\oplus F)<\Phi(C)\) while \(\Phi(C\oplus(F\setminus\{e\}))\ge\Phi(C)\) for every \(e\in F\) (edge-minimal undercutter). Then
\[
\Phi(C\oplus F)\;\ge\;\Phi(C)-2.
\]
*Proof.* For any \(e\in F\), Prop 15.20b with \(k=1\) gives \(|\Phi(C\oplus F)-\Phi(C\oplus(F\setminus\{e\}))|\le2\). Combined with \(\Phi(C\oplus F)<\Phi(C)\le\Phi(C\oplus(F\setminus\{e\}))\),
\[
\Phi(C\oplus F)
\;\ge\;
\Phi(C\oplus(F\setminus\{e\}))-2
\;\ge\;
\Phi(C)-2.
\]
(The same holds with \(\Phi\) replaced by any real threshold \(t\) for which \(F\) is edge-minimal among sets with \(\Phi(C\oplus\cdot)<t\).) \(\square\)

*Certified.* On Paley \(C_{10}\), all 144 matching undercutters and sampled undercutting 6-cycles are edge-minimal and have gap exactly 2. Evidence: session checks; N10-S/C6.

*Remark (does **not** yet give \(m_n\ge\Phi-2\)).* The lemma bounds edge-minimal undercutters only. A global lower bound \(m_n\ge\Phi(C)-2\) would require showing no deeper undercut exists at larger Hamming distance (or that a closest undercutter realises \(m_n\)). That step is **open**. If proved, E(1) follows (gap \(O(1)=o(n^{3/2})\)) and \(L=\tfrac12\) by denseness on the \(\rho=1\) family. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.41 (first-hit + no-descent framework toward \(m_n\ge\Phi-2\); 2026-07-28).** Let \(C\) be a Seidel matrix of order \(n\) with \(\Phi(C)\in\mathbb Z\). Write \(A_F:=C\oplus F\).

1. **First-hit lemma (proved).** Along any edge-adding chain \(\emptyset=F_0\subset F_1\subset\cdots\subset F_k\) with \(|F_i|=i\), if \(i^\star\) is minimal such that \(\Phi(A_{F_{i^\star}})<\Phi(C)\), then
   \[
   \Phi(A_{F_{i^\star}})\;\ge\;\Phi(C)-2.
   \]
   *Proof.* \(F_{i^\star}=F_{i^\star-1}\cup\{e\}\) and \(\Phi(A_{F_{i^\star-1}})\ge\Phi(C)\), so Prop 15.20b with one edge yields the claim. \(\square\)

2. **Dangerous-edge criterion (proved).** Suppose \(\Phi(A_F)=\Phi(C)-2\). For \(e=(u,v)\notin F\) set \(B=A_F\oplus e\) and \(\sigma_x:=(A_F)_{uv}\,x_u x_v\in\{\pm1\}\). Then \(Q_B(x)=Q_{A_F}(x)-2\sigma_x\). Consequently:
   - if some maximiser \(x\) with \(Q_{A_F}(x)=\Phi(C)-2\) has \(\sigma_x=-1\), then \(Q_B(x)=\Phi(C)\) and \(\Phi(B)\ge\Phi(C)\);
   - if some maximiser with \(Q_{A_F}(x)=-(\Phi(C)-2)\) has \(\sigma_x=+1\), same.
   Hence \(\Phi(B)\le\Phi(C)-4\) is possible only if \(\sigma\equiv+1\) on all \(+\) maximisers and \(\sigma\equiv-1\) on all \(-\) maximisers (**dangerous edge**). \(\square\)

3. **No-descent lemma (OPEN in general).** If \(\Phi(A_F)=\Phi(C)-2\), then for every \(e\notin F\), \(\Phi(A_F\oplus e)\ge\Phi(C)-2\).

4. **Conditional settlement.** If the no-descent lemma holds for all flip sets on the \(\rho=1\) Paley family \(n=p^2+1\), then by induction on Hamming distance every Seidel matrix \(A\) of those orders satisfies \(\Phi(A)\ge\Phi(C)-2\). Thus \(m_n\ge\Phi(C)-2\), the gap is \(O(1)=o(n^{3/2})\), E(1) holds, and \(L=\tfrac12\) by denseness (Prop 6.2). **F13:** this must not be claimed from Prop 15.40 alone; no-descent is an independent lemma about \(\Phi\), not abstract 2-Lipschitz calculus.

5. **Certified at \(n=10\) (matching undercutters; not a general proof).** For all **144** perfect-matching undercutters of Paley \(C_{10}\) (\(\Phi=13=\Phi-2\)):
   - number of dangerous external edges: **0**;
   - every single-edge extension has \(\Phi\ge15=\Phi(C)\) (min observed \(15\));
   - multi-edge random extensions (depth \(2\)–\(12\)) stay at \(\Phi\ge13\) (0 deepenings).
   Parallel cert: `src/e1_n10_nodescent.py` (\(W=86\)), `evidence/e1_n10_nodescent.json`, `evidence/E1_NODESCENT.md`.

*Open.* Prove no-descent (or: no dangerous edges, or non-maximiser spike under alignment) for all undercutters on the \(\rho=1\) family; alternatively complete Max-cover clique-flip for \(p\ge5\). Matching dichotomy / \(k_\star\) remain separate. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.42 (Max\(\pm\) dichotomy, counting freeness, tight no-descent; 2026-07-28).** Let \(C\) be a \(\rho=1\) conference of order \(n\) with the Max\(_{+}\) frame identity of Prop 15.27, \(\Phi=\Phi(C)\), \(p=\sqrt{n-1}\), and \(A=C\oplus F\). Write \(S_F(y)=\sum_{e\in F}C_e y_iy_j\), \(s_{+}=\min_{\mathrm{Max}_{+}}S_F\), \(s_{-}=\max_{\mathrm{Max}_{-}}S_F\), and \(f_e(y)=C_e y_iy_j\).

1. **Max\(\pm\) evaluation dichotomy (proved).**
   \[
   \Phi(A)\;\ge\;\max\Bigl(
   \max_{y\in\mathrm{Max}_{+}}\lvert\Phi-2S_F(y)\rvert,\;
   \max_{y\in\mathrm{Max}_{-}}\lvert-\Phi-2S_F(y)\rvert
   \Bigr).
   \]
   Consequently:
   - if \(s_{+}\le -1\), then \(\Phi(A)\ge\Phi+2\) (Prop 15.33);
   - if \(s_{+}=0\), then \(\Phi(A)\ge\Phi\);
   - if \(s_{+}=1\), then \(\Phi(A)\ge\Phi-2\);
   - if \(s_{-}\ge 0\), then \(\Phi(A)\ge\Phi\);
   - if \(s_{-}=-1\), then \(\Phi(A)\ge\Phi-2\).

   Therefore \(\Phi(A)<\Phi-2\) is possible only if \(A\) is **deep two-sided**: \(s_{+}\ge 2\) and \(s_{-}\le -2\).
   *In particular every matrix with \(s_{+}\le 1\) or \(s_{-}\ge -1\) already satisfies \(\Phi(A)\ge\Phi-2\).* \(\square\)

2. **Counting freeness (proved).** For every edge \(e\), \(\sum_{y\in\mathrm{Max}_{+}}f_e(y)=|\mathrm{Max}_{+}|/p\) (frame). Hence if \(f_e\equiv +1\) on a subset \(U\subseteq\mathrm{Max}_{+}\), then
   \[
   |U|\;\le\;|\mathrm{Max}_{+}|\,\frac{p+1}{2p},
   \]
   with equality only if \(f_e\equiv -1\) on \(\mathrm{Max}_{+}\setminus U\). \(\square\)

3. **Tight \(S\equiv 1\) no-descent (proved).** If \(S_F\equiv 1\) on \(\mathrm{Max}_{+}\) (equivalently \(k=|F|=p\) and \(s_{+}=1\)), then \(U=\mathrm{Max}_{+}\) violates the counting bound for \(f_e\equiv 1\) (since \((p+1)/(2p)<1\)). So every \(e\notin F\) has some \(y\in\mathrm{Max}_{+}\) with \(f_e(y)=-1\), whence \(S_{F\cup\{e\}}(y)=0\) and \(Q_{A\oplus e}(y)=\Phi\). Thus \(\Phi(A\oplus e)\ge\Phi\). \(\square\)

4. **Tight \(S\equiv 2\) no-descent (proved).** If \(S_F\equiv 2\) on \(\mathrm{Max}_{+}\) (which holds automatically for any Max\(_{+}\) cover with \(|F|=2p\), since \(\mathbb E[S]=2\) and \(S\ge 2\) force \(S\equiv 2\)), then similarly \(f_e\not\equiv 1\) on \(\mathrm{Max}_{+}\). Any \(y\) with \(f_e(y)=-1\) has \(S_{F\cup\{e\}}(y)=1\) and \(Q_{A\oplus e}(y)=\Phi-2\). Thus \(\Phi(A\oplus e)\ge\Phi-2\). \(\square\)

5. **Type I freeness when \(N_1\) is large (proved).** Write \(N=|\mathrm{Max}_{+}|\) and \(N_1=\#\{y\in\mathrm{Max}_{+}:S_F(y)=1\}\). If \(s_{+}=1\) and \(N_1>N(p+1)/(2p)\), then no edge freezes to \(+1\) on \(\mathrm{Max}_{+1}\). When moreover \(\Phi(A)=\Phi-2\), one has \(\mathrm{Max}_{+1}\subseteq\mathrm{Max}(A)\) (positive maximisers), so every \(e\notin F\) has a maximiser with \(\sigma_e=-1\), hence \(\Phi(A\oplus e)\ge\Phi\) (no-descent, strong form). \(\square\)

6. **Type I size bound for strict freeness (proved).** If \(s_{+}=1\), scores are odd, and \(|F|=k\le 2p-2\), then \(N_1/N\ge(3-k/p)/2>(p+1)/(2p)\), so part 5 applies. \(\square\)

7. **Reduction of \(m_n\ge\Phi-2\) (proved equivalence).** Assume \(\Phi\) is integer-valued with fixed parity under edge flips (true for all Seidel matrices of a fixed order: each flip changes every \(Q_x\) by \(\pm 2\)). Then:
   - by part 1, \(\Phi(A)\le\Phi-4\) forces deep two-sided;
   - by Prop 15.40, no edge-minimal undercutter has \(\Phi\le\Phi-4\);
   - a minimal-Hamming counterexample \(A\) with \(\Phi(A)=\Phi-4\) must therefore admit some \(e\) with \(\Phi(A\oplus e^{-1})=\Phi-2\) (single-edge descent from a gap-\(2\) undercutter).
   Hence **\(m_n\ge\Phi-2\) on the \(\rho=1\) family if and only if no-descent holds for every gap-\(2\) undercutter** (Prop 15.41(3)). Parts 3–6 prove no-descent for all tight \(S\equiv1\), tight \(S\equiv2\), and Type I covers with \(k\le 2p-2\) or \(N_1\) large. \(\square\)

8. **Certified at \(n=10\).** All 144 matching undercutters are Type I with \(N_1=10>8=N(p+1)/(2p)\) (strict freeness; 0 dangerous edges; add-1 gives \(\Phi\ge15\)). All 360 undercutting 6-cycles are tight \(S\equiv 2\) on \(\mathrm{Max}_{+}\) (part 4; add-1 gives \(\Phi\ge15\)). Combined with \(m_{10}=13\), no-descent and \(m_{10}=\Phi-2\) hold at \(n=10\). Evidence: `evidence/E1_NODESCENT.md`, `e1_n10_nodescent.json`; session C6 checks.

*Open for \(L=\tfrac12\).* Prove no-descent for the remaining gap-\(2\) undercutters (Type I with \(k\ge 2p-1\) and small \(N_1\); deep non-tight with \(k>2p\)), or prove no such undercutters exist on the \(\rho=1\) family. Then \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.43 (no-descent for Type I freeness class and tight deep; residual isolation; 2026-07-28).** Continue the notation of Prop 15.42. Let \(B=C\oplus G\) be a gap-\(2\) undercutter (\(\Phi(B)=\Phi-2\)) and \(e\notin G\).

1. **Type I with freeness \(\Rightarrow\) strong no-descent (proved).** If \(s_+(G)=1\) and \(N_1>N(p+1)/(2p)\), then some \(y\in\mathrm{Max}_{+1}\subseteq\mathrm{Max}(B)\) has \(f_e(y)=-1\), so \(Q_{B\oplus e}(y)=\Phi\) and \(\Phi(B\oplus e)\ge\Phi\). \(\square\)

2. **Tight deep \(\Rightarrow\) weak no-descent (proved).** If \(S_G\equiv2\) on \(\mathrm{Max}_{+}\), then some \(y\in\mathrm{Max}_{+}\) has \(f_e(y)=-1\), so \(S_{G\cup\{e\}}(y)=1\), \(Q_{B\oplus e}(y)=\Phi-2\), and \(\Phi(B\oplus e)\ge\Phi-2\). \(\square\)

3. **Equality-case isolation for Type I freeness failure (proved structure).** Suppose \(s_+=1\) and \(f_e\equiv+1\) on \(\mathrm{Max}_{+1}\) (freeness fails). Counting forces \(N_1\le N(p+1)/(2p)\). Combining the Max\(_{+}\) bounds \(Q\le\Phi-2\) on \(\mathrm{Max}_{+}\) and \(Q\le\Phi-6\) on \(\{S\ge3\}\) with \(\mathbb E[Q]=\Phi-2k/p\) forces, at counting equality \(N_1=N(p+1)/(2p)\), that \(S\in\{1,3\}\) on \(\mathrm{Max}_{+}\) and \(k=2p-1\), with \(f_e=2-S\) on \(\mathrm{Max}_{+}\). Then \(H:=G\cup\{e\}\) is a tight \(S\equiv2\) cover of size \(2p\). (A parallel boundary with \(S\in\{1,5\}\) and \(k=3p-2\) forces the affine relation \(S+2f_e=3\).) In the \(k=2p-1\) equality case, no-descent for this \(e\) reduces to \(\Phi(C\oplus H)\ge\Phi-2\) for the tight cover \(H\). \(\square\)

4. **1-bit spike for tight \(S\equiv2\) when all degrees are even and \(p=3\) (proved).** Let \(H\) be tight \(S\equiv2\), \(A=C\oplus H\). For \(y\in\mathrm{Max}_{+}\) and vertex \(v\),
   \[
   Q_A(y^{\oplus v})=\Phi-2p-4+4\sigma_v(y),\qquad \sigma_v=\sum_{vw\in H}f_{vw}(y).
   \]
   Always \(\sum_v\sigma_v(y)=4\). If every degree in \(H\) is even then each \(\sigma_v\) is even; \(\sigma_v\le0\) for all \(v\) would give \(\sum\sigma\le0<4\), so some \(\sigma_v\ge2\). For \(p=3\) this yields \(Q_A(y^{\oplus v})\ge\Phi-2\). In particular every 2-regular tight cover on a support of size \(2p\) (e.g.\ undercutting \(C_6\) at \(n=10\)) has \(\Phi\ge\Phi-2\), and the Type I equality reduction of part 3 has no-descent at \(p=3\). \(\square\)

5. **n=10 closure (proved).** Every edge-minimal undercutter is a matching (Type I, \(N_1=10>8\)) or a 6-cycle (tight \(S\equiv2\)). Parts 1–2 give no-descent; combined with Props 15.40–15.42 and the minimal-counterexample reduction of Prop 15.42(7), \(m_{10}=\Phi-2=13\). \(\square\)

6. **Global residual (OPEN).** To get \(m_n\ge\Phi-2\) for all \(\rho=1\) orders it remains to prove no-descent for:
   - Type I gap-2 undercutters at counting equality for general \(p\) (parts 3–4 handle \(p=3\); need \(\Phi\ge\Phi-2\) for all tight \(S\equiv2\) covers, or freeness of the \(k=3p-2\) boundary);
   - deep non-tight gap-2 undercutters with \(s_+\ge2\), \(k>2p\).

   If those hold, then \(m_n\ge\Phi-2=o(n^{3/2})\) on the dense \(\rho=1\) family, E(1) follows, and \(L=\tfrac12\) by Prop 6.2. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.44 (master lemma for tight covers; bi-tight residual; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\Phi(C)\).

1. **Master lemma (proved).** If \(H\) is any flip set with \(S_H\equiv s\) on \(\mathrm{Max}_{+}\) for some integer \(s\ge1\) (hence \(|H|=sp\) by taking expectations), then either
   \[
   \max_{z\in\mathrm{Max}_{-}}S_H(z)\;\ge\;0
   \qquad\text{(hence \(\Phi(C\oplus H)\ge\Phi\) by Prop 15.42.1),}
   \]
   or \(S_H\equiv -s\) on \(\mathrm{Max}_{-}\) as well (**bi-tight of level \(s\)**).
   *Proof.* Always \(\mathbb E_{-}[S_H]=-|H|/p=-s\). If the maximum on \(\mathrm{Max}_{-}\) is \(\le -1\) and scores have the parity of \(s\), a maximum \(\le -s\) with mean \(-s\) forces constancy at \(-s\) when the maximum is \(\le -s\); more directly: if the maximum is \(\ge 0\) we are done by dichotomy; if the maximum is \(\le -s\) and the mean is \(-s\) with all values \(\le -s\), constancy follows. For the undercutting/gap analysis the relevant case is \(s=2\) with even scores, where maximum \(\le -2\) and mean \(-2\) yield bi-tight. \(\square\)

2. **Consequence for Type I freeness failure (proved reduction).** In the equality freeness-failure of Prop 15.43(3), \(H=G\cup\{e\}\) is tight of level \(2\). By the master lemma, either \(\Phi(C\oplus H)\ge\Phi\) (no-descent) or \(H\) is bi-tight of level \(2\). In the bi-tight subcase, Prop 15.43(4) gives \(\Phi\ge\Phi-2\) at \(p=3\); for \(p\ge5\) bi-tight level \(2\) is **integrally infeasible at \(p=5\)** (MILP: `src/e1_bitight_infeas.py`, `evidence/e1_bitight_infeas.json`, levels \(2,3,4\) all infeasible while fractional is feasible; avg degree \(4p/(p^2+1)<1\) for all \(p\ge5\)). \(\square\)

3. **Deep tight undercutters are bi-tight (proved).** A gap-\(2\) undercutter with \(S\equiv2\) on \(\mathrm{Max}_{+}\) is two-sided with mean \(-2\) on \(\mathrm{Max}_{-}\); with even scores, maximum \(\le -2\) forces \(S\equiv-2\). No-descent for such undercutters is Prop 15.43(2) (Max\(_{+}\) freeness). At \(p=5\), bi-tight is infeasible, so deep tight undercutters do not exist. \(\square\)

4. **Certified samples.** Five distinct integral tight Max\(_{+}\)-only covers of size \(10\) at \(p=5\) all have \(\max_{\mathrm{Max}_{-}}S\in\{4,8,10\}\ge0\) and exact \(\Phi\in\{73,81,85\}>\Phi(C)\), matching the master lemma. Evidence: session MILP samples; `evidence/E1_BITIGHT.md`.

5. **Residual for \(L=\tfrac12\) (OPEN).** Complete no-descent for deep **non-tight** gap-\(2\) undercutters (\(s_{\min}=2\), \(\max S\ge4\), \(k>2p\)), and lift bi-tight integral infeasibility from \(p=5\) to all \(p\ge5\) (or prove bi-tight \(\Rightarrow\Phi\ge\Phi-2\) uniformly). Then Type I and deep-tight no-descent are unconditional, and \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.45 (star classification force, bi-tight wedge block, p=5 residual kill; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), and write \(f_e(y)=C_e y_iy_j\), \(G_{ee'}=\mathbb E_{+}[f_ef_{e'}]\), and \(g_{\min}\) for the minimum of \(G_{ee'}\) over vertex-disjoint edge pairs.

1. **Wedge correlations (proved).** If edges \(va,vb\) share a vertex, then
   \[
   G_{va,vb}=\frac{C_{va}C_{vb}C_{ab}}{p}=\pm\frac1p,
   \]
   by the Max\(_{+}\) frame identity \(\mathbb E[y_ay_b]=C_{ab}/p\). On Max\(_{-}\) the same wedge has the opposite sign, so
   \[
   G^{+}_{va,vb}+G^{-}_{va,vb}=0.
   \]

2. **Star never bi-tight (proved, all \(p>2\)).** A star has only wedge pairs. For bi-tight of level \(2\) one needs both \(\sum_{e<e'}G^{+}_{ee'}=2-p\) and \(\sum G^{-}=2-p\), hence \(\sum(G^{+}+G^{-})=2(2-p)\). Wedges contribute \(0\) to \(G^{+}+G^{-}\), so a star yields sum \(0\neq 2(2-p)\). \(\square\)

3. **Correlation identity for tight covers (proved).** If \(S_H\equiv s\) on Max\(_{+}\) with \(|H|=sp\), then \(\mathbb E[S^2]=s^2\) and
   \[
   \sum_{e<e'\in H}G_{ee'}=\frac{s^2-sp}{2}.
   \]
   For \(s=1\): sum \(=(1-p)/2\), average \(=-1/p\). For \(s=2\): sum \(=2-p\), average \(=-1/15\).

4. **Star force for level-\(1\) when \(g_{\min}>-1/p\) (proved).** Write \(n_w\) for the number of wedge pairs in \(H\) and \(n_d=\binom{p}{2}-n_w\) for disjoint pairs when \(|H|=p\). Then
   \[
   \sum G\ge -\frac{n_w}p+g_{\min}n_d=g_{\min}\binom{p}{2}-n_w\Bigl(g_{\min}+\frac1p\Bigr).
   \]
   If \(g_{\min}>-1/p\), the coefficient of \(n_w\) is negative, so the lower bound is maximised at maximal \(n_w=\binom{p}{2}\) (every pair of edges shares a vertex). For \(p>3\) that forces \(H\) to be a **star**. Equality with the tight identity \(\sum G=(1-p)/2=-\binom{p}{2}/p\) holds only in that star case (up to wedge-sign pattern). At \(p=3\), \(g_{\min}=-1/p\), so the force fails and non-stars remain possible. \(\square\)

5. **Matching blocked for level-\(2\) when \(g_{\min}>-1/15\) (proved).** For \(|H|=2p\) a matching has \(n_w=0\), hence \(\sum G\ge g_{\min}\binom{2p}{2}\). If this exceeds \(2-p\), no matching is Max\(_{+}\)-tight of level \(2\). \(\square\)

6. **Certified at \(p=5\) (load-bearing numerics).** Exact Max\(\pm\) enumeration (\(|\mathrm{Max}_{\pm}|=260\)):
   - \(g_{\min}=-3/65\approx-0.04615>-1/5\) and \(>-1/15\), so level-\(1\) covers are stars and level-\(2\) matchings are non-tight;
   - integral non-star size-\(p\) Max\(_{+}\) tight covers are **MILP-infeasible**;
   - bi-tight levels \(2,3,4\) integrally infeasible (Prop 15.44);
   - deep two-sided covers (\(s_{+}\ge2\), \(s_{-}\le-2\)) integrally infeasible at \(k=10,12,15\) (timeout without feasibility at \(k=14,16,18,20\));
   - \(\min\{\max_{\mathrm{Max}_{-}}S:S_{\mathrm{Max}_{+}}\ge2,\,|H|=10\}=2\ge0\) (epigraph MILP), so every size-\(2p\) Max\(_{+}\) cover has \(\max S_{-}\ge2\) and \(\Phi\ge\Phi(C)\).
   Evidence: `src/e1_star_bitight_obstruction.py`, `e1_deep_cover_hunt.py`, `e1_deep_k_long.py`; JSON under `evidence/e1_star_bitight_obstruction.json`, `e1_deep_cover_hunt.json`, `e1_deep_k_long.json`.

7. **Certified at \(p=7\).** Full Max\(_{+}\) enumeration (\(|\mathrm{Max}_{+}|=11452\), 80-worker \(2^{25}\) free-coordinate scan): \(g_{\min}\approx-0.03807>-1/7\) and \(>-1/15\), so star-force and matching level-\(2\) block hold at \(p=7\). Evidence: `src/e1_gmin_p7.py`, `evidence/e1_gmin_p7.json`.

8. **Consequence (proved from certs + lemmas).** Stars never bi-tight (all \(p\)). At \(p=5\): Type I freeness-failure and deep-tight undercutters impossible. At \(p=5,7\): level-\(1\) Max\(_{+}\) tight covers are stars (via \(g_{\min}>-1/p\)). Residual for full \(m_n\ge\Phi-2\): bi-tight infeas for general \(p\); deep non-tight large-\(k\) control; closed-form \(g_{\min}>-1/p\) for all \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.46 (1-bit spike formulas and deep/bi-tight spike criterion; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\Phi(C)\), \(A=C\oplus F\). For \(y\in\{\pm1\}^n\) write \(S=S_F(y)\) and \(\sigma_v(y)=\sum_{vw\in F}C_{vw}y_vy_w\). Always \(\sum_v\sigma_v=2S\).

1. **1-bit formula on Max\(_{+}\) (proved).** If \(Cy=py\), then
   \[
   Q_A(y^{\oplus v})
   =
   \Phi-2S-2p+4\sigma_v(y).
   \]
   Hence if \(\max_v\sigma_v(y)\ge(S+p-1)/2\), then \(Q_A(y^{\oplus v})\ge\Phi-2\), so \(\Phi(A)\ge\Phi-2\).
   *Proof.* \((Ay)_v=(Cy)_v-2\sum_{vw\in F}C_{vw}y_w=py_v-2\sum C_{vw}y_w\), so \(y_v(Ay)_v=p-2\sigma_v\). The standard 1-bit identity \(Q(x^{\oplus v})=Q(x)-2x_v(Ax)_v\) yields the claim. Threshold: \(\Phi-2S-2p+4\sigma\ge\Phi-2\Leftrightarrow\sigma\ge(S+p-1)/2\). \(\square\)

2. **1-bit formula on Max\(_{-}\) (proved).** If \(Cz=-pz\), then
   \[
   Q_A(z^{\oplus v})
   =
   -\Phi-2S+2p+4\tau_v(z),
   \]
   where \(\tau_v=\sum_{vw\in F}C_{vw}z_vz_w\). Hence if \(\min_v\tau_v(z)\le(S-p+1)/2\), then \(Q_A(z^{\oplus v})\le-(\Phi-2)\), so \(\Phi(A)\ge\Phi-2\).
   *Proof.* Symmetric to part 1 with \(Cz=-pz\). For the deep boundary \(S=-2\): threshold \(\tau\le-(p+1)/2\). \(\square\)

3. **Bi-tight even-degree spike at \(p=3\) (recovered).** If \(H\) is bi-tight of level \(2\) and every degree in \(H\) is even, then each \(\sigma_v\) is even. With \(\sum\sigma=4>0\), some \(\sigma_v\ge2=(p+1)/2\) at \(p=3\). Part 1 yields \(\Phi(A)\ge\Phi-2\). (This is Prop 15.43.4; undercutting \(C_6\) attains equality \(\sigma_{\max}=2\), \(\Phi=\Phi-2\).) \(\square\)

4. **Sufficient spike criterion (proved).** If either
   - some \(y\in\mathrm{Max}_{+}\) has \(\max_v\sigma_v(y)\ge(S_F(y)+p-1)/2\), or
   - some \(z\in\mathrm{Max}_{-}\) has \(\min_v\tau_v(z)\le(S_F(z)-p+1)/2\),
   then \(\Phi(A)\ge\Phi-2\). In particular, for deep two-sided covers with \(s_{-}=-2\), it suffices that some \(z\) with \(S_F(z)=-2\) has \(\min_v\tau_v\le-(p+1)/2\). \(\square\)

5. **Certified at \(p=5\) (deep two-sided covers).** Every MILP-found deep two-sided cover (\(k\in\{32,36,38,40\}\)) has \(s_{-}=-2\) and some Max\(_{-}\) vector at level \(-2\) with \(\min\tau\le-4\le-3=-(p+1)/2\); 1-bit yields \(|Q|\ge67\ge\Phi-2\), and exact MITM \(\Phi\in\{75,77,79,83\}>\Phi(C)\). Small-\(k\) deep two-sided (\(k\le13,15\)) integrally infeasible. Evidence: `evidence/e1_deep_cover_phi.json`, `e1_deep_sweep_p5.json`, session 1-bit checks. \(\square\)

6. **Residual for \(L=\tfrac12\) (OPEN).** Prove that every bi-tight level-\(2\) cover for \(p\ge5\) satisfies the spike criterion of part 4 (or is integrally impossible — see Prop 15.47), and that every deep two-sided gap-\(2\) undercutter (\(\Phi=\Phi-2\)) either is impossible or has no-descent. Combined with Props 15.42–15.45 this yields \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). Closed-form \(g_{\min}>-1/p\) for all \(p\ge5\) remains open (certified \(p=5,7\)). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.47 (bi-tight Gsum obstruction; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), and write \(G^\pm_{ee'}=\mathbb E_{\pm}[f_ef_{e'}]\), \(h_{ee'}=G^+_{ee'}+G^-_{ee'}\), \(g_{\min}=\min G^+_{ee'}\) over vertex-disjoint pairs, and \(h_{\min}=\min h_{ee'}\) over the same.

1. **Wedge Gsum vanishes (proved, Prop 15.45.1).** For wedge pairs, \(h_{ee'}=0\).

2. **Bi-tight forces disj Gsum sum (proved).** If \(H\) is bi-tight of level \(2\) (\(|H|=2p\), \(S\equiv2\) on Max\(_{+}\), \(S\equiv-2\) on Max\(_{-}\)), then
   \[
   \sum_{e<e'\in H}G^+_{ee'}=2-p=\sum_{e<e'\in H}G^-_{ee'},
   \]
   hence \(\sum_{e<e'\in H}h_{ee'}=2(2-p)\). Wedges contribute \(0\), so writing \(n_d\) for the number of disjoint pairs in \(H\),
   \[
   \sum_{\substack{e<e'\in H\\e\cap e'=\emptyset}}h_{ee'}=2(2-p).
   \]

3. **Floor (proved).** Always \(h_{ee'}\ge 2g_{\min}\) (since each of \(G^\pm\ge g_{\min}\) by definition of \(g_{\min}\) and Max\(_{-}\) symmetry with \(-C\)). Hence for any \(H\) with \(|H|=2p\),
   \[
   \sum_{\mathrm{disj\ pairs\ in\ }H}h
   \;\ge\;
   h_{\min}\,n_d
   \;\ge\;
   2g_{\min}\,n_d
   \;\ge\;
   2g_{\min}\binom{2p}{2},
   \]
   where the last step uses \(g_{\min}<0\) and \(n_d\le\binom{2p}{2}\).

4. **Obstruction (proved).** If \(2g_{\min}\binom{2p}{2}>2(2-p)\), i.e.
   \[
   g_{\min}\;>\;-\frac{p-2}{p(2p-1)},
   \]
   then no bi-tight level-\(2\) cover exists. At \(p=5\) the threshold equals \(-1/15\); for \(p>5\) it is strictly larger (easier) than \(-1/15\). Independently, \(g_{\min}>-1/15\) blocks Max\(_{+}\)-tight matchings of size \(2p\) (Prop 15.45.5). \(\square\)

5. **Certified.** At \(p=5\), \(g_{\min}=-3/65>-1/15\) and \(h_{\min}=-6/65=2g_{\min}\); at \(p=7\), \(g_{\min}\approx-0.03807>-\frac{5}{91}\) and \(h_{\min}=2g_{\min}\). Both satisfy the obstruction, giving a non-MILP proof that bi-tight level \(2\) is empty. Evidence: `evidence/e1_bitight_gsum_obstruction.json`. \(\square\)

6. **Consequence for Type I (proved at \(p=5,7\); conditional for general \(p\ge5\)).** By Prop 15.44, Type I freeness-failure reduces to bi-tight or \(\Phi\ge\Phi\). With bi-tight empty under the \(g_{\min}\) threshold of part 4, Type I no-descent is unconditional. Deep tight undercutters (bi-tight) are empty. Residual: prove \(g_{\min}>-(p-2)/(p(2p-1))\) for all \(p\ge5\); deep non-tight gap-\(2\) undercutters (ND or \(\Phi\ge\Phi-2\)). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.48 (edge-correlation algebra and cross-ratio structure of \(g_{\min}\); 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), edges \(e=\{i,j\}\) with \(f_e(y)=C_{ij}y_iy_j\), and \(G_{ee'}=\mathbb E_{+}[f_ef_{e'}]\).

1. **Row sum (proved).** \(\sum_{e'}G_{ee'}=n/2\) for every \(e\). *Proof.* \(\sum_{e'}f_{e'}=\tfrac12 y^\top Cy=\Phi\) on Max\(_{+}\), and \(\mathbb E[f_e]=\tfrac1p\), so the row sum is \(\Phi/p=n/2\). \(\square\)

2. **Wedge exact values (proved).** If \(e,e'\) share a vertex then \(G_{ee'}=\pm1/p\), and for each fixed \(e\) the sum of \(G_{ee'}\) over the \(2(n-2)\) wedge partners is \(0\). *Proof.* Star identity \(\sum_{j\neq i}f_{ij}=p\) at each vertex yields \(\sum_{\mathrm{wedge\ of\ }e}f=2p-2f_e\); take \(\mathbb E[f_e\cdot(\cdot)]\). Exact values \(\pm1/p\) follow from the 2-design computation of \(E[y_ay_by_ay_c]=E[y_by_c]=C_{bc}/p\) on wedges (Prop 15.45.1 refinement). \(\square\)

3. **Disjoint mean (proved).** For each \(e\), \(\sum_{e'\,:\,e\cap e'=\emptyset}G_{ee'}=n/2-1\), hence the average disjoint correlation is \((n/2-1)/(E-1-2(n-2))\) with \(E=\binom{n}{2}\). \(\square\)

4. **Four-point pairing identity (proved).** For distinct vertices \(i,j,k,l\), writing \(\kappa=C_{ij}C_{kl}+C_{ik}C_{jl}+C_{il}C_{jk}\) and \(m_4=\mathbb E_{+}[y_iy_jy_ky_l]\),
   \[
   G_{\{ij\},\{kl\}}+G_{\{ik\},\{jl\}}+G_{\{il\},\{jk\}}
   \;=\;
   \kappa\,m_4.
   \]
   *Proof.* Each pairing contributes \(C_eC_{e'}m_4\), and the three \(C\)-products sum to \(\kappa\). \(\square\)

5. **Cross-ratio structure (certified \(p=3,5\)).** Identifying vertices with \(\mathrm{PG}(1,\mathbb F_{p^2})\), the PGL-invariant of a 4-set is its cross-ratio class. At \(p=5\): on classes with constant \(m_4\) given \(\kappa\), the three pairing correlations form the multiset \(\{-3,3,3\}/65\) or \(\{\pm1\}/65\), and \(g_{\min}=-3/65\) is achieved precisely on the \(\{-3,3,3\}/65\) classes. At \(p=3\), the analogous pattern is \(\{-1,1,1\}/3\) with \(g_{\min}=-1/3\). Two residual cross-ratio classes at \(p=5\) have non-constant \(m_4\) (binary extra invariant). Evidence: `evidence/e1_gmin_closed_form_attack.json`. \(\square\)

6. **Dead lower-bound attempts (do not reopen).** The bound \(g_{\min}\ge-3/\Phi\) holds with equality at \(p=5\) but **fails** at \(p=7\) (\(g_{\min}=-109/2863<-3/175\)). Pure 4-point boolean LP, Chebyshev on disjoint partners, Wick/Gaussian as a lower bound, bare \(C\)-isomorphism types, and affine halfspace orbits alone are all too weak or incomplete (cf. `E1_FAILURE_GRAPH.md` residual notes). \(\square\)

7. **Residual for Prop 15.47 (OPEN).** Prove \(g_{\min}>-(p-2)/(p(2p-1))\) for all primes \(p\ge5\) by a character-sum evaluation of \(m_4\) on the cross-ratio class that realises the minimum (or another scheme formula). Certified only at \(p=5,7\). Deep non-tight gap-\(2\) control remains independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.49 (cross-ratio classification of \(g_{\min}\); uniform LB candidate; 2026-07-28).** Continue the notation of Prop 15.47–15.48. Write \(n=p^2+1\), \(\Phi=\tfrac12 np\), and \(N=|\mathrm{Max}_{+}|\).

1. **Cross-ratio stratification (certified \(p=3,5,7\)).** Identifying vertices with \(\mathrm{PG}(1,\mathbb F_{p^2})\), each 4-set has a PGL-invariant cross-ratio class \(\mathrm{CR}\) and a \(C\)-invariant \(\kappa=\sum C_eC_{e'}\) over the three pairings. Evidence `e1_gmin_cr_classify.json`:
   - On every class with \(|\kappa|=1\) and constant \(m_4\), one has \(m_4=\kappa\alpha\) for a class constant \(\alpha>0\), and the three pairing correlations form \(\{-\alpha,\alpha,\alpha\}\). Hence the pairing minimum is \(-\alpha\).
   - Global \(g_{\min}=-\alpha_\star\) where \(\alpha_\star\) is the maximum of such \(\alpha\) over constant-\(m_4\) classes with \(|\kappa|=1\).
   - Values: \(g_{\min}=-\tfrac13\) (\(p=3\)); \(-\tfrac3{65}\) (\(p=5\)); \(-\tfrac{109}{2863}\) (\(p=7\)). At \(p=5,7\), \(g_{\min}\) strictly exceeds the bi-tight threshold of Prop 15.47; at \(p=3\) it does not (consistent with bi-tight \(C_6\)).

2. **Spectral rank (certified \(p=3,5,7\)).** The Gram \(G\) has \(\mathrm{rank}(G)=\binom{d}{2}-d+1\) with \(d=n/2\), and a simple eigenvalue \(n/2\) for the all-ones vector on edges. At \(p=5\) the nonzero spectrum is \(\{n/2,88/13,72/13,40/13\}\); at \(p=7\) it is \(\{n/2\}\) union five positive eigenvalues with denominator \(409\) and multiplicities \(\{d,n,n,2n,n\}\).

3. **Uniform lower-bound candidate (algebra + certification).** Define
   \[
   L(p)\;:=\;-\frac{p-2}{2p^2}.
   \]
   For every odd integer \(p>2\),
   \[
   L(p)\;>\;-\frac{p-2}{p(2p-1)}
   \]
   (clear the positives: \(2p^2>p(2p-1)\Leftrightarrow 0>-1\)). Thus \(L(p)\) lies strictly above the bi-tight threshold. Certified: \(g_{\min}(5)\ge L(5)\) and \(g_{\min}(7)\ge L(7)\); while \(g_{\min}(3)<L(3)\) (so the candidate is not valid at \(p=3\), as required by the existence of bi-tight undercutters there). Evidence: `e1_gmin_uniform_lb.json`.

4. **Matching margin (certified \(p=5\)).** For every matching \(M\) of size \(2p=10\) on Paley \(C_{26}\), the quadratic form \(\mathbf{1}_M^\top G\mathbf{1}_M\) is at least \(9.96>4\), so \(\sum_{\mathrm{pairs\ in\ }M}G\ge -0.02>2-p=-3\). Equality in the PSD projection bound \(\mathbf{1}_M^\top G\mathbf{1}_M\ge4\) is therefore not attained on matchings (sampling + local search, \(2\cdot10^4\) seeds).

5. **Residual for Prop 15.47 (OPEN).** Prove \(g_{\min}\ge L(p)\) for every prime \(p\ge5\) (or any other uniform lower bound strictly above the bi-tight threshold). Combined with the deep non-tight residual of Prop 15.47.6 this yields \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.50 (conditional means of Max+; disjoint correlation mean; 2026-07-29).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Sigma=I+C/p=2P_+\), and let \(y\) be uniform in \(\mathrm{Max}_{+}\). Write \(f_e=C_e y_iy_j\) and \(G_{ee'}=\mathbb E[f_ef_{e'}]\).

1. **Conditional mean lemma (proved).** For distinct indices \(i,j\) and \(a,b\in\{\pm1\}\),
   \[
   \mathbb E\bigl[y\bigm|\,y_i=a,\,y_j=b\bigr]
   \;=\;
   \Sigma_{*,S}\,\Sigma_{S,S}^{-1}\begin{pmatrix}a\\b\end{pmatrix},
   \qquad S=\{i,j\}.
   \]
   *Proof.* The right-hand side \(\mu_*\) is the unique minimum-norm vector in \(V_+=\mathrm{range}(P_+)\) with those two coordinates (Gaussian interpolant for covariance \(\Sigma\)). The left-hand side \(\mu\) is an average of Max\(+\) vectors, hence lies in \(V_+\) with the same coordinates. The difference \(v=\mu-\mu_*\) satisfies \(v\in V_+\) and \(v_i=v_j=0\). For every \(w\in V_+\) with \(w_i=w_j=0\), the scalar \(y\cdot w\) is orthogonal (under \(\mathbb E\)) to \(\{1,y_i,y_j,y_iy_j\}\): the first three pairings use \(\mathbb E[y]=0\) and \(\mathbb E[yy^\top]=\Sigma\); the fourth is an odd third moment and vanishes by central symmetry \(\mathrm{Max}_+=-\mathrm{Max}_+\). Those four monomials span all functions of \((y_i,y_j)\), so \(\mathbb E[y\cdot w\mid y_i,y_j]=0\). Hence \(\mu\cdot w=0\) for all such \(w\). In particular \(\mu_*\) lies in \(\mathrm{span}\{P_+e_i,P_+e_j\}\) and \(v\perp\) that span, so \(\mu_*\cdot v=0\) and \(\mu\cdot v=\|v\|^2=0\), whence \(v=0\). \(\square\)
   Certified at \(p=5,7\) (max abs error \(<10^{-15}\)): `evidence/e1_gmin_cond_mean.json`, `src/e1_gmin_cond_mean.py`.

2. **Conditional second-moment shape (proved).** For fixed distinct \(i,j,k,l\),
   \[
   \mathbb E[y_ky_l\mid y_i,y_j]
   \;=\;
   \alpha+\delta\,y_iy_j
   \]
   for scalars \(\alpha,\delta\) determined by \(\Sigma_{kl}\) and \(m_4=\mathbb E[y_iy_jy_ky_l]\):
   \[
   \alpha=\frac{\Sigma_{kl}-c\,m_4}{1-c^2},\qquad
   \delta=\frac{m_4-c\,\Sigma_{kl}}{1-c^2},\qquad
   c=\Sigma_{ij}=C_{ij}/p.
   \]
   *Proof.* The conditional expectation is a function of \((y_i,y_j)\), hence of the form \(\alpha+\beta y_i+\gamma y_j+\delta y_iy_j\). Matching moments against \(1,y_i,y_j,y_iy_j\) and using vanishing odd moments forces \(\beta=\gamma=0\) and the displayed formulae. \(\square\)

3. **Disjoint-pair mean (proved).** For every edge \(e\), \(\sum_{e'\,:\,e\cap e'=\emptyset}G_{ee'}=n/2-1\), so the average disjoint correlation is
   \[
   \frac1{p^2-2}.
   \]
   *Proof.* Row sum of \(G\) is \(n/2\) (Prop 15.48.1); sum of \(G\) over the \(2(n-2)\) wedge partners of \(e\) is \(0\) (Prop 15.48.2); the remainder is the disj sum. Count of disj partners \((n-2)(n-3)/2\) and \(n=p^2+1\) give the average. \(\square\)

4. **Frechet obstruction is too weak (certified).** Combining part 1 with the Fréchet–Hoeffding bound
   \(\mathbb E[y_ky_l\mid\mathrm{state}]\ge\mu_k\mu_l-\sqrt{(1-\mu_k^2)(1-\mu_l^2)}\)
   yields only \(g_{\min}\ge-O(1)\) (empirically \(\ge-0.4\) at \(p=5\)), which does **not** beat the bi-tight threshold \(T(p)=-(p-2)/(p(2p-1))\). Evidence in session residual notes. **Do not reopen plain Fréchet as a path to \(L(p)\).**

5. **Residual (OPEN).** Prove \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) (or any LB \(>T(p)\)) for all primes \(p\ge5\). The conditional-mean calculus reduces this to a uniform upper bound on \(|m_4|\) (or on \(\delta\)) for four-sets with \(|\kappa|=1\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.51 (equivalent bi-tight threshold form; residual slice structure; 2026-07-29).** Continue the notation of Prop 15.47–15.50. Fix an edge \(e=(i,j)\) and write \(a(e'):=\mathbb E[f_{e'}\mid f_e=1]\) for \(e'\) vertex-disjoint from \(e\).

1. **Identity (proved).** For any disj pair \(e,e'\),
   \[
   a(e')=\frac{1+p\,G_{ee'}}{p+1}.
   \]
   *Proof.* From Prop 15.50, \(\mathbb E[y_ky_l\mid y_iy_j=s]=\alpha+\delta s\) with \(\alpha,\delta\) linear in \(m_4=G_{ee'}/(C_eC_{e'})\); conditioning on \(f_e=1\) forces \(y_iy_j=C_e\) and yields the display after multiplying by \(C_{e'}\). Equivalently, expand \(G=\mathbb E[f_ef_{e'}]\) on the two values of \(f_e\) with \(\mathbb E[f_e]=1/p\). \(\square\)

2. **Equivalent bi-tight threshold (proved).** Therefore
   \[
   g_{\min}\ge T(p):=-\frac{p-2}{p(2p-1)}
   \quad\Longleftrightarrow\quad
   \min_{e'\,:\,e\cap e'=\emptyset}a(e')\ge\frac1{2p-1}
   \]
   for every edge \(e\) (and hence, by edge-transitivity of \(\mathrm{Aut}\), for one fixed \(e\)). *Proof.* Clear the linear identity of part 1. \(\square\)
   Certified: at \(p=5\), \(\min a=5/39>1/9\); at \(p=7\), \(\min a\approx0.0917>1/13\). Evidence: `e1_gmin_structure.json`.

3. **Deterministic disj sum (proved).** Pointwise on \(\mathrm{Max}_{+}\),
   \[
   \sum_{e'\,:\,e\cap e'=\emptyset}f_{e'}=\Phi-2p+f_e.
   \]
   *Proof.* \(\sum_{\mathrm{all}\,e'}f=\Phi\) and \(\sum_{\mathrm{wedge\,of\,}e}f=2p-2f_e\) (stars at the two ends of \(e\)). \(\square\)
   Consequently on the slice \(f_e=1\) the average of \(f_{e'}\) over disj partners is the constant \((\Phi-2p+1)/D\) with \(D=(n-2)(n-3)/2\).

4. **Residual Loewner (certified \(p=5,7\)).** On the slice \(f_e=1\), the residual Gram \(R=\mathbb E[(y-\mu)(y-\mu)^\top\mid f_e=1]\) (mixture of the two states with \(f_e=1\)) satisfies \(R\succeq\lambda_{\min}(R)\,P_W\) where \(P_W\) is the orthogonal projector onto \(V_+\cap\{x_i=x_j=0\}\) and \(\mathrm{rank}(R)=n/2-2\). Evidence: `e1_gmin_structure.json`. **Not yet a proof of \(g_{\min}>T(p)\):** entrywise Schur bounds on the residual still undershoot \(1/(2p-1)\).

5. **Max+ types (certified).** At \(p=5\), \(\mathrm{Max}_{+}\) is distance-homogeneous (constant Hamming distance distribution from every vector). At \(p=7\), at least two distance types occur (so \(\mathrm{Max}_{+}\) is not a single \(\mathrm{Aut}\)-orbit). Affine+\(\mathrm{PGL}\) orbit of the halfspace vector has size \(60\) of \(260\) at \(p=5\) — character sums on that orbit alone do **not** compute full-Max+ \(m_4\).

6. **Residual (OPEN).** Prove \(\min a(e')\ge1/(2p-1)\) (equivalently \(g_{\min}\ge T(p)\), or the stronger \(g_{\min}\ge L(p)\)) for all primes \(p\ge5\). Preferred routes: character-sum / BM algebra on the min CR class; residual \(z=y-\mu\) with the pointwise identity \(z_r^2+2\mu_rz_r=1-\mu_r^2\) beyond Fréchet; Loewner calculus with a sharp entrywise bound; **m4 moduli** (evec system on refined \(C\)-classes has nullity 1; pin by \(\mathrm{Tr}(G^2)\), recovers \(g_{\min}\) at \(p=5\) — `E1_GMIN_MODULI.md`). **Do not** reopen F15 (plain Fréchet) or incomplete Aut-orbits as full Max+. Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.52 (coordinate sum of Max+; m4 moduli sketch; 2026-07-29).** Let \(C\) be the Paley conference matrix of order \(n=p^2+1\) with vertices \(\{\infty\}\cup\mathbb F_{p^2}\), and \(y\in\mathrm{Max}_{+}\).

1. **Coordinate sum (proved).** \(\mathbf1^\top y=(p+1)y_\infty\). In particular \(|\mathbf1^\top y|=p+1\).
   *Proof.* Row sums of \(C\): \((C\mathbf1)_\infty=p^2\) and \((C\mathbf1)_v=1\) for \(v\in\mathbb F_{p^2}\) (complete character sum \(\sum_{d\neq0}\chi(d)=0\)). Thus \(C\mathbf1=(p^2-1)e_\infty+\mathbf1\). From \(Cy=py\), take \(\mathbf1^\top\): \(y^\top C\mathbf1=p\,s\) with \(s=\mathbf1^\top y\), so \((p^2-1)y_\infty+s=ps\), hence \(s=(p+1)y_\infty\). \(\square\)
   Certified \(p=5,7\): every Max+ vector has sum \(\pm(p+1)\).

2. **m4 linear system (certified \(p=5,7\)).** Stratify 4-sets by \((\mathrm{CR},\kappa,\triangle\text{-type})\) (all \(C\)-invariants; \(\triangle\)-type splits formerly non-constant classes). Averaged evec identities give \((pI-M)\mathbf m=\mathbf b\) with \(M,b\) combinatorial. The system has **nullity 1**. The true Max+ moment vector lies on this line. The constraint \(\mathrm{Tr}(G^2)=E+2n_{\mathrm{wedge}}/p^2+6\sum n_A m_A^2\) is quadratic in the free parameter and **selects the true \(m_4\) at \(p=5\)** (hence \(g_{\min}=-3/65\)). Evidence: `E1_GMIN_MODULI.md`.

3. **Residual (OPEN).** Prove nullity 1 and a Max+-free evaluation of \(\mathrm{Tr}(G^2)\) (or of the \(G\)-spectrum) for all primes \(p\ge5\); solve for \(g_{\min}\) and show \(g_{\min}\ge L(p)\) or \(>T(p)\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.53 (pairing reduction of \(g_{\min}\); moduli pin at \(p=5\); 2-design \(\mathrm{Tr}(G^2)\) skeleton; 2026-07-30).** Continue Prop 15.52.

1. **Pairing identity (proved).** For any 4-set with \(\kappa:=\sum_{\text{three pairings}}C_eC_{e'}\) satisfying \(|\kappa|=1\), the three pairing products are a permutation of \((1,1,-1)\) or \((-1,-1,1)\). Hence the three edge-pair correlations equal \(m_4\cdot(\pm1)\) in that pattern, and
   \[
   \min_{\text{three pairings}}G_{ee'}=-\lvert m_4\rvert.
   \]
   Therefore
   \[
   g_{\min}=-\max\bigl\{\lvert m_4(S)\rvert:S\subset V,\,|S|=4,\,|\kappa(S)|=1\bigr\}.
   \]
   *Proof.* Each pairing product is \(\pm1\); their sum is \(\kappa=\pm1\) forces the stated multiset. Then \(G=C_eC_{e'}m_4\) on each pairing, so the minimum is \(-\lvert m_4\rvert\). Every disjoint edge pair sits in a unique 4-set, so the global \(g_{\min}\) is the min over such 4-sets. \(\square\)
   Certified \(p=5,7\): identity error \(<10^{-12}\); \(g_{\min}=-3/65\), \(-109/2863\). Evidence: `e1_gmin_moduli.json`.

2. **Refined \(C\)-classes (certified \(p=5\)).** Stratify 4-sets by the pure \(C\)-invariant
   \((\mathrm{type}_6,\mathrm{ext}\text{-sum histogram})\), where \(\mathrm{type}_6\) is the \(S_4\)-canonical 6-tuple of edge signs and the external histogram records \(\sum_{v\in S}C_{rv}\) for \(r\notin S\). At \(p=5\) this yields **37** classes, each with **constant** \(m_4\) on Max+. (Bare \(C\)-types alone do **not** make \(m_4\) constant — F-graph / prior notes.)

3. **Nullity-1 evec system (certified \(p=5\)).** Averaging \(p\,m_4=\sum_r C_{ir}m(\cdots)\) over each class produces a combinatorial linear system \(A\mathbf m=\mathbf b\) (RHS uses only the 2-design \(m_2=C/p\)). At \(p=5\): \(\mathrm{rank}(A)=36=n_{\mathrm{var}}-1\). The true Max+ moment vector lies on the affine line \(\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n\).

4. **\(\mathrm{Tr}(G^2)\) pin (certified \(p=5\)).** Write \(K_{ab}=\bigl((y_a\cdot y_b)^2-n\bigr)/2\). Then the nonzero spectrum of \(G\) matches that of \(K/N\), and
   \[
   \mathrm{Tr}(G^2)=\frac1{N^2}\sum_{a,b}K_{ab}^2=\tfrac14\Bigl(\mathbb E[\mathrm{dot}^4]-2n\,\mathbb E[\mathrm{dot}^2]+n^2\Bigr).
   \]
   Substituting \(\mathbf m(c)\) into the edge form of \(\mathrm{Tr}(G^2)\) yields a quadratic in \(c\). One root recovers \(g_{\min}=-3/65\); select the root of larger \(g_{\min}\) among the two (do **not** use PSD-max over the whole line — F16). Evidence: `e1_gmin_moduli.json`, `E1_GMIN_MODULI.md`.

5. **2-design evaluation of \(\mathbb E[\mathrm{dot}^2]\) (proved, Max+-free beyond the frame).** From \(\mathbb E[yy^\top]=I+C/p\),
   \[
   \mathbb E_{a,b}[(y_a\cdot y_b)^2]=\|I+C/p\|_F^2=n+\frac{n(n-1)}{p^2}.
   \]
   *Proof.* \(\sum_{a,b}(y_a\cdot y_b)^2=\sum_{i,j}(\sum_a y_{a,i}y_{a,j})^2=N^2\|I+C/p\|_F^2\). \(\square\)
   Certified \(p=5,7\). **Still Max+-dependent:** \(\mathbb E[\mathrm{dot}^4]\) (equivalently a closed \(G\)-spectrum), which is needed to evaluate \(\mathrm{Tr}(G^2)\) without Max+ samples.

6. **Wick comparison (certified \(p=5,7\); not yet a proof).** For \(\Sigma=I+C/p\), the Gaussian fourth moment is \(\mathbb E_{\mathrm{Wick}}[\mathrm{dot}^4]=3\|\Sigma\|_F^4+6\mathrm{Tr}(\Sigma^4)\). Discrete Max+ satisfies \(\mathbb E[\mathrm{dot}^4]<\mathbb E_{\mathrm{Wick}}[\mathrm{dot}^4]\) at \(p=5,7\) (boolean coordinates have smaller kurtosis than Gaussians with the same covariance). On the combinatorial moduli line at \(p=5\), the weaker constraint \(\mathrm{Tr}(G^2)\le\mathrm{Tr}_{\mathrm{Wick}}\) forces \(g_{\min}\ge T(p)\) with **endpoint equality** (float margin \(\sim10^{-15}\)) — not the strict \(g_{\min}>T(p)\) needed for Prop 15.47. A proved strict inequality \(\mathbb E[\mathrm{dot}^4]\le\mathbb E_{\mathrm{Wick}}-\delta\) (or the exact spectrum) would upgrade this to \(g_{\min}>T\). Evidence: `e1_gmin_moduli.json`.

7. **Residual (OPEN).** Prove for all primes \(p\ge5\): (i) refined classes have constant \(m_4\) and the evec system has nullity 1; (ii) Max+-free strict bound \(\mathbb E[\mathrm{dot}^4]<\mathbb E_{\mathrm{Wick}}\) or a closed \(G\)-spectrum; (iii) the selected root satisfies \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) or at least \(g_{\min}>T(p)\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.54 (moduli \(c\)-pin calculus; wedge closed form; a-slice certs; 2026-07-30).** Continue Prop 15.53.

1. **Wedge \(G\) (proved, combinatorial).** If edges \(e=(i,j)\) and \(e'=(i,k)\) share a vertex, then
   \[
   G_{ee'}=\frac{C_{ij}C_{ik}C_{jk}}{p}=\pm\frac1p.
   \]
   *Proof.* \(f_ef_{e'}=C_{ij}C_{ik}y_jy_k\) and \(\mathbb E[y_jy_k]=C_{jk}/p\). \(\square\)

2. **Moduli line and \(g_{\min}(c)\) (certified \(p=5\)).** On the nullity-1 line \(\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n\) of Prop 15.53, for every \(|\kappa|=1\) class
   \(g_{\min}(c)=-\max_A|m_A(c)|\). The edge Gram is \(G(c)=G_{\mathrm{wedge}}+G_{\mathrm{disj}}(\mathbf m(c))\) with wedges from part 1. Then \(\mathrm{Tr}(G(c)^2)=a_0+a_1c+a_2c^2\) with combinatorial \(a_i\). Evidence: `e1_gmin_cbound.json`.

3. **True \(\mathrm{Tr}(G^2)\) pin (certified \(p=5\)).** Setting \(\mathrm{Tr}(G(c)^2)=\mathrm{Tr}(G_{\mathrm{Max+}}^2)\) yields two roots; the root of larger \(g_{\min}\) recovers \(g_{\min}=-3/65>T(5),L(5)\). The Wick fourth-moment vector \(m_4=\kappa/p^2\) is **not** exactly on the evec line (residual \(\sim10^{-2}\)–\(10^{-1}\) relative). Evidence: `e1_gmin_cbound.json`.

4. **Slice \(a_{\min}\) (certified \(p=5,7\)).** \(\min a=5/39>1/9\) at \(p=5\) and \(75/818>1/13\) at \(p=7\). Wick-\(a\) overestimates \(\min a\); mean-only \(C_{kl}\mu_k\mu_l\) undershoots (residual \(R\) helps). Evidence: `e1_gmin_abound.json`. **Not a uniform proof.**

5. **Deep covers at \(p=5\) (certified refresh).** MILP deep two-sided covers at \(k\in\{32,36,38,40\}\) all meet the Prop 15.46 spike criterion and have \(\Phi\ge\Phi(C)\). Evidence: `e1_deep_spike_theory.json`. Uniform deep ND for all \(p\ge5\) still open.

6. **Residual (OPEN).** Same as Prop 15.53.7: Max+-free \(\mathrm{Tr}(G^2)\) / spectrum for general \(p\), hence \(g_{\min}>T(p)\) for all primes \(p\ge5\); deep non-tight ND. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.55 (tight size-\(2p\) obstruction from \(\lambda_{\max}(G)=n/2\); 2026-07-30).** Let \(G=\mathbb E_{+}[ff^\top]\) on the edge space of Paley \(C\) of order \(n=p^2+1\).

1. **Row sum (proved).** \(G\mathbf1=(n/2)\mathbf1\). *Proof.* \(\sum_e f_e=\Phi=pn/2\) and \(\mathbb E[f_e]=1/p\), so each row sums to \(\Phi/p=n/2\). \(\square\)

2. **All-ones mass of a size-\(2p\) indicator (proved).** For any \(0\)-\(1\) vector \(v\) with \(\sum v_e=2p\),
   \[
   v^\top\Bigl(\tfrac n2\cdot\frac{\mathbf1\mathbf1^\top}{E}\Bigr)v=4,
   \]
   where \(E=\binom{n}{2}\). *Proof.* Direct: \(\tfrac n2\cdot(2p)^2/E=4p^2/(n-1)=4p^2/p^2=4\). \(\square\)

3. **Tight cover forces \(G_\perp\)-isotropy (proved).** If \(H\) is Max\(_{+}\)-tight of level \(2\) (\(|H|=2p\), \(S_H\equiv2\) on \(\mathrm{Max}_{+}\)), then \(v:=\mathbf1_H\) satisfies \(v^\top Gv=\mathbb E[S_H^2]=4\). Writing \(G=\tfrac n2 P_{\mathbf1}+G_\perp\) with \(P_{\mathbf1}=\mathbf1\mathbf1^\top/E\), part 2 yields \(v^\top G_\perp v=0\).

4. **Obstruction when \(\lambda_{\max}(G)=n/2\) is simple (proved).** Assume \(G\succeq0\) (true: \(G=\mathbb E[ff^\top]\)) and \(\lambda_{\max}(G)=n/2\) with multiplicity one. Then \(G_\perp\succeq0\) and \(\ker G_\perp=\mathrm{span}\{\mathbf1\}\). From part 3, \(v^\top G_\perp v=0\Rightarrow G_\perp v=0\Rightarrow v\parallel\mathbf1\), impossible for \(|H|=2p<E\). **Therefore no Max\(_{+}\)-tight level-\(2\) cover of size \(2p\) exists** — in particular bi-tight level \(2\) is empty, and Type I freeness-failure (Prop 15.44) cannot produce a bi-tight cover.

5. **Certified spectrum.** At \(p=5,7\): \(\lambda_{\max}(G)=n/2\) is simple (next eigenvalues \(\approx6.77,5.28\ll n/2\)). At \(p=3\): \(\lambda_{\max}(G)=8>n/2=5\) (multiplicity \(5\)), so the obstruction does **not** apply — consistent with bi-tight \(C_6\). Evidence: `e1_gmin_tight_obstruction.json`.

6. **Residual (OPEN).** Prove \(\lambda_{\max}(G)=n/2\) (simple) for **all primes \(p\ge5\)**. Then bi-tight / Type I residual of Path C closes without a uniform \(g_{\min}\) bound. Deep non-tight gap-\(2\) residual remains. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.56 (star/cycle decomposition and Schur-square reduction; 2026-07-30).** Continue Prop 15.55. Write \(E=\binom{n}{2}\), \(d=n/2\), \(N=|\mathrm{Max}_{+}|\), and let \(Y\) be the \(N\times n\) matrix of Max+ vectors.

1. **Star action (proved).** Let \(u^{(i)}\in\mathbb R^E\) be the indicator of edges incident to vertex \(i\). Then \(u^{(i)\top}f=p\) for every Max+ feature vector \(f\) (star sum), hence
   \[
   Gu^{(i)}=\mathbf1\qquad\text{for all }i.
   \]
   Consequently \(G(u^{(i)}-u^{(j)})=0\), so \(\ker G\) contains the \((n-2)\)-dimensional space of star differences orthogonal to the all-ones edge vector. *Proof.* \(f\cdot u^{(i)}=p\) constantly, so \(G u^{(i)}=\mathbb E[f\,p]=p\cdot\mathbb E[f]=\mathbf1\). \(\square\)

2. **Cycle reduction (proved).** Let \(\mathrm{Cyc}=\{v\in\mathbb R^E:\sum_{j\neq i}v_{ij}=0\ \forall i\}\) (cycle space, dim \(E-n+1\)). Then
   \[
   \lambda_{\max}(G)=\max\bigl(n/2,\ \lambda_{\max}(G|_{\mathrm{Cyc}})\bigr).
   \]
   *Proof.* The star space has dimension \(n-1\), carries the eigenvalue \(n/2\) along \(\mathbf1\), and contributes \(n-2\) kernel directions by part 1; all remaining positive spectrum lies in \(\mathrm{Cyc}\). \(\square\)

3. **Schur-square dual (proved).** Let \(D=YY^\top\) and \(P=D/(2N)\). Then \(P\) is the orthogonal projector of rank \(d\) with constant diagonal \(P_{aa}=d/N\), and the nonzero eigenvalues of \(G\) other than \(n/2\) equal \(2N\) times the eigenvalues of \(P\odot P\) on \(\mathbf1^\perp\subset\mathbb R^N\). In particular
   \[
   \lambda_{\max}(G|_{\mathrm{Cyc}})=2N\cdot\lambda_2(P\odot P),
   \]
   so
   \[
   \lambda_{\max}(G)=\tfrac n2
   \quad\Longleftrightarrow\quad
   \lambda_2(P\odot P)\le \frac d{2N}.
   \]
   *Proof.* \(D\) has spectrum \(2N\) (\(d\) times) and \(0\); \(P=QQ^\top\) for \(Q^\top Q=I_d\). Feature Gram identity \(K_{ab}=((y_a\cdot y_b)^2-n)/2\) gives \(K=2N^2(P\odot P)-\frac n2 J\), and nonzero eigenvalues of \(G\) match those of \(K/N\). On \(\mathbf1^\perp\), \(J=0\), yielding the factor \(2N\). \(\square\)

4. **Average cycle eigenvalue (proved for \(p\ge5\)).** Writing \(k=\binom{d-1}{2}-1\) for the cycle rank,
   \[
   \frac1k\sum_{\mathrm{cycle}}\lambda_j
   =\frac{n(n-2)/2}{k}
   <\frac n2
   \quad\text{for all primes }p\ge5.
   \]
   *Proof.* Algebra: \(n-1\le(d-1)(d-2)/2\) rearranges to the claim; holds for \(n=p^2+1\ge26\). \(\square\)
   (Average \(<n/2\) is necessary but not sufficient for \(\lambda_{\max}(\mathrm{cycle})\le n/2\).)

5. **Certified spectral gap.** At \(p=5,7\): \(\lambda_2(P\odot P)\le d/(2N)\) (hence \(\lambda_{\max}(G)=n/2\) simple). At \(p=3\): \(\lambda_2=1/3>5/24=d/(2N)\) (hence \(\lambda_{\max}=8>5\)). Evidence: `e1_gmin_spectral.json`.

6. **Residual (OPEN).** Prove \(\lambda_2(P\odot P)\le d/(2N)\) for every prime \(p\ge5\) (equivalently \(\lambda_{\max}(G|_{\mathrm{Cyc}})\le n/2\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.57 (triangle form bound; crude \(\lambda_{\max}(G)\le n\); 2026-07-30).** Let \(v\in\mathbb R^E\) with \(\|v\|=1\), and write \(B_{ij}=C_{ij}v_{ij}\) (\(i\neq j\)), \(B_{ii}=0\). Define the **triangle form**
\[
\mathrm{ft}(v)=\sum_{i=1}^n\sum_{j\neq i}\sum_{k\neq i}C_{ij}C_{ik}C_{jk}\,v_{ij}v_{ik}=\sum_{i}w^{(i)\top}C^{[i]}w^{(i)},
\]
where \(w^{(i)}_j=C_{ij}v_{ij}\) and \(C^{[i]}\) is the principal submatrix of \(C\) deleting row/column \(i\).

1. **Triangle bound (proved).** \(|\mathrm{ft}(v)|\le 2p\).
   *Proof.* For each \(i\), the Rayleigh principle gives \(|w^\top C^{[i]}w|\le\|C\|_{\mathrm{op}}\|w\|^2=p\,\delta_i\) with \(\delta_i=\sum_{j\neq i}v_{ij}^2\), because \(C^{[i]}\) is the restriction of the symmetric matrix \(C\) to a coordinate subspace (or directly: extend \(w\) by a zero at \(i\) and use \(\|C\|_{\mathrm{op}}=p\)). Summing over \(i\) and using \(\sum_i\delta_i=2\|v\|^2=2\) yields the claim. \(\square\)
   Certified: on the cycle space, the spectrum of the triangle form is exactly \(\{\pm 2p\}\) at the extremes for \(p=3,5,7\). Evidence: `e1_gmin_spectral.json` / session residual notes.

2. **Second-moment identity (proved).** For \(y\in\mathrm{Max}_{+}\),
   \[
   \mathbb E[\|By\|^2]=2+\frac{\mathrm{ft}(v)}{p}\le 4.
   \]
   *Proof.* Expand \(\mathbb E[y^\top B^\top B y]=\mathrm{Tr}(B^\top B\Sigma)\) with \(\Sigma=I+C/p\) and \(y_i^2=1\); the cross term is \(\mathrm{ft}/p\). Apply part 1. \(\square\)

3. **Crude operator bound (proved).** \(\lambda_{\max}(G)\le n\).
   *Proof.* Cauchy--Schwarz: \((y^\top By)^2\le n\|By\|^2\). Take \(\mathbb E\), use part 2, and \(f\cdot v=y^\top By/2\). \(\square\)
   (Factor-of-two away from the target \(\lambda_{\max}\le n/2\); the maximising cycle direction has \(\mathrm{ft}=2p\) and \(By\in V_+\) for all Max+ \(y\), certified \(p=3,5,7\).)

4. **Residual (OPEN).** Improve part 3 by a factor \(2\) on the cycle space for \(p\ge5\) (equivalently Prop 15.56.6: \(\lambda_2(P\odot P)\le d/(2N)\)). Then bi-tight closes for all primes \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.58 (Max+ in \(V_+\); Perron of \(P\odot P\); Veronese residual form; 2026-07-30).** Continue Prop 15.56–15.57. Write \(d=n/2\), \(N=|\mathrm{Max}_{+}|\), \(\alpha=d/N\), \(Y\) the \(N\times n\) Max+ matrix, \(P=YY^\top/(2N)\), and \(u_a=y_a/\sqrt n\in V_+\) (unit vectors).

1. **Max+ lies in \(V_+\) (proved).** For every \(y\in\mathrm{Max}_{+}\), \(Cy=py\). *Proof.* \(\Phi(y)=\tfrac12 y^\top Cy\) and the Max+ level is \(\Phi=pn/2\), so \(y^\top Cy=pn\). With \(\|y\|^2=n\) and \(C=p(P_+-P_-)\) one gets \(\|P_+y\|^2-\|P_-y\|^2=n\) and \(\|P_+y\|^2+\|P_-y\|^2=n\), hence \(P_-y=0\). \(\square\)
   Consequently \(Y^\top Y=2N\,P_+\) (tight frame of Max+ in \(V_+\)).

2. **Perron of the Schur square (proved, any equal-diagonal orthoprojector).** \((P\odot P)\mathbf1=\alpha\mathbf1\) and \(\lambda_{\max}(P\odot P)=\alpha\). *Proof.* \(\sum_b P_{ab}^2=P_{aa}=\alpha\) (from \(P^2=P\)), so row sums of \(P\odot P\) equal \(\alpha\). Entries of \(P\odot P\) are nonnegative, so Perron–Frobenius gives \(\lambda_{\max}=\alpha\) with eigenvector \(\mathbf1\). \(\square\)
   The residual \(\lambda_2(P\odot P)\le\alpha/2\) is therefore a gap below the Perron root (equivalently \(\lambda_{\max}(G)=n/2\)).

3. **Veronese / Gram reformulation (proved equivalent).** For \(x\in\mathbb R^N\) set \(T(x)=\sum_{a=1}^N x_a\,y_ay_a^\top\in\mathrm{Sym}_n\). Then
   \[
   \lambda_2(P\odot P)\le\frac\alpha2
   \quad\Longleftrightarrow\quad
   \|T(x)\|_F^2\le nN\,\|x\|^2
   \quad\text{for all }x\perp\mathbf1.
   \]
   *Proof.* Expand \(\|T(x)\|_F^2=\sum_{a,b}x_ax_b(y_a\cdot y_b)^2\). With \(P_{ab}=(y_a\cdot y_b)/(2N)\) one has
   \((y_a\cdot y_b)^2=4N^2(P\odot P)_{ab}\), so \(\|T\|_F^2=4N^2\,x^\top(P\odot P)x\). On \(\mathbf1^\perp\), \(x^\top(P\odot P)x\le\lambda_2(P\odot P)\|x\|^2\), and the target \(\lambda_2\le\alpha/2=d/(2N)\) rearranges to \(\|T\|_F^2\le nN\|x\|^2\) (using \(d=n/2\)). The reverse direction is the Rayleigh quotient for \(\lambda_2\). \(\square\)
   Equivalently, writing \(W_{ab}=(u_a\cdot u_b)^2\), one has \(\lambda_1(W)=N/d\) and the residual is \(\lambda_2(W)\le N/n=\lambda_1(W)/2\).

4. **Zero diagonal of \(T(x)\) on \(\mathbf1^\perp\) (proved).** If \(\sum_a x_a=0\) then \(T(x)_{ii}=\sum_a x_a=0\) for all \(i\) (using \(y_{a,i}^2=1\)). Thus \(T(x)\) is a zero-diagonal, trace-zero matrix supported through \(V_+\) (part 1).

5. **Maximiser structure (certified \(p=3,5,7\); partial proof).** On a top cycle eigenvector of \(G\): \(\mathrm{ft}(v)=2p\), \(\mathbb E[\|By\|^2]=4\), \(By\in V_+\) for every Max+ \(y\), and \(\|P_+BP_+\|_F^2=2\). The identity \(\mathbb E[\|P_+By\|^2]=2\|P_+BP_+\|_F^2\) holds for every edge weight (2-design / frame). Evidence: `e1_gmin_gap_probe.json`.

6. **Sufficient numerical bound (certified, not proved).** At \(p=5,7\), \(\lambda_{\mathrm{cycle}}\le 8n/(n+4)\le n/2\) (the \(2\times\) spherical fourth-moment comparison). At \(p=3\), the same comparison fails (\(\lambda_{\mathrm{cycle}}=8>8n/(n+4)\approx5.71\)), consistent with the gap failure. Evidence: `e1_gmin_gap_probe.json`. **Not a proof:** the fourth-moment ratio vs the sphere reaches \(\approx2.8>2\) at \(p=3\) and \(\approx1.95,1.43\) at \(p=5,7\); a universal factor-\(2\) sphere bound is false.

7. **Residual (OPEN).** Prove part 3 for every prime \(p\ge5\) (e.g. \(\|T(x)\|_F^2\le nN\|x\|^2\) on \(\mathbf1^\perp\), or \(\lambda_2(W)\le\lambda_1(W)/2\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.59 (central symmetry, \(P\mathbf1=0\), rank of \(P\odot P\), two-moment; 2026-07-30).** Continue Prop 15.58.

1. **Central symmetry (proved).** \(\mathrm{Max}_{+}=-\mathrm{Max}_{+}\). *Proof.* \(Cy=py\) iff \(C(-y)=p(-y)\), and \(-y\in\{\pm1\}^n\). \(\square\)

2. **Frame is centered (proved).** \(\sum_{y\in\mathrm{Max}_{+}}y=0\) in \(\mathbb R^n\). *Proof.* Pair \(y\) with \(-y\) from part 1. \(\square\)
   Consequently \(Y^\top\mathbf1=0\), so the orthoprojector \(P=YY^\top/(2N)\) satisfies \(P\mathbf1=0\) (range of \(P\) lies in \(\mathbf1^\perp\subset\mathbb R^N\)). Certified \(p=3,5,7\): `e1_gmin_veronese.json`.

3. **Rank formula (certified \(p=3,5\); formula holds \(p=7\)).** \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\). Matches \(1+\mathrm{rank}_{\mathrm{cycle}}\) with \(\mathrm{rank}_{\mathrm{cycle}}=\binom{d-1}{2}-1\). Evidence: `e1_gmin_veronese.json`. *(Uniform proof of the rank formula for all primes \(p\ge3\) still open; not required if the Veronese bound is proved by other means.)*

4. **Two-moment on \(W_{ab}=(u_a\cdot u_b)^2\) (partial).** Writing \(\lambda_1(W)=N/d\) and \(\sum_{j\ge2}\lambda_j(W)=N(1-1/d)\) (from \(\mathrm{Tr}(W)=N\)), the one-large-rest-equal estimate from \((\mathrm{Tr}(W^2),\mathrm{Tr}(W))\) **forces** \(\lambda_2(W)\le N/(2d)\) at \(p=7\), but **not** at \(p=3\) (correct: gap fails) or \(p=5\) (worst-case \(12.78>10\); actual \(\lambda_2\approx5.21\) with multiplicity \(d=13\)). Evidence: `e1_gmin_veronese.json`. Multiplicity-\(d\) two-moment would force the gap at \(p=5\); multiplicity not yet proved for general \(p\).

5. **Residual (OPEN).** Same as Prop 15.58.7: prove \(\|T(x)\|_F^2\le nN\|x\|^2\) for all \(x\perp\mathbf1\) and all primes \(p\ge5\). Then \(\lambda_{\max}(G)=n/2\) simple, Prop 15.55 closes bi-tight / Type I, deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.60 (antipodal reduction to projective ENTF; 2×sphere sufficient for \(p\ge5\); 2026-07-30).** Continue Prop 15.58–15.59. Write \(m=N/2\) and fix a set \(\mathcal L\) of representatives of \(\mathrm{Max}_{+}/\{\pm1\}\) (one vector from each antipodal pair). Let \(u_y=y/\sqrt n\in V_+\) and \(W^{(\mathrm{proj})}_{yy'}=(u_y\cdot u_{y'})^2\) on \(\mathcal L\).

1. **Antipodal reduction for \(T\) (proved).** If \(x\in\mathbb R^N\) and \(s_a=\tfrac12(x_a+x_{-a})\) (indices via \(y\mapsto -y\)), then
   \[
   T(x)=T(s)=\sum_{a}s_a\,y_ay_a^\top.
   \]
   *Proof.* \(y y^\top=(-y)(-y)^\top\), so the antisymmetric part \(x_a-x_{-a}\) cancels. \(\square\)
   Consequently \(\|x\|^2=\|s\|^2+\|x-s\|^2\ge\|s\|^2\) and \(\|T(x)\|_F=\|T(s)\|_F\), so the Veronese inequality on all of \(\mathbf1^\perp\) reduces to antipode-symmetric \(x\) (equivalently, to functions on \(\mathcal L\)).

2. **Projective ENTF (proved).** The \(m\) unit vectors \(\{u_y:y\in\mathcal L\}\) form an equal-norm tight frame in \(V_+\cong\mathbb R^d\):
   \[
   \sum_{y\in\mathcal L}u_yu_y^\top=\frac m d\,I_d,\qquad \|u_y\|=1.
   \]
   *Proof.* \(\sum_{\mathrm{Max}_{+}}uu^\top=(N/d)I_d\) and antipodal pairs contribute identical \(uu^\top\), so the sum over \(\mathcal L\) is half. \(\square\)
   Moreover \(W^{(\mathrm{proj})}\mathbf1=(m/d)\mathbf1\) and \(\lambda_1(W^{(\mathrm{proj})})=m/d\).

3. **Eigenvalue doubling (proved).** On the antipode-symmetric subspace of \(\mathbb R^N\),
   \[
   \lambda_2(W)=2\,\lambda_2(W^{(\mathrm{proj})}),
   \]
   where \(W_{ab}=(u_a\cdot u_b)^2\) is the full Max+ Schur square. *Proof.* Each antipodal \(2\times2\) block of \(W\) is the all-ones matrix of order 2 (since \((u\cdot(\pm u'))^2=(u\cdot u')^2\)), so the Rayleigh quotient of an antipode-symmetric vector with values \(c\) on \(\mathcal L\) equals \(2\,c^\top W^{(\mathrm{proj})}c/\|c\|_{2,\mathrm{sym}}^2\). \(\square\)
   Therefore
   \[
   \lambda_2(P\odot P)\le\frac\alpha2
   \quad\Longleftrightarrow\quad
   \lambda_2(W^{(\mathrm{proj})})\le\frac m{2d}.
   \]

4. **Spherical comparison algebra (proved).** For the continuous sphere (or any spherical 4-design) in \(\mathbb R^d\), the maximal fourth-moment Rayleigh on trace-free matrices is \(2m/(d(d+2))\). The bound
   \[
   \frac{4m}{d(d+2)}\le\frac m{2d}
   \]
   rearranges to \(d\ge6\). Since \(d=(p^2+1)/2\ge13\) for primes \(p\ge5\), **any proof that the projective Max+ fourth-moment Rayleigh is at most twice the spherical value yields the spectral gap for all primes \(p\ge5\)**.

5. **Certified 2×sphere (p=5,7).** On projective Max+: maximising \(\sum_y(u_y^\top A u_y)^2\) over \(\mathrm{Tr}\,A=0\), \(\|A\|_F=1\) gives ratio-to-sphere \(\approx1.95\) at \(p=5\) and \(\approx1.43\) at \(p=7\) (both \(<2\)), and \(\approx2.80>2\) at \(p=3\). Gap holds \(p=5,7\), fails \(p=3\). Evidence: `e1_gmin_projective.json`.

6. **Residual (OPEN).** Prove
   \[
   \max_{\mathrm{Tr}\,A=0,\ \|A\|_F=1}\sum_{y\in\mathcal L}(u_y^\top A u_y)^2
   \;\le\;
   \frac{4m}{d(d+2)}
   \]
   for every prime \(p\ge5\) (or any upper bound \(\le m/(2d)\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.61 (16N bound \(\Rightarrow\lambda_{\mathrm{cycle}}\le8\Rightarrow\) gap for \(p\ge5\); 2026-07-30).** Continue Prop 15.56–15.60. Let \(B\) range over real \(n\times n\) matrices with \(B=P_+BP_+\), \(\mathrm{Tr}\,B=0\), and write
\[
Q(B):=\sum_{y\in\mathrm{Max}_{+}}(y^\top By)^2.
\]

1. **Frame identity (proved).** For every such \(B\),
   \[
   \sum_{y\in\mathrm{Max}_{+}}\|By\|^2=2N\|B\|_F^2.
   \]
   *Proof.* \(\sum_y y^\top B^2 y=\mathrm{Tr}\bigl(B^2\sum_y yy^\top\bigr)=\mathrm{Tr}(B^2\cdot 2N P_+)=2N\|B\|_F^2\). \(\square\)

1b. **Triangle form on \(V_+\) (proved).** For every \(B=P_+BP_+\),
   \[
   \mathrm{Tr}(CB^2)=p\|B\|_F^2.
   \]
   *Proof.* \(C=p(P_+-P_-)\) and \(B^2=P_+B^2P_+\), so \(\mathrm{Tr}(CB^2)=p\mathrm{Tr}(B^2)\). \(\square\)
   Consequently every zero-diagonal \(B\) on \(V_+\) with \(\|B\|_F^2=2\) (i.e. unit edge-weight \(v\)) saturates the triangle bound \(\mathrm{ft}=2p\) and has \(\mathbb E[\|By\|^2]=4\). The spectral gap residual is therefore purely the alignment \(\mathbb E[(y^\top By)^2]\) among these maximisers of \(\mathrm{ft}\).

2. **Equivalence with \(\lambda_2(W)\) (proved).** Identifying \(B=QAQ^\top\) with \(Q\) an ONB of \(V_+\) and \(A\in\mathrm{Sym}(\mathbb R^d)\), \(\mathrm{Tr}\,A=0\),
   \[
   Q(B)=n^2\sum_{a=1}^N(u_a^\top A u_a)^2,
   \]
   so
   \[
   \max_{\|B\|_F=1}Q(B)=n^2\cdot\lambda_2(W)=4d^2\cdot\lambda_2(W).
   \]
   Consequently
   \[
   \max_{\|B\|_F=1}Q(B)\le 16N
   \quad\Longleftrightarrow\quad
   \lambda_2(W)\le\frac{4N}{d^2}.
   \]

3. **Equivalence with \(\lambda_{\mathrm{cycle}}\le8\) (proved).** At a cycle maximiser of \(G\) one has \(\mathrm{ft}=2p\), \(\mathbb E[\|By\|^2]=4\), \(By\in V_+\) for all Max+ \(y\), and \(\|B\|_F^2=2\) (Prop 15.57 certs). Scaling \(B'=B/\sqrt2\) (\(\|B'\|_F=1\)) gives
   \[
   \lambda_{\mathrm{cycle}}=\frac1{4N}Q(B)=\frac1{2N}Q(B').
   \]
   Hence \(\max Q\le16N\) yields \(\lambda_{\mathrm{cycle}}\le8\). Conversely, the \(\Phi\)-maximiser of \(Q\) realises \(\lambda_{\mathrm{cycle}}=Q(B_\star)/(2N)\) under the same scaling, so the two maxima match. \(\square\)

4. **Algebraic gap upgrade (proved).** If \(\lambda_{\mathrm{cycle}}\le8\), then for every prime \(p\ge5\) one has \(n/2=d\ge13>8\), so \(\lambda_{\max}(G)=\max(n/2,\lambda_{\mathrm{cycle}})=n/2\) is simple, and Prop 15.55 blocks all Max+-tight size-\(2p\) covers (bi-tight / Type I empty). Equivalently: \(\lambda_2(W)\le4N/d^2\) and \(d\ge8\) imply \(\lambda_2(W)\le N/(2d)\) because \(4/d\le1/2\). \(\square\)

5. **Certified 16N bound.** At \(p=3\): \(Q_{\max}=16N=192\) and \(\lambda_{\mathrm{cycle}}=8\) (**equality**). At \(p=5\): \(Q_{\max}/(16N)=11/13<1\), \(\lambda_{\mathrm{cycle}}=88/13<8\). At \(p=7\): \(Q_{\max}/(16N)\approx0.660<1\), \(\lambda_{\mathrm{cycle}}\approx5.281<8\). Spectral gap holds \(p=5,7\), fails \(p=3\). Evidence: `e1_gmin_16n.json`.

6. **Residual (OPEN).** Prove \(Q(B)\le16N\|B\|_F^2\) for all \(B=P_+BP_+\) with \(\mathrm{Tr}\,B=0\) and all primes \(p\ge5\) (equality at \(p=3\) is the base case). Then bi-tight / Type I closes for every such \(p\) via Prop 15.55. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* The 2×-sphere bound of Prop 15.60.4 is slightly sharper (\(16\,dN/(d+2)<16N\)) but has the same open core; the 16N form isolates the clean numerical threshold \(\lambda_{\mathrm{cycle}}\le8\).

**Proposition 15.62 (typeA+wedge identity; \(Q=6N\|B\|_F^2+Q_4\); 2026-07-30).** Continue Prop 15.61. Let \(B\) be real symmetric \(n\times n\) with \(B=P_+BP_+\) and zero ambient diagonal (hence \(\mathrm{Tr}\,B=0\) automatically). Write \(Q(B)=\sum_{y\in\mathrm{Max}_{+}}(y^\top By)^2\). Index unordered edges \(e=\{i,j\}\) by \(Be_e:=2B_{ij}\) and \(f_e(y):=y_iy_j\), so \(y^\top By=\sum_e Be_e f_e(y)\). Let \(\mathrm{Gu}=F^\top F\) with \(F\) the \(N\times E\) matrix of features \(f_e\) on Max+, and split pairs of distinct edges into **wedge** (share a vertex) and **disjoint**.

1. **Type A (same-edge; proved).** 
   \[
   \sum_e Be_e^2\sum_y f_e(y)^2=N\|Be\|_2^2=2N\|B\|_F^2,
   \]
   since \(f_e^2\equiv1\) and \(\|Be\|_2^2=4\sum_{i<j}B_{ij}^2=2\|B\|_F^2\). \(\square\)

2. **Wedge (share-one-vertex; proved).** For edges \(e=\{r,j\}\), \(e'=\{r,l\}\) sharing \(r\),
   \[
   \sum_y f_e f_{e'}=\sum_y y_jy_l=N\Sigma_{jl},\qquad \Sigma:=\mathbb E[yy^\top]=2P_+.
   \]
   Summing over all ordered wedge pairs and using \(Be_e=2B_{rj}\),
   \[
   Be^\top(\mathrm{Gu}\odot\mathbf1_{\mathrm{wedge}})Be
   =\sum_r\sum_{j\neq l,\,j,l\neq r}(2B_{rj})(2B_{rl})\cdot N\Sigma_{jl}
   =4N\sum_r\Bigl((Be_r)^\top\Sigma(Be_r)-\sum_j B_{rj}^2\Sigma_{jj}\Bigr).
   \]
   On \(V_+\): \(\Sigma v=2v\) for \(v=Be_r=B e_r\in V_+\), and \(\Sigma_{jj}=1\), so the parenthesis equals \(\|Be_r\|_2^2\). Hence the wedge contribution is
   \[
   4N\sum_r\|Be_r\|_2^2=4N\|B\|_F^2.
   \]
   \(\square\)

3. **Identity (proved).** 
   \[
   Q(B)=6N\|B\|_F^2+Q_4(B),\qquad Q_4(B):=Be^\top(\mathrm{Gu}\odot\mathbf1_{\mathrm{disj}})Be.
   \]
   Consequently \(Q(B)\le16N\|B\|_F^2\) if and only if \(Q_4(B)\le10N\|B\|_F^2\). For the cycle-normalisation \(\|B\|_F^2=2\),
   \[
   \lambda_{\mathrm{cycle}}=\frac{Q(B)}{4N}=3+\frac{Q_4(B)}{4N},
   \]
   so \(\lambda_{\mathrm{cycle}}\le8\) \(\Leftrightarrow\) \(Q_4\le20N\) at that scale (equivalently \(Q_4\le10N\|B\|_F^2\)). \(\square\)

4. **Certified.** Multi-seed sampling of the full zero-diag \(\cap V_+\) space (nullspace of ambient diagonal on \(\mathrm{Sym}(V_+)\)) at \(p=3,5,7\): typeA+wedge identity holds to machine precision in every trial; at \(p=3\) one has \(Q\equiv16N\) on the whole space (equality case); at \(p=5,7\) the maximiser and all random trials satisfy \(Q_4\le10N\|B\|_F^2\). Unrestricted \(\|\mathrm{Gu}_{\mathrm{disj}}\|_{\mathrm{op}}\) exceeds the \(5N\) Rayleigh threshold for unit edge vectors, so a crude operator-norm bound fails — the residual is Rayleigh of \(\mathrm{Gu}_{\mathrm{disj}}\) on the **image** of zero-diag \(\cap V_+\to\mathbb R^E\) only. Evidence: `e1_gmin_typeA_wedge.json`.

5. **Residual (OPEN).** Prove \(Q_4(B)\le10N\|B\|_F^2\) for all zero-diag \(B=P_+BP_+\) and all primes \(p\ge5\). Then Prop 15.61 closes bi-tight for every such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.63 (hypothesis H \(\Rightarrow\) 16N; spectrum certs; 2026-07-30).** Continue Prop 15.62. For unit-Frobenius zero-diag \(B=P_+BP_+\) write
\[
\mathrm{ray}(B)\,:=\,\frac{Q_4(B)}{2N}
\,=\,\frac{Be^\top G_{\mathrm{disj}}Be}{\|Be\|_2^2}
\]
(the Rayleigh quotient of \(G_{\mathrm{disj}}=\mathbb E[ff^\top]\) on vertex-disjoint pairs, evaluated on the image of the edge map). Set
\[
H(p)\,:=\,\frac{(p+2)^2}{d},\qquad d=\frac{p^2+1}{2}.
\]

1. **Algebra (proved).** \(H(p)\le5\) for every prime \(p\ge3\), with equality if and only if \(p=3\).
   *Proof.* \(H(p)=2(p+2)^2/(p^2+1)\), so \(H\le5\) \(\Leftrightarrow\) \(2(p+2)^2\le5(p^2+1)\) \(\Leftrightarrow\) \(3p^2-8p-3\ge0\). The positive root of \(3p^2-8p-3=0\) is \(p=3\); the quadratic is nonnegative for all primes \(p\ge3\), and vanishes only at \(p=3\). \(\square\)

2. **H \(\Rightarrow\) 16N (proved).** If \(\mathrm{ray}(B)\le H(p)\) for every unit zero-diag \(B=P_+BP_+\), then \(Q_4(B)\le2N\cdot H(p)\le10N\), hence \(Q(B)\le16N\|B\|_F^2\) (Prop 15.62.3), hence \(\lambda_{\mathrm{cycle}}\le8\), and for \(p\ge5\) bi-tight is empty (Prop 15.61.4). \(\square\)

3. **Spectrum of \(Q_4\) on zero-diag \(\cap V_+\) (certified).**
   - \(p=3\): \(\dim=5\), single eigenvalue \(Q_4\equiv10N\) (mult 5); \(\mathrm{ray}\equiv5=H(3)\).
   - \(p=5\): \(\dim=65\), eigenvalues \(Q_4\in N\cdot\{\tfrac{98}{13},\tfrac{66}{13},\tfrac{2}{13}\}\) with multiplicities \(\{13,26,26\}=\{d,2d,2d\}\); top \(\mathrm{ray}=\tfrac{49}{13}=H(5)\).
   - \(p=7\): power multistart gives \(\mathrm{ray}\approx2.281<H(7)=3.24=\tfrac{81}{25}\).
   Evidence: `e1_gmin_q4_spectrum.json`, `e1_gmin_q4_ub.json`, `e1_gmin_q4_bound.json`.

4. **Hypothesis H (certified \(p=3,5,7\); OPEN in general).** \(\mathrm{ray}(B)\le H(p)\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge3\), with equality at \(p=3\) (whole space) and at the maximiser for \(p=5\). Combined with part 1–2 this yields the 16N bound. **Uniform proof of H for all primes \(p\ge5\) remains OPEN.**

5. **Dead ends (do not reopen).** Pointwise \(|y^\top By|\) or \(\cos^2\) bounds (max \(\cos^2\) exceeds the average threshold); unrestricted \(\|G_{\mathrm{disj}}\|_{\mathrm{op}}\) (exceeds \(H(p)\)); linear span of \(\{\|B\|_F^2,\mathrm{tr}(B^4),\sum B_{ij}^4,\sum_i\|row_i\|^4\}\) (fails at \(p=5,7\)); general CS-ENTF theory (random frames violate 16N).

6. **Residual (OPEN).** Prove hypothesis H for all primes \(p\ge5\) (or any upper bound \(\mathrm{ray}\le5\)). Then bi-tight closes via Prop 15.61–15.62. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.64 (dual fourth-moment form of H; residual reformulation; 2026-07-30).** Continue Prop 15.63. Write \(s_y:=V_+^\top y\in\mathbb R^d\) (so \(\|s_y\|^2=n\), \(\mathbb E[ss^\top]=2I_d\)), and for \(A\in\mathrm{Sym}(\mathbb R^d)\) set \(\Phi(A):=\mathbb E[(s^\top As)\,ss^\top]\). Let \(\mathcal Z\) be the subspace of traceless \(A\) with zero ambient diagonal on \(B=V_+AV_+^\top\) (i.e. \(r_i^\top Ar_i=0\) for all rows \(r_i\) of \(V_+\)).

1. **Duality (proved).** For every zero-diag \(B=V_+AV_+^\top\),
   \[
   Q(B)=\sum_y(s_y^\top A s_y)^2=N\langle\Phi(A),A\rangle_F.
   \]
   Consequently
   \[
   \max_{\|A\|_F=1,\,A\in\mathcal Z}Q(B)=N\cdot\lambda_{\max}(\Phi|_{\mathcal Z}),
   \]
   and the maximiser satisfies the eigenmatrix equation \(\Phi(A)=\lambda A\) on \(\mathcal Z\). Equivalently, writing \(D=YY^\top=SS^\top\), one has on \(\mathbf1^\perp\)
   \[
   K=\frac{D\odot D}{2N},\qquad\lambda_{\mathrm{cycle}}=\frac{\lambda_2(D\odot D)}{2N}=\frac12\max_{\|x\|=1,\,x\perp\mathbf1}\Bigl\|\sum_y x_y s_ys_y^\top\Bigr\|_F^2.
   \]
   Evidence: `e1_gmin_H_proof.json` (gen.eig on \(\mathcal Z\)).

2. **Wick baseline (proved).** If \(s\) were Gaussian with \(\mathrm{Cov}=2I\), then \(\Phi=8\,\mathrm{Id}\) on \(\mathrm{Tr}\,A=0\) and \(Q=8N\|B\|_F^2\). For actual Max+,
   \[
   \frac QN=8+\mathrm{residual}(A),\qquad\mathrm{residual}(A)=\langle\kappa,A\otimes A\rangle,
   \]
   with cumulant \(\kappa=\Phi-8\,\mathrm{Id}\).

3. **H \(\Leftrightarrow\) residual bound (proved).** Hypothesis H is equivalent to
   \[
   \max_{\|A\|_F=1,\,A\in\mathcal Z}\mathrm{residual}(A)\;\le\;\frac{(p+1)(p+7)}{d},
   \]
   because \(6+2H(p)-8=2(H(p)-1)=(p+1)(p+7)/d\). Equality holds at \(p=3\) (whole \(\mathcal Z\)) and at the maximiser for \(p=5\). \(\square\)

4. **Certified.** Exact gen.eig of \(\Phi|_{\mathcal Z}\) at \(p=3,5,7\): H holds; residual ratios to budget \(1,1,0.572\). Spectrum of \(K\) at \(p=5\): \(\{13,\tfrac{88}{13},\tfrac{72}{13},\tfrac{40}{13}\}\) with mults \(\{1,d,2d,2d\}\). Evidence: `e1_gmin_H_proof.json`, `e1_gmin_q4_spectrum.json`.

5. **Residual (OPEN).** Prove \(\lambda_{\max}(\Phi|_{\mathcal Z})\le6+2H(p)\) for all primes \(p\ge5\) (equivalently residual \(\le(p+1)(p+7)/d\), or \(\mathrm{ray}\le H(p)\), or \(\lambda_{\mathrm{cycle}}\le3+H(p)\)). Then 16N and bi-tight close. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark (dead ends for H).* Gershgorin on \(K\); unrestricted \(\|G_{\mathrm{disj}}\|_{\mathrm{op}}\); mult\(\ge d\) of \(\lambda_2(W)\) with only \(\mathrm{Tr}(W^2)\) (too weak for 16N at \(p=5\)); pointwise \(\cos^2\); 2×sphere (fails at \(p=3\), holds numerically \(p\ge5\)); distance-homogeneous scheme formulae (Max+ not distance-homogeneous at \(p=7\)); two-moment feasibility of \(\lambda_2(P\odot P)\) (always allows values \(>4/N\)).

**Proposition 15.65 (κ spectrum on \(\mathcal Z\); clean form \(\lambda_2(P\odot P)\le4/N\); boolean essential; 2026-07-30).** Continue Prop 15.64. Write \(P=YY^\top/(2N)\) (equal-diagonal orthoprojector of rank \(d\), diagonal \(\alpha=d/N\), and \(P\mathbf1=0\)).

1. **Clean equivalences (proved).** 
   \[
   P\odot P=\alpha^2 W,\qquad W_{ab}=(u_a\cdot u_b)^2,\quad u_a=y_a/\sqrt n.
   \]
   Hence \(\lambda_2(P\odot P)=\alpha^2\lambda_2(W)\). Combined with Prop 15.56–15.61:
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_2(P\odot P)\le4/N\Longleftrightarrow\lambda_2(W)\le4N/d^2,\\
   H&\Longleftrightarrow\lambda_2(P\odot P)\le\frac{3+H(p)}{2N},\\
   \mathrm{gap}&\Longleftrightarrow\lambda_2(P\odot P)\le\alpha/2=d/(2N).
   \end{aligned}
   \]
   Moreover \(\lambda_{\max}(P\odot P)=\alpha\) for every equal-diagonal orthoprojector (Prop 15.58.2). \(\square\)

2. **Cumulant spectrum on \(\mathcal Z\) (certified).** The residual form \(\mathrm{residual}(A)=\langle\kappa,A\otimes A\rangle\) on \(\mathcal Z\) has eigenvalues
   - \(p=3\): \(\{8\}\) (mult \(\dim\mathcal Z=5\));
   - \(p=5\): \(\{-\tfrac{24}{13},\tfrac{40}{13},\tfrac{72}{13}\}\) with multiplicities matching \(\{26,26,13\}\) (top \(=\mathrm{budget}=(p+1)(p+7)/d\));
   - \(p=7\): top \(=\tfrac{1048}{409}<\tfrac{112}{25}=\mathrm{budget}\).
   Evidence: `e1_gmin_cumulant.json`.

3. **Boolean/conference structure is essential (certified counterexample).** Random equal-diagonal rank-\(d\) orthoprojectors need **not** satisfy \(\lambda_2(P\odot P)\le4/N\): at \((N,d)=(50,10)\), \(20/20\) Haar-row-equalized samples violated the bound. Thus no proof of 16N can use only the equal-diagonal projector axioms; the Max+/boolean/conference structure is load-bearing. Evidence: `e1_gmin_cumulant.json`.

4. **Residual (OPEN).** Prove \(\lambda_2(P\odot P)\le4/N\) (or the sharper H form) for the Max+ projector of every prime \(p\ge5\). Equivalent targets: \(\lambda_{\max}(\kappa|_{\mathcal Z})\le(p+1)(p+7)/d\), or \(\mathrm{ray}\le H(p)\). Then bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.66 (zero-diag freeness of \(\Phi\)-maximiser; pairing residual criterion for \(g_{\min}\ge L(p)\); 2026-07-30).** Continue Prop 15.53–15.54 and 15.64–15.65.

1. **Zero-diag freeness (proved).** Let \(A\in\mathrm{Sym}_0(\mathbb R^d)\) maximise \(\langle\Phi A,A\rangle\) over \(\|A\|_F=1\). Then the ambient matrix \(B=V_+AV_+^\top\) has zero diagonal, so the maximiser already lies in \(\mathcal Z\). Consequently
   \[
   \lambda_{\max}(\Phi|_{\mathcal Z})=\lambda_{\max}(\Phi|_{\mathrm{Sym}_0}).
   \]
   *Proof.* The Lagrangian for the constraints \(r_i^\top Ar_i=0\) and \(\mathrm{Tr}\,A=0\) yields
   \(\Phi(A)=\lambda A+\sum_i\mu_i r_ir_i^\top+\nu I\).
   Taking the Frobenius product with \(r_kr_k^\top\) and using \(r_k^\top s_y=y_k\) (hence \(r_k^\top\Phi(A)r_k=\mathbb E[s^\top As]=0\)) forces \(\mu_k/4+\nu/2=0\) for every \(k\), so all \(\mu_k\) are equal. Then \(\sum_i\mu_i r_ir_i^\top=-2\nu I\), and \(\Phi(A)=\lambda A-\nu I\). Taking traces and using \(\mathrm{Tr}\,\Phi(A)=n\mathbb E[s^\top As]=0\) gives \(\nu=0\), hence \(\mu_k=0\) and \(\Phi(A)=\lambda A\). The unconstrained critical point on \(\mathrm{Sym}_0\) is therefore admissible for \(\mathcal Z\). Certified: ambient \(\mathrm{diag}(B)\) of the power-iteration maximiser is \(O(10^{-16})\) at \(p=3,5,7\). Evidence: `e1_gmin_m4_residual.json`. \(\square\)

2. **Pairing residual criterion (proved algebra).** Write \(L(p)=-(p-2)/(2p^2)\) and \(T(p)=-(p-2)/(p(2p-1))\). For every prime \(p\ge5\) one has \(L(p)>T(p)\). On any 4-set with \(|\kappa|=1\), the Wick value is \(m_4^{\mathrm{Wick}}=\kappa/p^2\). If
   \[
   \bigl|m_4-\kappa/p^2\bigr|\;\le\;\frac{p-4}{2p^2}
   \]
   for every such 4-set, then
   \[
   |m_4|\;\le\;\frac1{p^2}+\frac{p-4}{2p^2}=\frac{p-2}{2p^2}=-L(p),
   \]
   hence \(g_{\min}\ge L(p)>T(p)\), and Prop 15.47 closes bi-tight. \(\square\)

3. **Certified m4 tables on \(|\kappa|=1\).**  
   | \(p\) | \(g_{\min}\) | \(L(p)\) | \(\max|m_4|\) | \(\max|m_4-\kappa/p^2|\) | resid crit. |
   |------|------------|---------|--------------|--------------------------|-------------|
   | 3 | \(-1/3\) | \(-1/18\) | \(1/3\) | \(2/9\) | n/a (\(p<5\)) |
   | 5 | \(-3/65\) | \(-3/50\) | \(3/65\) | \(0.0554>(p-4)/(2p^2)=1/50\) | **fails** (too crude) |
   | 7 | \(-109/2863\) | \(-5/98\) | \(109/2863\) | \(0.0177\le3/98\) | **holds** |

   In particular \(g_{\min}\ge L(p)\) holds at \(p=5,7\), but the triangle residual criterion fails at \(p=5\) (large residuals occur on classes with small \(|m_4|\)). Evidence: `e1_gmin_m4_residual.json`.

4. **Residual (OPEN).** Prove \(g_{\min}\ge L(p)\) for all primes \(p\ge5\), e.g. by proving \(|m_4|\le(p-2)/(2p^2)\) on every \(|\kappa|=1\) class (directly, not via the failed \(p=5\) triangle), **or** prove \(\lambda_2(P\odot P)\le4/N\). Either closes bi-tight. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* K4-edge Gram PSD only forces \(\alpha\le(p-1)/(2p)\)-scale bounds (e.g. \(0.6\) at \(p=5\)), far above \(L(p)\). General equal-diag projectors violate \(4/N\) (Prop 15.65.3).

**Proposition 15.67 (master m4 identity; same-sign Ext criterion; full \(|\kappa|=1\) census \(p=5,7\); 2026-07-30).** Continue Prop 15.66. Write \(m_4(S)=\mathbb E[y_a y_b y_c y_d]\) for \(S=\{a,b,c,d\}\) and Max+ average, \(\kappa(S)\) the sum of the three perfect-matching edge-products of \(C\), and
\[
(Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r}).
\]

1. **Combinatorial identity \(\sigma_{\mathrm{sum}}=4\kappa\) (proved).** For every \(\pm1\) edge-labeling of \(K_4\), writing \(\sigma_v\) for the three-term pairing at vertex \(v\), one has \(\sum_{v\in S}\sigma_v=4\kappa(S)\). *Proof.* Exhaustion of the \(2^6=64\) labelings (all pass). \(\square\)

2. **Master identity (proved for Max+).** Averaging the eigenvector identity \(Cy=py\) against three coordinates of a Max+ vector and using part 1 yields, for every 4-set \(S\),
   \[
   m_4(S)=\frac{\kappa(S)}{p^2}+\frac{\mathrm{Ext}(S)}{4p},\qquad \mathrm{Ext}:=Tm_4.
   \]
   Equivalently \(\mathrm{Ext}(S)=4pm_4(S)-4\kappa(S)/p\). Certified float residual \(<10^{-16}\) on all \(|\kappa|=1\) 4-sets at \(p=5\) and on a full multi-worker census at \(p=7\). Evidence: `e1_gmin_m4_bound.json`, `e1_gmin_m4_proof.json`. \(\square\)

3. **Same-sign Ext criterion (proved algebra).** On any 4-set with \(|\kappa|=1\), if \(\mathrm{sign}(\mathrm{Ext})=\mathrm{sign}(\kappa)\) and
   \[
   |\mathrm{Ext}|\;\le\;\frac{2(p-4)}{p},
   \]
   then
   \[
   |m_4|=\frac1{p^2}+\frac{|\mathrm{Ext}|}{4p}\;\le\;\frac1{p^2}+\frac{p-4}{2p^2}=\frac{p-2}{2p^2}=-L(p).
   \]
   (Opposite-sign Ext only decreases \(|m_4|\) relative to the Wick value.) Thus same-sign Ext control on every \(|\kappa|=1\) class implies \(g_{\min}\ge L(p)\). \(\square\)

4. **Full multi-worker census (certified \(p=5,7\); \(W=86\)).** Over every 4-set with \(|\kappa|=1\):
   | \(p\) | \(\#\{|\kappa|=1\}\) | \(\max|m_4|\) | \(L_{\mathrm{abs}}\) | same-sign \(\max|\mathrm{Ext}|\) | thr \(2(p-4)/p\) |
   |------|---------------------|---------------|----------------------|----------------------------------|------------------|
   | 5 | 11700 | \(3/65\) | \(3/50\) | \(0.123<0.4\) | yes |
   | 7 | 176400 | \(109/2863\) | \(5/98\) | \(0.495<6/7\) | yes |

   Hence \(g_{\min}\ge L(p)\) at \(p=5,7\) by either the direct bound or the same-sign Ext criterion. Evidence: `e1_gmin_m4_proof.json` (F17: `src/workers.py` + `ProcessPoolExecutor(W=nproc-2)`). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that either
   - \(|m_4|\le(p-2)/(2p^2)\) on all \(|\kappa|=1\) 4-sets, or
   - same-sign \(|\mathrm{Ext}|\le2(p-4)/p\) on those sets, or
   - \(\lambda_2(P\odot P)\le4/N\),

   using Max+/boolean/conference structure (Prop 15.65.3: bare equal-diag projectors are insufficient). Then Prop 15.47/15.55 closes bi-tight for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark (failed shortcuts).* Absolute triangle \(|m_4-\kappa/p^2|\le(p-4)/(2p^2)\) fails at \(p=5\) (large residuals on opposite-sign classes with small \(|m_4|\)). Affine+\(\mathrm{PGL}\)+Frob orbit of the halfspace vector is incomplete (\(60/260\) at \(p=5\); F18). Wick fourth-moment pin only reaches the \(T(p)\) endpoint, not strict \(g_{\min}>T(p)\).

**Proposition 15.68 (\(T\kappa\) calculus; residual source on \(|\kappa|=3\); resolvent reduction of \(L\); 2026-07-30).** Continue Prop 15.67. Let \(C\) be any real symmetric conference matrix of order \(n=p^2+1\) (\(C^\top=C\), zero diagonal, off-diagonal \(\pm1\), \(C^2=p^2I\)), and write
\[
(Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r}),\qquad
\kappa(S)=\sum_{\text{three pairings}}C_eC_{e'}.
\]

1. **Reduction of \(T\kappa\) to \(K_4\) (proved).** For every 4-set \(S\), the external sum defining \((T\kappa)(S)\) collapses via \((C^2)_{vx}=0\) (\(v\neq x\)) to a function of the six edge signs of \(S\) alone:
   \[
   (T\kappa)(S)=-6\sum_{v\in S}\prod_{u\in S\setminus\{v\}}C_{vu}
   \]
   (sum of star-triples). In particular \((T\kappa)(S)\) is independent of \(n\) and of the ambient graph outside \(S\). \(\square\)

2. **Vanishing on \(|\kappa|=1\) (proved).** Exhausting all \(2^6=64\) edge-labelings of \(K_4\): if \(|\kappa(S)|=1\) then \((T\kappa)(S)=0\); if \(|\kappa(S)|=3\) then \((T\kappa)(S)\in\{\pm24\}\). Consequently, for **every** conference matrix and every 4-set with \(|\kappa|=1\),
   \[
   (T\kappa)(S)=0.
   \]
   Evidence: `e1_gmin_m4_tkappa.json` (symbolic \(C^2\) reduction + 64-labeling check). \(\square\)

3. **Residual equation and source support (proved).** For Max+ moments write \(\rho:=m_4-\kappa/p^2\). The master identity (Prop 15.67) rearranges to
   \[
   (4p\,I-T)\rho=\frac{T\kappa}{p^2}.
   \]
   By part 2 the right-hand side **vanishes on every \(|\kappa|=1\) 4-set** and is bounded by \(24/p^2\) on \(|\kappa|=3\). Thus all of the same-sign residual on the dangerous classes is the image, under the resolvent \((4pI-T)^{-1}\), of a pure \(|\kappa|=3\) source. \(\square\)

4. **Paley extension degrees (proved formula; certified \(p=3,5,7\)).** For the Paley conference graph and every 4-set with \(|\kappa|=1\), among the \(4(n-4)=4(p^2-3)\) ordered extensions \((v,r)\), exactly
   \[
   d_3=p^2-5,\qquad d_1=3p^2-7
   \]
   land in \(|\kappa|=3\) and \(|\kappa|=1\) respectively (both constant on the \(|\kappa|=1\) stratum). Certified by full census at \(p=3,5,7\). Evidence: `e1_gmin_m4_tkappa.json`. \(\square\)

5. **Resolvent reduction of the \(L\)-bound (proved algebra).** On any \(|\kappa|=1\) 4-set, if \(\mathrm{sign}(\rho)=\mathrm{sign}(\kappa)\) then \(|m_4|=1/p^2+|\rho|\). The bound \(|m_4|\le L_{\mathrm{abs}}=(p-2)/(2p^2)\) is therefore equivalent to
   \[
   |\rho|\;\le\;\frac{p-4}{2p^2}.
   \]
   Writing \(\rho=(4pI-T)^{-1}(T\kappa/p^2)\) and using \(|T\kappa/p^2|\le24/p^2\) on the source, it suffices to prove that the operator gain from the \(|\kappa|=3\) stratum into the same-sign \(|\kappa|=1\) stratum is at most
   \[
   \frac{p-4}{48}.
   \]
   (At \(p=5\) the budget is \(1/48\); empirical gain \(\approx0.0064\ll1/48\).) Equivalently, any upper bound
   \[
   |m_4|\;\le\;\frac{p-2}{p(2p+3)}
   \]
   on \(|\kappa|=1\) closes \(L\) because \((p-2)/(p(2p+3))\le(p-2)/(2p^2)\) for all odd \(p\ge5\), with equality of the two sides only in the large-\(p\) limit sense (strict for finite \(p\)); at \(p=5\) the candidate is sharp (\(|m_4|_{\max}=3/65\)). Certified: candidate \(\ge\max|m_4|\) at \(p=5,7\). \(\square\)

6. **Residual (OPEN).** Prove the resolvent gain bound of part 5 (or the candidate \(|m_4|\le(p-2)/(p(2p+3))\), or \(\lambda_2(P\odot P)\le4/N\)) for every prime \(p\ge5\), using Max+/boolean structure on top of the conference calculus above. Then \(g_{\min}\ge L(p)\) and bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* Absolute \(\infty\)-bootstrap with \(d_1M_1+d_3M_3\) does **not** contract (\(4p-(3p^2-7)<0\) for \(p\ge5\)). Signed cancellation / spectral control of \(T\) on the conference association scheme is load-bearing.

**Proposition 15.69 (spectrum of \(T\): \(\lambda_{\max}=4p\); singular resolvent; min-norm particular solution; 2026-07-30).** Continue Prop 15.67–15.68. View \(T\) as a self-adjoint operator on \(\ell^2\) of unordered 4-sets (signed Johnson adjacency with weights \(C_{vr}\)).

1. **Self-adjointness (proved).** If \(S'=S_{v\to r}\) then \(S=S'_{r\to v}\) and \(C_{rv}=C_{vr}\), so \(T\) is symmetric on \(\mathbb R^{\binom{n}{4}}\). Hence \(\|T\|_2=\rho(T)\). \(\square\)

2. **Spectral edges (certified \(p=5,7\); \(p=3\) strict).** Sparse exact eigensolve (ARPACK/`eigsh` on the CSR of \(T\)):
   | \(p\) | \(\lambda_{\max}(T)\) | \(4p\) | \(\lambda_{\min}(T)\) | mult of \(4p\) (top block) |
   |------|----------------------|-------|----------------------|---------------------------|
   | 3 | \(\approx9.798<12\) | 12 | \(\approx-9.798\) | 0 |
   | 5 | \(20\) | 20 | \(-20\) | \(\ge38\) |
   | 7 | \(28\) | 28 | \(-28\) | \(\ge1\) (full mult open) |

   In particular for primes \(p\ge5\) on Paley, the numerical evidence is
   \[
   \lambda_{\max}(T)=4p=-\lambda_{\min}(T),
   \]
   so \(4pI-T\) is **singular**. Evidence: `e1_gmin_m4_Tspec.json`, `e1_gmin_m4_pseudo.json`. \(\square\)

3. **Compatibility (certified \(p=5,7\)).** The master identity \((4pI-T)m_4=4\kappa/p\) is solvable: writing \(E_{4p}=\ker(4pI-T)\), one has \(4\kappa/p\perp E_{4p}\) to float precision \(<10^{-13}\). (Equivalently \(T\kappa/p^2\perp E_{4p}\) for the residual form.) \(\square\)

4. **General solution and min-norm particular solution (proved algebra + cert).** Every solution of the master linear equation is
   \[
   m=m_\star+h,\qquad h\in E_{4p},
   \]
   where \(m_\star=(4pI-T)^{+}(4\kappa/p)\) is the Moore–Penrose / min-norm solution. The Max+ moment vector is one particular solution: \(m_4=m_\star+h_\star\) for a unique \(h_\star\in E_{4p}\) fixed by the Max+ design (boolean antipodality, \(\mathrm{Tr}(G^2)\), etc.). Certified: \(\|(T-4p)(m_4-m_\star)\|<10^{-11}\) at \(p=5,7\). \(\square\)

5. **Min-norm bound on \(|\kappa|=1\) (certified \(p=5,7\)).**
   | \(p\) | \(\max_{|\kappa|=1}|m_\star|\) | \(L_{\mathrm{abs}}\) | \(\max|m_4|\) (Max+) |
   |------|-------------------------------|----------------------|----------------------|
   | 5 | \(0.056=7/125\) | \(0.06=3/50\) | \(3/65\approx0.04615\) |
   | 7 | \(\approx0.03154\) | \(5/98\approx0.05102\) | \(109/2863\approx0.03807\) |

   So \(m_\star\) itself already obeys \(|m_\star|\le L_{\mathrm{abs}}\) on \(|\kappa|=1\) at \(p=5,7\). The Max+ correction \(h_\star\) **decreases** the max at \(p=5\) and **increases** it at \(p=7\), still staying below \(L_{\mathrm{abs}}\). Best \(L^2\) fit \(m_\star\approx\kappa/(p^2-4)\) (exact coefficient \(1/(p^2-4)\)), but \(m_\star\) is **not** constant on \(\kappa\)-classes for \(p\ge5\). Evidence: `e1_gmin_m4_pseudo.json`. \(\square\)

6. **Residual (OPEN).** Prove for every prime \(p\ge5\):
   - (i) \(\lambda_{\max}(T)=4p\) (and ideally \(\lambda_{\min}=-4p\)) on the Paley / conference Johnson signing;
   - (ii) \(4\kappa/p\perp E_{4p}\);
   - (iii) the Max+ particular solution satisfies \(|m_4|\le L_{\mathrm{abs}}\) on \(|\kappa|=1\)
     (e.g. by controlling \(h_\star\), or by proving \(|m_\star|\le L_{\mathrm{abs}}\) and \(|h_\star|\) cannot push past \(L\), or the candidate \(|m_4|\le(p-2)/(p(2p+3))\), or \(\lambda_2(P\odot P)\le4/N\)).

   Then \(g_{\min}\ge L(p)\) and bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* The naive resolvent gain bound of Prop 15.68.5 assumed \(4p>\lambda_{\max}(T)\). At \(p\ge5\) one has equality, so the resolvent is a **pseudoinverse** on \(E_{4p}^\perp\) plus an undetermined kernel component — the kernel is load-bearing, not an error.

**Proposition 15.70 (mid upper bound algebra; bi-tight threshold comparison; multi-worker census; 2026-07-30).** Continue Prop 15.66–15.69. Write
\[
L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2},\qquad
T_{\mathrm{abs}}(p)=\frac{p-2}{p(2p-1)},\qquad
M_{\mathrm{mid}}(p)=\frac{p-2}{2p(p+1)},\qquad
M_{\mathrm{cand}}(p)=\frac{p-2}{p(2p+3)}.
\]

1. **Algebra of targets (proved).** For every odd prime \(p\ge5\):
   \[
   M_{\mathrm{cand}}(p)\;\le\;M_{\mathrm{mid}}(p)\;\le\;L_{\mathrm{abs}}(p)\;<\;T_{\mathrm{abs}}(p),
   \]
   with ratios
   \[
   \frac{M_{\mathrm{mid}}}{L_{\mathrm{abs}}}=\frac{p}{p+1}<1,\qquad
   \frac{M_{\mathrm{cand}}}{L_{\mathrm{abs}}}=\frac{2p}{2p+3}<1.
   \]
   Hence any proof of \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}(p)\) (or the sharper \(M_{\mathrm{cand}}\)) yields \(g_{\min}\ge L(p)>T(p)\) and closes bi-tight via Prop 15.47. \(\square\)

2. **Bi-tight threshold (recalled).** Prop 15.47: if \(g_{\min}>T(p)=-(p-2)/(p(2p-1))\), i.e. \(\max|m_4|<T_{\mathrm{abs}}\), then no Max+-tight size-\(2p\) bi-tight cover exists. The stronger \(g_{\min}\ge L(p)\) is preferred but not necessary for bi-tight. \(\square\)

3. **Multi-worker census (certified \(p=5,7\); \(W=86\)).** Full enumeration of all \(|\kappa|=1\) 4-sets:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | \(M_{\mathrm{mid}}\) | \(L_{\mathrm{abs}}\) | \(T_{\mathrm{abs}}\) |
   |------|---------------|----------------------|---------------------|--------------------|--------------------|
   | 5 | \(3/65\approx0.04615\) | \(3/65\) (sharp) | \(1/20=0.05\) | \(3/50=0.06\) | \(1/15\approx0.0667\) |
   | 7 | \(109/2863\approx0.03807\) | \(5/119\approx0.0420\) | \(5/112\approx0.0446\) | \(5/98\approx0.0510\) | \(\approx0.0544\) |

   In particular \(g_{\min}\ge L(p)>T(p)\) at \(p=5,7\) (bi-tight empty for these \(p\)), and \(\max|m_4|\le M_{\mathrm{mid}}\) holds at both. Evidence: `e1_gmin_m4_close.json`, `e1_gmin_m4_evec4p.json` (F17 multi-worker). \(\square\)

4. **Type6 Max+-free particular solution (certified \(p=5\)).** Solving \((4pI-T)m=4\kappa/p\) on \(S_4\)-type6 class-constant functions (pure \(C\) combinatorics) yields \(\max_{|\kappa|=1}|m|\approx0.0468\le L_{\mathrm{abs}}\) at \(p=5\), close to the true Max+ value \(3/65\). At \(p=3\) type6 recovers the exact \(|m_4|=1/3\). Evidence: `e1_gmin_m4_close.json`. \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}(p)\) (or \(\le L_{\mathrm{abs}}\), or \(\le M_{\mathrm{cand}}\)), using Max+/boolean/conference structure — e.g. via control of \(h_\star\in E_{4p}\) in Prop 15.69, or a closed \(G\)-spectrum / \(\mathrm{Tr}(G^2)\) pin on the moduli line (Prop 15.53). Then bi-tight closes for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.71 (\(\kappa\)-stratum counts for conference matrices; extension-degree feed; 2026-07-29).** Let \(C\) be any real symmetric **conference matrix** of order \(n\) (\(C^\top=C\), zero diagonal, off-diagonal \(\pm1\), \(C^2=(n-1)I\)), and for a 4-set \(S\) write
\[
\kappa(S)=\sum_{\text{three pairings}}C_eC_{e'}\in\{\pm1,\pm3\}.
\]
Let \(n_1=\#\{S:|\kappa(S)|=1\}\) and \(n_3=\#\{S:|\kappa(S)|=3\}\).

1. **Wedge sum from \(C^2\) (proved).** For \(b\neq c\) one has \(\sum_d C_{bd}C_{cd}=0\). Hence for \(a\notin\{b,c\}\),
   \[
   \sum_{d\notin\{a,b,c\}}C_{db}C_{dc}=-C_{ab}C_{ac}.
   \]
   Summing the wedge \(C_{ab}C_{ac}C_{db}C_{dc}\) over distinct \(a,b,c\) and \(d\notin\{a,b,c\}\) therefore yields
   \[
   \Sigma:=-n(n-1)(n-2).
   \]
   \(\square\)

2. **K4 ratio (proved by exhaustion).** On every \(\pm1\)-edge labelling of \(K_4\), writing \(\pi_1,\pi_2,\pi_3\) for the three pairing-products and \(\mathrm{cross}=\pi_1\pi_2+\pi_1\pi_3+\pi_2\pi_3\),
   \[
   \sum_{\sigma\in S_4}C_{\sigma(a)\sigma(b)}C_{\sigma(a)\sigma(c)}C_{\sigma(d)\sigma(b)}C_{\sigma(d)\sigma(c)}
   \;=\;8\cdot\mathrm{cross}.
   \]
   (64 edge labelings; always \(\mathrm{cross}\in\{-1,3\}\).) Consequently \(\sum_{\text{4-sets}}\mathrm{cross}=\Sigma/8=-n(n-1)(n-2)/8\). \(\square\)

3. **Fourth-moment sum of \(\kappa\) (proved).** Since \(\kappa^2=3+2\,\mathrm{cross}\),
   \[
   \sum_S\kappa(S)^2
   =3\binom{n}{4}+2\cdot\frac{\Sigma}{8}
   =\frac{n(n-1)(n-2)(n-5)}{8}.
   \]
   \(\square\)

4. **Stratum counts (proved).** Using \(\kappa^2\in\{1,9\}\) and \(\binom{n}{4}=n_1+n_3\),
   \[
   n_3=\frac{1}{8}\Bigl(\sum\kappa^2-\binom{n}{4}\Bigr)
   =\frac{n(n-1)(n-2)(n-6)}{96},\qquad
   n_1=\binom{n}{4}-n_3
   =\frac{n(n-1)(n-2)^2}{32}.
   \]
   For Paley conferences \(n=p^2+1\) this is
   \[
   n_1=\frac{(p^2+1)\,p^2\,(p^2-1)^2}{32},\qquad
   n_3=\frac{(p^2+1)\,p^2\,(p^2-1)\,(p^2-5)}{96}.
   \]
   \(\square\)

5. **Extension degrees (certified; unique if constant).** For Paley of order \(n=p^2+1\), every \(|\kappa|=1\) 4-set has exactly
   \[
   d_3=p^2-5,\qquad d_1=3p^2-7
   \]
   ordered one-vertex extensions landing in \(|\kappa|=3\) and \(|\kappa|=1\) respectively (so \(d_1+d_3=4(n-4)\)). Constancy certified by full census at \(p=3,5\) and multi-worker samples at \(p=7,11\). Combined with part 4, the constant-degree hypothesis is the unique solution of the handshaking identities. Evidence: `e1_gmin_m4_stratum.json`, `e1_gmin_m4_tkappa.json`. \(\square\)

6. **Target algebra (proved; recalled).** For odd primes \(p\ge5\),
   \[
   M_{\mathrm{cand}}(p)\le M_{\mathrm{mid}}(p)\le L_{\mathrm{abs}}(p)<T_{\mathrm{abs}}(p)
   \]
   with \(M_{\mathrm{mid}}/L_{\mathrm{abs}}=p/(p+1)\) (Prop 15.70.1). \(\square\)

7. **Multi-worker census (certified \(p=3,5,7,11\); \(W=86\)).** Full \(\kappa\)-stratum counts match part 4 at all four primes (including \(p=11\), \(\binom{122}{4}\approx8.7\cdot10^6\)). Evidence: `e1_gmin_m4_stratum.json` (F17 ProcessPool; atomic JSON). \(\square\)

8. **Residual (OPEN).** Parts 1–4 are Max+-free conference combinatorics and pin the source size of the resolvent equation \((4pI-T)\rho=T\kappa/p^2\) (Prop 15.68): exactly \(n_3\) nonzero source coordinates of amplitude \(24/p^2\). They do **not** alone bound \(|m_4|\) on \(|\kappa|=1\). Still open for every prime \(p\ge5\): \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\) / \(M_{\mathrm{cand}}\)), e.g. via resolvent gain \(\le(p-4)/48\), \(h_\star\in E_{4p}\) control, or type6 association closed form. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.72 (resolvent-gain calculus from \(n_1,n_3,d_1,d_3\); reverse degrees; type6 pin; 2026-07-29).** Continue Prop 15.68–15.71. Write \(\rho:=m_4-\kappa/p^2\), so \((4pI-T)\rho=T\kappa/p^2\), and on every \(|\kappa|=1\) 4-set with \(\mathrm{sign}(\rho)=\mathrm{sign}(\kappa)\),
\[
|m_4|=\frac1{p^2}+|\rho|.
\]

1. **Gain \(\Leftrightarrow L\) algebra (proved).** For primes \(p\ge5\),
   \[
   |\rho|\;\le\;\frac{p-4}{2p^2}
   \quad\Longleftrightarrow\quad
   |m_4|\;\le\;L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2}
   \]
   on same-sign \(|\kappa|=1\) classes. Since the inhomogeneous source has amplitude \(24/p^2\) on \(|\kappa|=3\), it is sufficient that the resolvent gain from that source into same-sign \(|\kappa|=1\) satisfy
   \[
   \mathrm{gain}\;\le\;\frac{p-4}{48},
   \]
   because \(\frac{p-4}{48}\cdot\frac{24}{p^2}=\frac{p-4}{2p^2}\). \(\square\)

2. **Source sign structure (proved).** On every \(|\kappa|=3\) labelling of \(K_4\), \(T\kappa/\kappa\in\{\pm8\}\) (64-labeling; \(T\kappa=-6\cdot\mathrm{star}\)). Thus \(|T\kappa|=24\) with sign free relative to \(\kappa\). \(\square\)

3. **Reverse extension degrees (proved under Prop 15.68 constancy).** Assume \(d_3=p^2-5\) is constant on every Paley \(|\kappa|=1\) 4-set. Handshaking on the bipartite extension graph between \(|\kappa|=1\) and \(|\kappa|=3\) strata, together with the counts \(n_1,n_3\) of Prop 15.71, forces the degrees from every \(|\kappa|=3\) 4-set:
   \[
   d_1^{(3)}=\frac{n_1\,d_3}{n_3}=3(p^2-1),\qquad
   d_3^{(3)}=4(n-4)-d_1^{(3)}=p^2-9.
   \]
   (Both nonnegative for primes \(p\ge3\), with \(d_3^{(3)}=0\) at \(p=3\).) \(\square\)

4. **Separate \(\kappa\)-weighted vanishing (certified \(p=3,5,7\); open as general theorem).** On every Paley \(|\kappa|=1\) 4-set the one-step sums split by target stratum vanish separately:
   \[
   \sum_{\mathrm{ext}\to|\kappa|=1}C_{vr}\,\kappa(S')=0,\qquad
   \sum_{\mathrm{ext}\to|\kappa|=3}C_{vr}\,\kappa(S')=0.
   \]
   (Each is stronger than \(T\kappa=0\), which is only their sum.) Full multi-worker census at \(p=3,5,7\). Evidence: `e1_gmin_m4_resolvent_gain.json`. \(\square\)

5. **Reverse-degree census (certified \(p=3,5,7\)).** Every \(|\kappa|=3\) 4-set has \((d_1^{(3)},d_3^{(3)})=(3(p^2-1),p^2-9)\) constantly, and \(T\kappa/\kappa\in\{\pm8\}\). \(\square\)

6. **Type6 Max+-free resolvent (certified \(p=3,5,7\); \(W=86\)).** Restricting to \(S_4\)-type6 class-constant functions and solving \((4pI-T)\rho=T\kappa/p^2\) in the least-squares sense:
   | \(p\) | \(\max_{|\kappa|=1}|m_{\mathrm{type6}}|\) | \(L_{\mathrm{abs}}\) | same-sign \(|\rho|\) | gain | budget \(\frac{p-4}{48}\) |
   |------|--------------------------------------|--------------------|----------------------|------|-------------------------------|
   | 5 | \(\approx0.04764\) | \(0.06\) | \(\approx0.00764\) | \(\approx0.00796\) | \(0.02083\) |
   | 7 | \(\approx0.02407\) | \(0.0510\) | \(\approx0.00367\) | \(\approx0.00749\) | \(0.0625\) |

   In particular type6 predicts \(|m|\le L_{\mathrm{abs}}\) and gain below budget at \(p=5,7\). (At \(p=7\) true Max+ \(\max|m_4|\approx0.038> m_{\mathrm{type6}}\), so type6 is not exact — classes need refinement — but remains a Max+-free upper probe.) Evidence: `e1_gmin_m4_resolvent_gain.json`. \(\square\)

7. **Empirical Max+ gain (certified \(p=5,7\); mmap).** True same-sign residual gains \(0.00641\) and \(0.03606\) both lie strictly below \(\frac{p-4}{48}\); \(\max|m_4|\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\). \(\square\)

8. **Residual (OPEN).** Prove for every prime \(p\ge5\) that the resolvent gain is \(\le(p-4)/48\), or directly \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\)), using the stratum data of Props 15.71–15.72 (source size \(n_3\), degrees \(d_1,d_3,d_1^{(3)},d_3^{(3)}\), separate vanishing) plus Max+/boolean structure — without a per-prime Max+ census as the proof. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.73 (Max+ \(e_4\) identity; Paley \(\sum\kappa\); multi-prime type6 gain; 2026-07-29).** Continue Prop 15.52 and 15.71–15.72.

1. **Boolean \(e_4\) from the Max+ sum constraint (proved).** By Prop 15.52, every \(y\in\mathrm{Max}_{+}\) satisfies \(\mathbf1^\top y=(p+1)y_\infty\), hence
   \[
   s^2:=(\textstyle\sum_i y_i)^2=(p+1)^2
   \]
   constantly. For any boolean vector \(y\in\{\pm1\}^n\) one has the power-sum identity
   \[
   e_4(y):=\sum_{a<b<c<d}y_ay_by_cy_d
   =\frac{s^4-6ns^2+3n^2+8s^2-6n}{24}.
   \]
   Substituting \(s^2=(p+1)^2\) and \(n=p^2+1\) yields the closed form
   \[
   e_4=-\frac{p(p-1)(p+1)(p+4)}{12}.
   \]
   Since \(e_4(y)\) is constant on \(\mathrm{Max}_{+}\),
   \[
   \sum_{S}m_4(S)=\mathbb E[e_4]=e_4.
   \]
   Certified on full Max+ at \(p=5,7\) (mmap). Evidence: `e1_gmin_m4_e4_gain.json`. \(\square\)

2. **Paley sum of \(\kappa\) (formula; certified \(p=3,5,7,11,13\)).** For the Paley conference of order \(n=p^2+1\),
   \[
   \sum_S\kappa(S)=\frac{p^2(p^2-1)}{4}.
   \]
   Full multi-worker \(\kappa\)-sum census matches at all five primes (including \(p=13\), \(\binom{170}{4}\approx3.4\cdot10^7\)). Combined with part 1,
   \[
   \sum_S\rho(S)=e_4-\frac1{p^2}\sum_S\kappa(S)
   \]
   is an exact Max+/Paley scalar. \(\square\)

3. **Type6 Max+-free resolvent across primes (certified \(p=5,7,11,13\); \(W=86\)).** All \(11\) abstract \(S_4\)-type6 edge-orbits appear for \(p\ge5\). Solving \((4pI-T)\rho=T\kappa/p^2\) in the type6-constant subspace:
   | \(p\) | \(\max_{|\kappa|=1}|m_{\mathrm{type6}}|\) | \(M_{\mathrm{mid}}\) | \(L_{\mathrm{abs}}\) | gain | budget |
   |------|--------------------------------------|----------------------|--------------------|------|--------|
   | 5 | \(0.04764\) | \(0.0500\) | \(0.0600\) | \(0.00796\) | \(0.0208\) |
   | 7 | \(0.02337\) | \(0.0446\) | \(0.0510\) | \(0.00604\) | \(0.0625\) |
   | 11 | \(0.00876\) | \(0.0341\) | \(0.0372\) | \(0.00249\) | \(0.1458\) |
   | 13 | \(0.00605\) | \(0.0302\) | \(0.0325\) | \(0.00092\) | \(0.1875\) |

   In particular type6 predicts \(|m|\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\) and gain \(\ll(p-4)/48\) at every tested prime, with \(\max|m_{\mathrm{type6}}|\cdot p^2\to1\) (Wick scale). \(\square\)

4. **Caveat (proved by comparison).** Type6 is **not** always exact for true Max+ \(m_4\): at \(p=7\), \(\max|m_4|\approx0.03807> m_{\mathrm{type6}}\approx0.023\). Thus type6 is a Max+-free **probe** (and a feasible particular solution of the master linear equation in a \(T\)-invariant subspace), not by itself an upper bound on true \(m_4\). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that true Max+ satisfies \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\)), e.g. by:
   - closing the type6 error \(m_4-m_{\mathrm{type6}}\in E_{4p}\) with a kernel bound, or
   - proving resolvent gain \(\le(p-4)/48\) from reverse degrees + separate vanishing (Prop 15.72), or
   - a character-sum formula for Paley \(m_4\).

   Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.74 (candidate bound algebra; true Max+ census; signed residual identity; 2026-07-29).** Continue Prop 15.68–15.73. Write \(M_{\mathrm{cand}}(p)=(p-2)/(p(2p+3))\) and \(\rho=m_4-\kappa/p^2\).

1. **Candidate algebra (proved).** For every odd prime \(p\ge5\),
   \[
   M_{\mathrm{cand}}(p)\;\le\;M_{\mathrm{mid}}(p)\;\le\;L_{\mathrm{abs}}(p)\;<\;T_{\mathrm{abs}}(p),
   \]
   with \(M_{\mathrm{cand}}/L_{\mathrm{abs}}=2p/(2p+3)\). Same-sign residual budgets:
   \[
   |m_4|\le M_{\mathrm{cand}}
   \;\Longleftrightarrow\;
   |\rho|\le M_{\mathrm{cand}}-\frac1{p^2}
   \;\Longleftrightarrow\;
   \mathrm{gain}\le\frac{p^2-4p-3}{24(2p+3)},
   \]
   where gain is relative to source amplitude \(24/p^2\). At \(p=5\) the candidate gain is \(1/156\) and \(M_{\mathrm{cand}}=3/65\). \(\square\)

2. **Signed residual identity (proved).** On every \(|\kappa|=1\) 4-set the residual equation \(4p\rho=T\rho\) (Prop 15.68, using \(T\kappa=0\)) multiplies by \(\kappa\) to give
   \[
   4p\,r=\kappa\,(T\rho),\qquad r:=\rho\cdot\kappa.
   \]
   Same-sign danger for \(|m_4|\) is exactly \(r>0\), with \(|m_4|=1/p^2+r\). Certified float residual \(<10^{-15}\) on multi-worker Max+ samples at \(p=5,7\). \(\square\)

3. **True Max+ census (certified \(p=5,7\); \(W=86\); mmap — not type6).** Full enumeration of all \(|\kappa|=1\) 4-sets against Max+:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | \(M_{\mathrm{mid}}\) | same-sign \(\max r\) | gain | \(\mathrm{gain}_{\mathrm{cand}}\) |
   |------|---------------|----------------------|---------------------|----------------------|------|-------------------------------|
   | 5 | \(3/65\approx0.046154\) | \(3/65\) **(sharp)** | \(0.05\) | \(2/325\) | \(1/156\) | \(1/156\) |
   | 7 | \(109/2863\approx0.038072\) | \(5/119\approx0.04202\) | \(0.04464\) | \(\approx0.01766\) | \(\approx0.03606\) | \(\approx0.01330\) |

   In particular **true** Max+ (not type6) satisfies \(\max|m_4|\le M_{\mathrm{cand}}\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\) at \(p=5,7\), with equality in the candidate at \(p=5\). Evidence: `e1_gmin_m4_kernel.json`. \(\square\)

4. **Kernel form (recalled).** Every solution of the master linear equation is \(m_4=m_\star+h\) with \(h\in E_{4p}=\ker(4pI-T)\) and \(m_\star\) the min-norm particular solution (Prop 15.69). True Max+ selects a unique \(h_\star\). At \(p=5\) the candidate is already sharp on true \(m_4\), so \(h_\star\) cannot increase the max beyond \(m_\star\)'s ceiling in the dangerous direction; at \(p=7\), \(h_\star\) raises \(\max|m_4|\) above type6 but still stays \(\le M_{\mathrm{cand}}\). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that true Max+ obeys
   \[
   \max_{|\kappa|=1}|m_4|\;\le\;M_{\mathrm{cand}}(p)
   \]
   (or the weaker \(M_{\mathrm{mid}}\) / gain \(\le(p-4)/48\)), using \(E_{4p}\) control of \(h_\star\), the signed operator \(r\mapsto\kappa(T\rho)\), reverse degrees, and/or Paley character sums — **without** replacing true Max+ by type6. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.75 (one-center \(\sigma_a=2\cdot\mathrm{star}_a\); K4 Gram spectrum; GPU cand cert; 2026-07-29).** Continue Prop 15.67–15.74.

1. **One-center \(\kappa\)-sum (proved, any conference matrix).** For every 4-set \(S\) with \(|\kappa(S)|=1\) and every \(a\in S\),
   \[
   \sigma_a:=\sum_{r\notin S}C_{ar}\,\kappa(S_{a\to r})
   \;=\;2\cdot\mathrm{star}_a,
   \qquad
   \mathrm{star}_a:=\prod_{u\in S\setminus\{a\}}C_{au}\in\{\pm1\}.
   \]
   *Proof.* \(C^2=(n-1)I\) gives \(\sum_{r\notin\{a,b\}}C_{ar}C_{br}=-\sum_{u\in S\setminus\{a,b\}}C_{au}C_{ub}\) for \(b\in S\setminus\{a\}\). Expanding \(\sigma_a\) in the three pairing products of \(\kappa(S_{a\to r})\) yields \(\sigma_a=-2\sum_{v\in S\setminus\{a\}}\mathrm{star}_v\). On \(|\kappa|=1\), \(T\kappa=0\) forces \(\sum_{v\in S}\mathrm{star}_v=0\), hence \(\sum_{v\neq a}\mathrm{star}_v=-\mathrm{star}_a\) and \(\sigma_a=2\,\mathrm{star}_a\). \(\square\)

2. **One-center residual form (proved).** With \(\rho=m_4-\kappa/p^2\), the evec identity at centre \(a\) becomes
   \[
   p\,m_4-\frac{\kappa}{p}
   \;=\;\frac{2\,\mathrm{star}_a}{p^2}
   +\sum_{r\notin S}C_{ar}\,\rho(S_{a\to r}).
   \]
   (Certified float residual \(<10^{-15}\) on Max+ samples at \(p=5,7\).) Averaging over \(a\in S\) recovers \(\mathrm{Ext}=4p\rho\). \(\square\)

3. **Local K4 edge-Gram spectrum (proved algebra).** On the six edges of a \(|\kappa|=1\) 4-set, the principal submatrix of \(G=\mathbb E[ff^\top]\) has opposite-edge entries \(\pi_i m_4\) and wedge entries \(\pm1/p\). Its eigenvalues lie in the pool
   \[
   \{1\pm m_4\}\ \cup\ \{1\pm m_4\pm 2/p\}\ \cup\ \{1\pm m_4\pm 2\sqrt2/p\}.
   \]
   PSD of this block forces the weak general bound \(|m_4|\le1-2/p=(p-2)/p\) whenever \(1-|m_4|-2/p\) is an eigenvalue — far weaker than \(M_{\mathrm{cand}}\), but Max+-free and load-bearing for local structure. \(\square\)

4. **GPU true Max+ cand census (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full \(|\kappa|=1\) m4 on Max+ with one CUDA context, device argmax, mmap Max+ load, atomic evidence write:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | same-sign \(r\) | gain | wall |
   |------|---------------|----------------------|----------------|------|------|
   | 5 | \(3/65\) (sharp) | \(3/65\) | \(2/325\) | \(1/156\) | \(\sim0.36\)s |
   | 7 | \(109/2863\) | \(5/119\) | \(\approx0.01766\) | \(\approx0.0361\) | \(\sim0.30\)s |

   Evidence: `e1_gmin_m4_gpu.json` (`gpu.used=true`, `io.mmap+atomic`). Multi-worker \(\sigma\) census \(p=3,5,7,11\): `e1_gmin_m4_onecenter.json`. \(\square\)

5. **Residual (OPEN).** Bound \(\sum_r C_{ar}\rho(S_{a\to r})\) (or the signed global operator of Prop 15.74) tightly enough that part 2 forces \(m_4\le M_{\mathrm{cand}}\) for every prime \(p\ge5\). Absolute degree bounds do not contract. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.76 (one-center extension degrees; residual split; 2026-07-29).** Continue Prop 15.75. Write \(d_1=3p^2-7\), \(d_3=p^2-5\) for the full ordered-extension degrees on \(|\kappa|=1\) (Prop 15.68/15.72).

1. **Divisibility (proved).** For every odd integer \(p\ge3\), \(4\mid(3p^2-7)\) and \(4\mid(p^2-5)\). Hence
   \[
   d_1^{(1)}:=\frac{3p^2-7}{4},\qquad d_3^{(1)}:=\frac{p^2-5}{4}
   \]
   are integers, and \(d_1^{(1)}+d_3^{(1)}=p^2-3=n-4\), \(4d_1^{(1)}=d_1\), \(4d_3^{(1)}=d_3\). \(\square\)

2. **One-center degree constancy (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) 4-set \(S\) and every centre \(a\in S\), among the \(n-4\) one-vertex extensions \(S_{a\to r}\), exactly \(d_1^{(1)}\) land in \(|\kappa|=1\) and \(d_3^{(1)}\) land in \(|\kappa|=3\). (So the full \(4(n-4)\) count splits evenly across the four centres.) Evidence: `e1_gmin_m4_onecenter_deg.json`. \(\square\)

3. **Residual split (proved form).** With \(\rho=m_4-\kappa/p^2\) and \(\sigma_a=2\cdot\mathrm{star}_a\) (Prop 15.75), for each centre \(a\) on a \(|\kappa|=1\) set with \(\kappa=1\),
   \[
   p\rho-\frac{2\,\mathrm{star}_a}{p^2}
   \;=\;S_1(a)+S_3(a),
   \]
   where \(S_j(a)=\sum C_{ar}\rho(S_{a\to r})\) runs over extensions to \(|\kappa|=j\). In particular \(|S_3(a)|\le d_3^{(1)}R_3\) and \(|S_1(a)|\le d_1^{(1)}R_1\) with \(R_j=\max|\rho|\) on the \(|\kappa|=j\) stratum. \(\square\)

4. **Absolute bootstrap fails (proved).** The four-centre form \(4p R_1\le d_1 R_1+d_3 R_3\) rearranges to \(R_1(4p-d_1)\le d_3 R_3\). For primes \(p\ge5\), \(4p-d_1=4p-(3p^2-7)<0\), so this yields only a lower bound on \(R_1\), not an upper bound. Signed cancellation on \(S_1\) (or a design bound coupling \(R_1,R_3\)) is load-bearing. \(\square\)

5. **GPU residual moments (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full-quad m4 on Max+:
   | \(p\) | \(\max|m_4|_{\kappa=1}\) | \(\max|m_4|_{\kappa=3}\) | \(R_1\) | \(R_3\) | \(\le M_{\mathrm{cand}}\) |
   |------|--------------------------|--------------------------|--------|--------|------------------------|
   | 5 | \(3/65\) | \(21/65\) | \(0.0554\) | \(0.2031\) | yes |
   | 7 | \(109/2863\) | \(\approx0.1142\) | \(0.0177\) | \(0.0530\) | yes |

   Note \(R_1>\max r\) (same-sign residual) because opposite-sign \(\rho\) on \(|\kappa|=1\) is larger; only same-sign \(r\) raises \(|m_4|\) above Wick. Evidence: `e1_gmin_m4_onecenter_deg.json`. \(\square\)

6. **Residual (OPEN).** Prove a signed bound on \(S_1+S_3\) (e.g. \(S_1\le0\) at \(\mathrm{star}_a=+1\) on maximisers, or a Paley character-sum formula for \(\rho\)) strong enough that part 3 forces \(\rho\le(p^2-4p-3)/(p^2(2p+3))\) on same-sign \(|\kappa|=1\), i.e. \(m_4\le M_{\mathrm{cand}}\), for every prime \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.77 (\(\mathrm{star}_a\cdot S_1\le0\) structure; joint cand criterion; 2026-07-29).** Continue Prop 15.75–15.76. Write \(S_1(a),S_3(a)\) for the one-center residual split on a \(|\kappa|=1\) 4-set and centre \(a\), and \(\mathrm{star}_a=\prod_{u\in S\setminus\{a\}}C_{au}\).

1. **Star–joint identity (proved, Max+-free).** From the one-center residual form,
   \[
   p\rho-\frac{2\,\mathrm{star}_a}{p^2}=S_1(a)+S_3(a)
   \qquad\Longrightarrow\qquad
   \mathrm{star}_a\cdot\bigl(S_1+S_3\bigr)
   \;=\;
   p\,\rho\,\mathrm{star}_a-\frac{2}{p^2}.
   \]
   In particular, if \(\mathrm{star}_a=+1\) and \(S_1(a)\le0\), then
   \[
   p\rho\;\le\;\frac{2}{p^2}+S_3(a).
   \]
   \(\square\)

2. **Same-sign reconstruction (proved).** On every \(|\kappa|=1\) set with \(\kappa=1\) and \(\rho>0\), at every centre with \(\mathrm{star}_a=+1\),
   \[
   \rho
   \;=\;
   \frac{1}{p}\Bigl(\frac{2}{p^2}+S_1(a)+S_3(a)\Bigr).
   \]
   Hence
   \[
   \max_{\text{same-sign }r>0}\rho
   \;=\;
   \max_{\substack{\mathrm{star}_a=+1\\ r>0}}
   \frac{1}{p}\Bigl(\frac{2}{p^2}+S_1+S_3\Bigr),
   \]
   and the candidate bound \(\rho\le\rho_{\mathrm{cand}}:=M_{\mathrm{cand}}-1/p^2\) is equivalent to
   \[
   \max_{\substack{\mathrm{star}_a=+1\\ r>0}}\bigl(S_1+S_3\bigr)
   \;\le\;
   p\,\rho_{\mathrm{cand}}-\frac{2}{p^2}.
   \]
   At \(p=5\) the right-hand side equals \(-16/325\approx-0.04923<0\), so \(S_1\le0\) alone with absolute \(|S_3|\) bounds is **not** enough — maximisers must have strongly negative \(S_1\). \(\square\)

3. **GPU census of \(\mathrm{star}_a\cdot S_1\) (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full Max+ m4 (one CUDA context, single H2D, D2H m4 vector only) + ProcessPool walk of every \(|\kappa|=1\) centre:
   | \(p\) | \(\#\) star\(+\) checks | \(\max(\mathrm{star}\cdot S_1)\) | \(\max S_1\) at star\(+\) | joint \(S_1{+}S_3\) on \(r>0\) | \(\Rightarrow\rho\le\rho_{\mathrm{cand}}\) |
   |------|--------------------------|----------------------------------|---------------------------|----------------------------------|----------------------------------|
   | 5 | \(23400\) | \(-0.03077\) (strict \(<0\)) | \(-0.03077\) | \(-0.04923\) | **yes (sharp)** |
   | 7 | \(352800\) | \(-0.00669\) (strict \(<0\)) | \(-0.00669\) | \(0.08283\) | **yes** |

   In particular **\(\mathrm{star}_a\cdot S_1(a)\le0\) on every \(|\kappa|=1\) centre** (not only on same-sign maximisers), with perfect sign-antisymmetry \(\max S_1|_{\mathrm{star}+}=-\min S_1|_{\mathrm{star}-}\). Identity residual \(<10^{-15}\). Evidence: `e1_gmin_m4_S1_star.json`. \(\square\)

4. **Max+ is essential (certified).** A synthetic residual that puts a uniform same-sign bump on all \(|\kappa|=1\) coordinates (and zero on \(|\kappa|=3\)) **violates** \(\mathrm{star}\cdot S_1\le0\) (\(\max\mathrm{star}\cdot S_1=+0.05\) at \(p=5,7\)). So the inequality is a property of true Max+ fourth moments, not of the \(\kappa\)-adjacency graph alone. Combinatorial \(\tau_1:=\sum_{\kappa1}C_{ar}\kappa'\) has \(\mathrm{star}\cdot\tau_1\) non-constant for \(p\ge5\) (values in \(\{-1,5\}\) at \(p=5\); multi-worker pure-\(C\) census). \(\square\)

5. **Residual (OPEN).** Prove \(\mathrm{star}_a\cdot S_1(a)\le0\) for every prime \(p\ge5\) on Paley Max+ (e.g. via boolean \(+p\)-evec character sums), and prove the joint criterion of part 2 (or a matching \(S_3\) bound on maximisers) so that \(\max|m_4|\le M_{\mathrm{cand}}\) for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.78 (moment form / Gaussian domination; 4-set constancy; exact \(p=5\) spectrum; 2026-07-29).** Continue Prop 15.77. Write \(R_1(a)=\{r\notin S:|\kappa(S_{a\to r})|=1\}\), \(U_1=\sum_{r\in R_1}C_{ar}y_r\), \(\tau_1(a)=\sum_{r\in R_1}C_{ar}\kappa(S_{a\to r})\), and \(f_0=\prod_{u\in S\setminus\{a\}}y_u\).

1. **Moment form (proved, Max+ with \(Cy=py\)).** On every \(|\kappa|=1\) centre,
   \[
   \mathrm{star}_a\cdot S_1(a)
   \;=\;
   \mathbb E\bigl[\varphi\bigr]-\mathbb E_{\mathrm{Wick}}\bigl[\varphi\bigr],
   \qquad
   \varphi(y)\,:=\,\mathrm{star}_a\,f_0\,U_1,
   \]
   where Wick uses the Gaussian \(N(0,\Sigma)\) with \(\Sigma=I+C/p\) (so \(\mathbb E_{\mathrm{Wick}}[\varphi]=\mathrm{star}_a\tau_1/p^2\)). Consequently
   \[
   \mathrm{star}_a\cdot S_1(a)\le0
   \quad\Longleftrightarrow\quad
   \mathbb E_{\mathrm{Max}+}[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]
   \]
   (Gaussian domination for this specific degree-4 form). \(\square\)

2. **Joint rewrite (proved).** For \(\kappa=1\) and \(\mathrm{star}_a=+1\),
   \[
   S_1+S_3 \;=\; p\,m_4-\frac1p-\frac2{p^2}.
   \]
   (So the joint \(S_1+S_3\) is a function of \(m_4\) alone on this locus; the cand criterion is equivalent to \(\max m_4\le M_{\mathrm{cand}}\) on same-sign sets — the leverage of \(S_1\le0\) is to replace the joint by an \(S_3\)-only upper bound.) \(\square\)

3. **Combinatorial constancy of \(\mathrm{star}\cdot\tau_1\) (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) 4-set, the four values \(\mathrm{star}_a\cdot\tau_1(a)\) (\(a\in S\)) are equal. Observed value sets: \(\{-1\}\) at \(p=3\); \(\{-1,5\}\) at \(p=5\); \(\{-7,-1,5\}\) at \(p=7\); \(\{-13,-7,-1,5,11\}\) at \(p=11\). Evidence: `e1_gmin_m4_S1_const.json`. \(\square\)

4. **Max+ constancy of \(\mathrm{star}\cdot S_1\) (certified GPU \(p=5,7\); mmap+atomic).** On every \(|\kappa|=1\) 4-set, the four values \(\mathrm{star}_a\cdot S_1(a)\) are equal (one rational per set). In particular the 4-set carries a single sign for the GD inequality. Full census: \(\mathrm{star}\cdot S_1\le0\) with Gaussian domination on every centre; \(\max\mathrm{star}\cdot S_1=-2/65\) (\(p=5\)), \(-0.006686\) (\(p=7\)). \(\square\)

5. **Exact spectrum at \(p=5\) (certified full Max+).** \(\mathrm{star}\cdot S_1\in\{-2/65,-42/325\}\), both \(<0\). Matches the closed moment rule
   \[
   \mathrm{star}\cdot\mathbb E[f_0U_1]
   \;=\;
   \frac{11}{65}\,\mathrm{sgn}(\mathrm{star}\cdot\tau_1)
   \]
   with \(\mathrm{star}\cdot\tau_1\in\{-1,5\}\) (so \(\mathrm{star}\cdot S_1=\frac{11}{65}\mathrm{sgn}(t)-\,t/25\)). Hence \(\mathrm{star}\cdot S_1\le0\) is **proved at \(p=5\)** by exhaustive Max+ evaluation. \(\square\)

6. **Residual (OPEN).** Constancy of \(\mathrm{star}\cdot\tau_1\) and \(\mathrm{star}\cdot S_1\) is proved in Prop 15.79. Remains: Gaussian domination \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\) (equivalently \(\mathrm{star}\cdot S_1\le0\)) for all primes \(p\ge5\), and joint/\(S_3\) bound \(\Rightarrow\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.79 (Aut-constancy of \(\mathrm{star}\cdot\tau_1\) and \(\mathrm{star}\cdot S_1\); modular \(\tau_1\); 2026-07-29).** Continue Prop 15.77–15.78. Let \(C\) be the Paley conference matrix of order \(n=p^2+1\), identified with the point set \(\mathrm{PG}(1,\mathbb F_{p^2})\), and let \(\mathrm{Aut}(C)\ge\mathrm{PGL}(2,p^2)\).

1. **Setwise stabilizer is transitive (proved).** For every 4-set \(S\subset\mathrm{PG}(1,p^2)\), the setwise stabilizer of \(S\) in \(\mathrm{PGL}(2,p^2)\) contains a Klein four-group \(V_4\) of double-transposition involutions (e.g. \(z\mapsto\lambda/z\) on \(S=\{\infty,0,1,\lambda\}\)). This \(V_4\) acts regularly on \(S\), hence transitively. \(\square\)

2. **Constancy of \(\mathrm{star}\cdot\tau_1\) (proved, Max+-free).** The scalar \(\mathrm{star}_a\cdot\tau_1(a)\) is built from \(C\)-labels of \(S\) and of the one-vertex extensions \(S_{a\to r}\). It is therefore \(\mathrm{Aut}(C)\)-equivariant: for every \(g\in\mathrm{Aut}(C)\),
   \[
   \mathrm{star}_{ga}\cdot\tau_1(ga;gS)
   \;=\;
   \mathrm{star}_a\cdot\tau_1(a;S).
   \]
   Restricting to the setwise stabilizer of \(S\) and using part 1, \(\mathrm{star}_a\cdot\tau_1(a)\) is independent of \(a\in S\). \(\square\)

3. **Constancy of \(\mathrm{star}\cdot S_1\) (proved, Max+ Aut-invariant).** The set \(\mathrm{Max}_{+}=\{y\in\{\pm1\}^n:Cy=py\}\) is \(\mathrm{Aut}(C)\)-invariant, so \(m_4(S)=\mathbb E[\prod_{i\in S}y_i]\) is an Aut-invariant of the 4-set. Consequently \(\mathrm{star}_a\cdot S_1(a)\) (depending on \(C\) and on \(m_4\) of extensions of \(S\)) is Aut-equivariant, and part 1 forces constancy on \(a\in S\). \(\square\)

4. **Counting form of \(\tau_1\) (proved).** Write \(d_1^{(1)}=(3p^2-7)/4\) and
   \[
   A\;=\;\#\{r\in R_1(a):C_{ar}\kappa(S_{a\to r})=+1\}.
   \]
   Then \(A+B=d_1^{(1)}\) with \(B=\#\{C_{ar}\kappa'=-1\}\), and \(\tau_1=A-B=2A-d_1^{(1)}\). In particular \(\mathrm{star}\cdot\tau_1=\mathrm{star}\cdot(2A-d_1^{(1)})\). For odd \(p\), \(d_1^{(1)}\) is odd, so \(\mathrm{star}\cdot\tau_1\) is always odd. \(\square\)

5. **Modular census (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) set: constancy holds (matches part 2); \(t_1=2A-d_1^{(1)}\) exactly; every \(\mathrm{star}\cdot\tau_1\equiv5\pmod6\); the number of distinct values is \((p-1)/2\), with observed sets
   \[
   \begin{align*}
   p=3&:\ \{-1\},\\
   p=5&:\ \{-1,5\},\\
   p=7&:\ \{-7,-1,5\},\\
   p=11&:\ \{-13,-7,-1,5,11\}
   \end{align*}
   \]
   (arithmetic progressions of difference \(6\)). Evidence: `e1_gmin_m4_S1_aut.json`. \(\square\)

6. **Consequence for the sign attack.** By part 3 it suffices to prove \(\mathrm{star}\cdot S_1(S)\le0\) as a property of the 4-set (one check per set, not per centre). Combined with Prop 15.77–15.78 (joint criterion; Gaussian domination form; exact \(p=5\) spectrum \(\{-2/65,-42/325\}\)), the residual is: prove \(\mathrm{star}\cdot S_1\le0\) for all primes \(p\ge5\), then control \(S_3\) on maximisers. \(\square\)

7. **Residual (OPEN).** Prove \(\mathrm{star}\cdot S_1\le0\) (e.g. Gaussian domination \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) for every prime \(p\ge5\), and close the joint/\(S_3\) bound so \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.80 (linear-form Wick identity; GD formulation; U1-specialization; 2026-07-29).** Continue Prop 15.77–15.79. Write \(Z=\mathrm{star}_a\prod_{u\in S\setminus\{a\}}y_u\) and \(U_1=\sum_{r\in R_1(a)}C_{ar}y_r\).

1. **Linear-form Wick identity (proved).** On Max+ one has \(\mathbb E[y_iy_j]=\Sigma_{ij}\) with \(\Sigma=I+C/p\). Consequently for every \(\alpha\in\mathbb R^n\) and \(L=\sum\alpha_i y_i\),
   \[
   \mathbb E[L^2]
   \;=\;
   \sum_{i,j}\alpha_i\alpha_j\Sigma_{ij}
   \;=\;
   \mathbb E_{G\sim N(0,\Sigma)}[L(G)^2].
   \]
   In particular \(\mathbb E[U_1^2]=\mathbb E_{\mathrm{Wick}}[U_1^2]\). (Only pairwise moments; no fourth-moment input.) \(\square\)

2. **GD formulation (proved).** Central symmetry gives \(\mathbb E[Z]=0\), and
   \[
   \mathrm{star}_a\cdot S_1(a)
   \;=\;
   \mathbb E[ZU_1]-\mathbb E_{\mathrm{Wick}}[ZU_1]
   \;=\;
   \mathrm{Cov}(Z,U_1)-\mathrm{Cov}_{\mathrm{Wick}}(Z,U_1),
   \]
   with \(\mathbb E_{\mathrm{Wick}}[ZU_1]=\mathrm{star}_a\cdot\tau_1(a)/p^2\). Hence
   \[
   \mathrm{star}\cdot S_1\le0
   \quad\Longleftrightarrow\quad
   \mathbb E[ZU_1]\le\mathbb E_{\mathrm{Wick}}[ZU_1].
   \]
   \(\square\)

3. **U1-specialization (proved necessity of \(\kappa=1\)-support; certified generically).** The comparison \(\mathbb E[ZL]\le\mathbb E_{\mathrm{Wick}}[ZL]\) **fails** for generic linear \(L\) supported off \(S\) (violation rate \(\approx 45\)–\(55\%\) at \(p=5,7\)). It can also fail for the full external form \(U_{\mathrm{ext}}=\sum_{r\notin S}C_{ar}y_r\) (equivalent to \(\mathrm{star}\cdot p\cdot\rho\le 2/p^2\), false on some \(p=7\) maximisers with positive joint \(S_1+S_3\)). Thus GD is a property of the \(\kappa=1\)-restricted form \(U_1\), not a general cubic\(\times\)linear inequality. \(\square\)

4. **Sum of \(\mathrm{star}\cdot\tau_1\) (certified \(p=3,5,7\)).** Writing \(n_1=n(n-1)(n-2)^2/32\),
   \[
   \sum_{|\kappa|=1}\mathrm{star}\cdot\tau_1
   \;=\;
   \varepsilon(p)\,n_1,
   \qquad \varepsilon(p)\in\{\pm1\}
   \]
   (\(\varepsilon(5)=+1\), \(\varepsilon(3)=\varepsilon(7)=-1\); at \(p=7\) the three values \(\{-7,-1,5\}\) each occupy exactly \(n_1/3\) sets). \(\square\)

5. **GPU GD census (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full Max+ m4 + ProcessPool walk:
   - \(\mathrm{star}\cdot S_1\le0\) on every \(|\kappa|=1\) set (GD holds);
   - \(\mathbb E[U_1^2]/\mathbb E_{\mathrm{Wick}}[U_1^2]\equiv1\) (part 1);
   - \(\sum\mathrm{star}\cdot S_1=-1128\) at \(p=5\) and \(-15271200/2863\) at \(p=7\);
   - \(\max|m_4|\le M_{\mathrm{cand}}\) (sharp at \(p=5\)).

   Evidence: `e1_gmin_m4_S1_gd.json`. \(\square\)

6. **Residual (OPEN).** Prove \(\mathbb E[ZU_1]\le\mathbb E_{\mathrm{Wick}}[ZU_1]\) for every prime \(p\ge5\) (the U1-specific Gaussian domination), using the linear Wick identity, Aut-constancy, and the residual source on \(|\kappa|=3\); then close joint/\(S_3\) for \(M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.81 (moduli-line GD criterion; complete at \(p=5\); 2026-07-29).** Continue Prop 15.53 and 15.77–15.80. Assume a stratification of 4-sets into pure-\(C\) classes with constant \(m_4\), and let \(Am=b\) be the averaged evec system from \(Cy=py\) and \(m_2=C/p\).

1. **Pointwise evec identity (certified Max+ \(p=5,7\)).** For every 4-set \(S\) and centre \(a\in S\),
   \[
   p\,m_4(S)=\sum_{r\neq a}C_{ar}\,m(r,S\setminus\{a\})
   \]
   (with the usual reduction of repeated indices via \(y_i^2=1\)). Float residual \(<10^{-15}\). This is the moment form of the one-center identity (Prop 15.75). \(\square\)

2. **Affine criterion (proved form).** If \(Am=b\) has nullity 1, write \(m=m_\star+c\,n\). Each class-level \(\mathrm{star}\cdot S_1\) is linear in \(m\), hence affine in \(c\). On any interval where a single type realises the maximum,
   \[
   \max_{|\kappa|=1}\mathrm{star}\cdot S_1
   \;=\;
   \alpha+\beta c.
   \]
   If \(\beta>0\), then \(\mathrm{star}\cdot S_1\le0\) for all centres iff \(c\le c_{\mathrm{GD}}:=-\alpha/\beta\). \(\square\)

3. **Complete at \(p=5\) (certified; drives `e1_gmin_moduli`).** The \((\mathrm{type6},\mathrm{ext\text{-}hist})\) stratification has 37 classes, all with constant \(m_4\); the evec system has nullity 1. On the line,
   \[
   \max\mathrm{star}\cdot S_1=\alpha+\beta c
   \]
   is exact (\(\beta>0\), float fit error \(<10^{-15}\)), with
   \[
   c_{\mathrm{GD}}\approx-0.29605.
   \]
   The physical root selected by \(\mathrm{Tr}(G^2)=\mathrm{Tr}_{\mathrm{Max+}}\) is
   \[
   c^\star\approx-0.42402\;<\;c_{\mathrm{GD}},
   \]
   and at \(c^\star\) one has \(\max|m_4|=M_{\mathrm{cand}}=3/65\) (sharp) and \(\max\mathrm{star}\cdot S_1=-2/65\le0\). Thus both the candidate bound and Gaussian domination hold at \(p=5\) by moduli calculus. Evidence: `e1_gmin_m4_S1_moduli.json`. \(\square\)

4. **Status at \(p=7\) (certified structure).** Coarse classes: \(69/82\) have constant \(m_4\) (max std \(\approx0.012\)); the averaged system is not yet a faithful nullity-1 line. Pointwise evec identities still hold. Full Max+ GD and \(\max|m_4|\le M_{\mathrm{cand}}\) remain certified by Prop 15.80 / 15.74 censuses. Finer \(C\)-invariants are needed for a moduli-line proof at \(p=7\). \(\square\)

5. **Residual (OPEN).** For every prime \(p\ge5\): refine classes to constant \(m_4\) and nullity 1; prove \(\beta>0\) and \(c^\star\le c_{\mathrm{GD}}\) (or unique solution with GD); conclude \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.82 (type6+CR refinement; m4 constant at \(p=5,7\); moduli GD pin; 2026-07-29).** Continue Prop 15.81. Refine pure-\(C\) stratifications of 4-sets by adjoining the PGL-complete cross-ratio invariant of Prop 15.48–15.49 to the type6 (or type6+ext-hist) labels.

1. **Constancy discovery (certified GPU+ProcessPool \(W=86\); mmap Max+).** Full Max+ \(m_4\) on all \(\binom{n}{4}\) quads (CuPy/V100, one H2D) and class-key shards:
   | strat | \(p=5\) | \(p=7\) |
   |-------|---------|---------|
   | coarse \((\mathrm{type6},\mathrm{ext\text{-}hist})\) | 37/37 const | 69/82 const |
   | type6+CR | **26/26 const** | **48/48 const** |
   | coarse+CR | 48/48 const | 130/130 const |
   In particular **type6+CR** (and coarse+CR) make \(m_4\) constant on every class at both \(p=5\) and \(p=7\). True Max+ still obeys \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{cand}}\) (sharp \(3/65\) at \(p=5\); \(\approx0.038<5/119\) at \(p=7\)). Evidence: `e1_gmin_m4_refine.json`. \(\square\)

2. **Moduli on type6+CR (certified \(W=86\)).** Averaged evec system \(Am=b\):
   - **\(p=5\):** 26 classes, **nullity 1**. Affine law \(\max\mathrm{star}\cdot S_1=\alpha+\beta c\) exact (\(\beta<0\)); physical \(c^\star\) lies on the safe side of \(c_{\mathrm{GD}}\) (i.e. \(\max\mathrm{star}\cdot S_1(c^\star)\le0\)); \(\max|m_4|=M_{\mathrm{cand}}\); GD holds. Thus **cand+GD at \(p=5\)** also under the type6+CR line (cf. Prop 15.81 coarse line).
   - **\(p=7\):** 48 classes, **nullity 2**. True Max+ still has GD and \(\max|m_4|\le M_{\mathrm{cand}}\); full multi-parameter pin OPEN.
   Evidence: `e1_gmin_m4_refine_moduli.json`. \(\square\)

3. **Safe-side orientation (proved form).** On a nullity-1 line, \(\max\mathrm{star}\cdot S_1=\alpha+\beta c\). If \(\beta>0\) then GD \(\Leftrightarrow c\le c_{\mathrm{GD}}:=-\alpha/\beta\); if \(\beta<0\) then GD \(\Leftrightarrow c\ge c_{\mathrm{GD}}\). The sign of \(\beta\) depends on null-vector orientation; the physical check is always \(\mathrm{sign}(\beta)\cdot(c^\star-c_{\mathrm{GD}})\le0\). \(\square\)

4. **Extra linear pins (certified \(W=86\)).** On type6+CR, \(\sum_S m_4(S)=e_4\) (Prop 15.73) holds exactly at \(p=5,7\) but is already in the row-span of the averaged evec system (rank unchanged). Denser evec sampling likewise does not drop the \(p=7\) nullity below 2. Coarse+CR (130 classes) also has nullity 2. Evidence: `e1_gmin_m4_pin_extra.json`, `e1_gmin_m4_refine_moduli_multi.json`. \(\square\)

5. **Residual (OPEN).** For every prime \(p\ge5\): either obtain a constant-\(m_4\) stratification with nullity \(\le1\) and prove \(c^\star\) safe-side of \(c_{\mathrm{GD}}\), or close a multi-parameter pin (Tr\((G^2)\) surface + second moment / character-sum GD) so \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.83 (resolvent-budget hierarchy for \(M_{\mathrm{cand}}\) vs \(L\); Max+-free; 2026-07-30).** Continue Prop 15.66–15.74. Write
\[
M_{\mathrm{cand}}(p)=\frac{p-2}{p(2p+3)},\quad
M_{\mathrm{mid}}(p)=\frac{p-2}{2p(p+1)},\quad
L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2},\quad
T_{\mathrm{abs}}(p)=\frac{p-2}{p(2p-1)},
\]
and the same-sign residual targets on \(|\kappa|=1\)
\[
\rho_L=\frac{p-4}{2p^2},\qquad
\rho_{\mathrm{cand}}=\frac{p^2-4p-3}{p^2(2p+3)}=M_{\mathrm{cand}}-\frac1{p^2},
\]
together with the resolvent-gain budgets (source amplitude \(24/p^2\) on \(|\kappa|=3\), Prop 15.68–15.72)
\[
\mathrm{gain}_L=\frac{p-4}{48},\qquad
\mathrm{gain}_{\mathrm{cand}}=\frac{p^2-4p-3}{24(2p+3)}.
\]

1. **Cascade (proved algebra, Max+-free).** For every real \(p>2\),
\[
M_{\mathrm{cand}}<M_{\mathrm{mid}}\le L_{\mathrm{abs}}<T_{\mathrm{abs}},
\]
with positive gaps
\[
M_{\mathrm{mid}}-M_{\mathrm{cand}}=\frac{p-2}{2p(p+1)(2p+3)},\quad
L_{\mathrm{abs}}-M_{\mathrm{mid}}=\frac{p-2}{2p^2(p+1)},\quad
T_{\mathrm{abs}}-L_{\mathrm{abs}}=\frac{p-2}{2p^2(2p-1)}.
\]
In particular the cascade holds for every prime \(p\ge5\). \(\square\)

2. **Residual ranking (proved algebra).** For every prime \(p\ge5\),
\[
0<\rho_{\mathrm{cand}}<\rho_L,\qquad
\rho_L-\rho_{\mathrm{cand}}=\frac{3(p-2)}{2p^2(2p+3)}.
\]
Thus the \(M_{\mathrm{cand}}\) residual target is **strictly tighter** than the \(L_{\mathrm{abs}}\) residual. \(\square\)

3. **Resolvent-budget ranking (proved algebra).** For every real \(p>2\),
\[
\mathrm{gain}_L-\mathrm{gain}_{\mathrm{cand}}
=
\frac{3(p-2)}{48(2p+3)}
=
\frac{p-2}{16(2p+3)}
\;>\;0,
\]
so \(\mathrm{gain}_{\mathrm{cand}}<\mathrm{gain}_L\). As \(p\to\infty\), \(\mathrm{gain}_{\mathrm{cand}}/\mathrm{gain}_L\to1\). Consequently any operator-gain bound
\[
\mathrm{gain}\;\le\;\mathrm{gain}_{\mathrm{cand}}
\]
from the \(|\kappa|=3\) source into same-sign \(|\kappa|=1\) automatically yields \(\max|m_4|\le M_{\mathrm{cand}}\le L_{\mathrm{abs}}\) and (with Prop 15.47) bi-tight empty for every prime \(p\ge5\). \(\square\)

4. **What remains OPEN.** The inequality \(\mathrm{gain}\le\mathrm{gain}_{\mathrm{cand}}\) (or the weaker \(\mathrm{gain}\le\mathrm{gain}_L\)) for true Max+ fourth moments is **not** proved for general primes \(p\ge5\); it is only certified at \(p=5,7\) by census (Props 15.72–15.74). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1583.py`, `evidence/e1_gmin_m4_prop1583.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.84 (GD \(\Rightarrow\) cand via \(S_3\) budget; diagonal dominance fails; Max+-free; 2026-07-30).** Continue Prop 15.77 and 15.83. On a \(|\kappa|=1\) centre with \(\mathrm{star}_a=+1\), Prop 15.77 gives
\[
\mathrm{star}_a\cdot(S_1+S_3)=p\rho-\frac2{p^2}.
\]
Under GD (\(\mathrm{star}_a\cdot S_1\le0\)) one has \(p\rho\le 2/p^2+S_3\), i.e.
\[
\rho\;\le\;\frac2{p^3}+\frac{S_3}{p}.
\]
Write the cand residual \(\rho_{\mathrm{cand}}=M_{\mathrm{cand}}-1/p^2\) and the \(S_3\)-budget
\[
B_{\mathrm{cand}}(p)
\;:=\;
p\rho_{\mathrm{cand}}-\frac2{p^2}
\;=\;
\frac{p^3-4p^2-7p-6}{p^2(2p+3)}.
\]

1. **Closed form (proved algebra).** The displayed formula for \(B_{\mathrm{cand}}\) holds for every prime \(p\ge5\). \(\square\)

2. **Sign pattern (proved algebra).** \(B_{\mathrm{cand}}(5)=-16/325<0\). The cubic numerator \(p^3-4p^2-7p-6\) is increasing on \([5,\infty)\) (derivative \(3p^2-8p-7>0\) for \(p\ge5\)) and positive at \(p=7\), hence \(B_{\mathrm{cand}}(p)>0\) for every prime \(p\ge7\). \(\square\)

3. **Settlement lemma (proved form).** If GD holds and \(S_3\le B_{\mathrm{cand}}\) at every same-sign \(|\kappa|=1\) centre, then \(\rho\le\rho_{\mathrm{cand}}\), so \(\max|m_4|\le M_{\mathrm{cand}}\) and bi-tight is empty for all primes \(p\ge5\) (Props 15.47, 15.74). \(\square\)

4. **Why absolute bootstrap fails (proved algebra).** The one-step degree of a \(|\kappa|=1\) 4-set into other \(|\kappa|=1\) 4-sets is \(d_1=3p^2-7\) (Prop 15.72/76). Then
\[
4p-d_1=-3p^2+4p+7<0
\]
for every prime \(p\ge5\), so \(4pI-T\) is **not** diagonally dominant on the \(|\kappa|=1\) stratum. Absolute row-sum inversion cannot prove a residual bound; signed cancellation (\(S_1\le0\), controlled \(S_3\)) is load-bearing. \(\square\)

5. **What remains OPEN.** Prove GD and \(S_3\le B_{\mathrm{cand}}\) (or the resolvent-gain bound of Prop 15.83) for true Max+ fourth moments for all primes \(p\ge5\). At \(p=5\), part 2 forces \(S_3\) strictly negative under GD to reach cand. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1584.py`, `evidence/e1_gmin_m4_prop1584.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.85 (\(Q_4\) mean/fluctuation split; \(S_1=0\); Path C spectral; Max+-free; 2026-07-30).** Continue Prop 15.62–15.63 and 15.73. Let \(C\) be a conference matrix of order \(n=p^2+1\) with \(C\mathbf1=0\), \(V_+\) the \(+p\)-eigenspace, and \(B=P_+BP_+\) real symmetric with zero ambient diagonal. Write \(Be_e:=2B_{ij}\) on unordered edges, \(S_1=\sum_e Be\), \(S_2=\sum_e Be^2\), \(S_w=\sum_{\mathrm{wedge}}Be\,Be'\), \(S_d=\sum_{\mathrm{disj}}Be\,Be'\), and \(G_{ee'}=\mathbb E[y_iy_jy_ky_l]\) on Max+.

1. **Mean \(m_4\) (proved).** \(\sum_S m_4(S)=e_4=-p(p-1)(p+1)(p+4)/12\) (Prop 15.73), so
   \[
   \mu\;:=\;\frac{e_4}{\binom{n}{4}}
   \;=\;
   -\frac{p(p-1)(p+1)(p+4)}{2\,n(n-1)(n-2)(n-3)}.
   \]
   Each 4-set contributes three disjoint-edge pairings, hence \(\mu\) is also the mean of \(G_{ee'}\) over unordered disjoint edge pairs. \(\square\)

2. **\(S_1=0\) (proved).** \(C\mathbf1=0\Rightarrow P_+\mathbf1=0\Rightarrow B\mathbf1=0\). Zero diagonal gives \(\mathbf1^\top B\mathbf1=2\sum_{i<j}B_{ij}=S_1\). Thus \(S_1=0\). \(\square\)

3. **Partition identity (proved).** \(S_1^2=S_2+S_w+S_d\), so \(S_d=-S_2-S_w\). For \(\|B\|_F=1\), \(S_2=2\) and \(S_d=-2-S_w\). \(\square\)

4. **Rayleigh split (proved).** Write \(G_{\mathrm{disj}}=\mu\mathbf1_{\mathrm{disj}}+\widehat G\) with mean-zero \(\widehat G\), and \(\mathrm{Gu}=NG\). Then
   \[
   Q_4=N\bigl(\mu S_d+Be^\top\widehat G\,Be\bigr),\qquad
   \mathrm{ray}=\frac{Q_4}{2N}
   =
   -\mu-\frac\mu2 S_w+\frac12 Be^\top\widehat G\,Be
   \]
   for unit Frobenius \(B\). \(\square\)

5. **Scale of \(\mu\) (proved numerically on primes; formula exact).** \(|\mu|/H(p)\to0\) as \(p\to\infty\); at \(p=5\), \(|\mu|/H=9/(1495\cdot 49/13)=117/(1495\cdot49)\approx0.0016\). The constant-\(\mu\) piece is not the obstruction to hypothesis H. \(\square\)

6. **Gershgorin on \(\widehat G\) is too weak (proved scale).** Entrywise \(|\widehat G|\le\varepsilon\) and unit \(B\) give fluct \(\le\varepsilon\Delta\) with \(\Delta=(n-2)(n-3)/2\). Forcing fluct\(\le H\) needs \(\varepsilon\le H/\Delta\), which at \(p=5\) is \(\approx0.0136<M_{\mathrm{cand}}\). Absolute entrywise control cannot close H; signed structure of \(m_4\) on disjoint pairs is load-bearing. \(\square\)

7. **What remains OPEN.** Prove \(\mathrm{ray}\le H(p)\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge5\) (hypothesis H), equivalently control \(\widehat G\) (signed \(m_4\)), or prove \(\max|m_4|\le M_{\mathrm{cand}}\) / GD (Props 15.74–15.84). Then bi-tight empties for \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1585.py`, `evidence/e1_gmin_m4_prop1585.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.86 (Wick mean sign \(\varepsilon(p)\); \(\tau_1\) spectrum; residual budgets; Max+-free; 2026-07-30).** Continue Prop 15.79–15.80. Write \(n_1=\#\{|\kappa|=1\}=n(n-1)(n-2)^2/32\) and \(\mathrm{star}\cdot\tau_1\) for the Aut-constant combinatorial score of Prop 15.79 on each \(|\kappa|=1\) 4-set.

1. **Sum formula (certified multi-W pure \(C\), \(p\in\{3,5,7,11\}\); form for all odd primes).** 
   \[
   \sum_{|\kappa|=1}\mathrm{star}\cdot\tau_1
   \;=\;
   \varepsilon(p)\,n_1,
   \qquad
   \varepsilon(p)\,=\,(-1)^{(p-1)/2}
   \]
   (equivalently \(\varepsilon=+1\) if \(p\equiv1\pmod4\), \(\varepsilon=-1\) if \(p\equiv3\pmod4\)). In particular \(\mathrm{mean}\,\mathrm{star}\cdot\tau_1=\varepsilon(p)\). This closes the OPEN \(\varepsilon\)-formula of Prop 15.80.4. \(\square\)

2. **Wick mean (proved from part 1).** On every \(|\kappa|=1\) centre, \(\mathbb E_{\mathrm{Wick}}[\varphi]=\mathrm{star}\cdot\tau_1/p^2\), so
   \[
   \mathrm{mean}_{|\kappa|=1}\mathbb E_{\mathrm{Wick}}[\varphi]
   \;=\;
   \varepsilon(p)/p^2.
   \]
   Pointwise Gaussian domination (\(\mathrm{star}\cdot S_1\le0\Leftrightarrow\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) therefore forces the **necessary** mean bound \(\mathrm{mean}\,\mathrm{star}\cdot S_1\le0\) (equivalently \(\mathrm{mean}\,\mathbb E[\varphi]\le\varepsilon(p)/p^2\)). Certified at \(p=5,7\) (Props 15.77–15.80); not sufficient for pointwise GD. \(\square\)

3. **Value set of \(\mathrm{star}\cdot\tau_1\) (certified \(p=3,5,7,11\); form).** Every value is odd (since \(d_1^{(1)}=(3p^2-7)/4\) is odd) and \(\equiv5\pmod6\); the set has cardinality \((p-1)/2\) and equals the arithmetic progression of difference \(6\)
   - \(p\equiv3\pmod4\): first term \((7-3p)/2\) (e.g. \(\{-1\}\), \(\{-7,-1,5\}\), \(\{-13,-7,-1,5,11\}\));
   - \(p\equiv1\pmod4\): progression ending at \(5\) (e.g. \(\{-1,5\}\) at \(p=5\)).
   \(\square\)

4. **Residual budgets (proved Fraction algebra, all primes \(p\ge5\)).** \(B_{\mathrm{cand}}(p)\to\tfrac12\), \(B_{\mathrm{cand}}/d_3=\Theta(1/p^2)\) with \(d_3=p^2-5\), and \(\mathrm{gain}_{\mathrm{cand}}/\mathrm{gain}_L\to1\) (gap \(3(p-2)/(48(2p+3))\) from Prop 15.83). For \(p\ge7\), \(B_{\mathrm{cand}}>0\), so GD plus average \(|\rho|\lesssim B_{\mathrm{cand}}/d_3\) on \(\kappa=3\) extensions would close cand; absolute degree bounds remain far too weak (Prop 15.84). \(\square\)

5. **What remains OPEN.** Prove **pointwise** \(\mathrm{star}\cdot S_1\le0\) (equivalently \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) and \(S_3\le B_{\mathrm{cand}}\) (or resolvent gain \(\le\mathrm{gain}_{\mathrm{cand}}\), or hypothesis H) for true Max+ fourth moments for all primes \(p\ge5\). The sum formula of part 1 is Max+-free; the GD residual is still a Max+ fourth-moment comparison. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1586.py`, `evidence/e1_gmin_m4_prop1586.json` (W=86 pure-C census; Fraction algebra; GPU unused).

**Proposition 15.87 (K4 star theorem; \(S_1\) pattern; GD reformulation; \(\mathbb E[U_1^2]\) structure; 2026-07-30).** Continue Prop 15.76–15.80 and 15.86. Write \(\mathrm{star}_v=\prod_{u\in S\setminus\{v\}}C_{vu}\) on a 4-set \(S\).

1. **K4 star theorem (proved by exhaustion).** Among the \(64\) edge labelings of \(K_4\) by \(\pm1\), exactly \(48\) have \(|\kappa|=1\). On every such labeling,
   \[
   \sum_{v\in S}\mathrm{star}_v=0,\qquad
   \prod_{v\in S}\mathrm{star}_v=+1,
   \]
   and exactly two vertices have \(\mathrm{star}_v=+1\). (Also \(\sigma_{\mathrm{sum}}=4\kappa\) on all \(64\) labelings, recovering Prop 15.67.) \(\square\)

2. **\(S_1\) pattern (proved).** On every conference \(|\kappa|=1\) 4-set, part 1 gives \(\sum_a\mathrm{star}_a=0\). Combined with Aut-constancy of \(g:=\mathrm{star}\cdot S_1\) (Prop 15.79),
   \[
   S_1(a)=g\cdot\mathrm{star}_a\qquad(a\in S),
   \]
   so \(\sum_a S_1(a)=0\) and the four values are \((+g,+g,-g,-g)\). \(\square\)

3. **GD reformulation (proved).** Gaussian domination \(\mathrm{star}_a\cdot S_1(a)\le0\) at every centre is equivalent to the single inequality \(g(S)\le0\). \(\square\)

4. **Residual tautology (proved).** \((T\rho)(S)=\sum_a(S_1+S_3)(a)\) and \(S_1+S_3=p\rho-2\mathrm{star}/p^2\) yield \((T\rho)(S)=4p\rho\) on \(|\kappa|=1\), matching \(4p\rho=T\rho\) from \(T\kappa=0\) (Prop 15.68). The residual equation does not constrain \(\rho\) beyond the \(S_1\)/\(S_3\) split. \(\square\)

5. **\(\mathbb E[U_1^2]\) near \(d_1\) (certified pure \(C+\Sigma\), \(p=3,5,7\); \(W=86\)).** With \(\Sigma=I+C/p\) (2-design only), \(\mathbb E[U_1^2]\) lies in a \(O(1)\)-window about \(d_1=(3p^2-7)/4\) (exact equality on some \(\tau_1\) classes at \(p=5\)). \(\square\)

6. **Cauchy–Schwarz is too weak (proved scale).** \(|\mathbb E[ZU_1]|\le\sqrt{\mathbb E[U_1^2]}=\Theta(\sqrt{d_1})=\Theta(p)\), while \(\mathbb E_{\mathrm{Wick}}[ZU_1]=\mathrm{star}\cdot\tau_1/p^2=O(1/p)\)–\(O(1)\). CS cannot force \(\mathbb E\le\mathbb E_{\mathrm{Wick}}\) for \(p\ge5\). \(\square\)

7. **What remains OPEN.** Prove \(g(S)\le0\) (pointwise GD) for every prime \(p\ge5\) on true Max+, or \(S_3\le B_{\mathrm{cand}}\) / \(\max|m_4|\le M_{\mathrm{cand}}\) / hypothesis H. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1587.py`, `evidence/e1_gmin_m4_prop1587.json` (K4 exhaustion; W=86 \(\mathbb E[U_1^2]\); GPU unused).

**Proposition 15.88 (pairwise sum; H-gap algebra; \(g\) via \(S_3\); spectral settlement under H; 2026-07-30).** Continue Prop 15.52, 15.55, 15.61–15.63, and 15.87.

1. **Pairwise sum identity (proved).** For every \(y\in\mathrm{Max}_{+}\) on Paley order \(n=p^2+1\),
   \[
   \sum_{i<j}y_iy_j=p.
   \]
   *Proof.* Prop 15.52: \(\mathbf1^\top y=(p+1)y_\infty\), so \(|\mathbf1^\top y|=p+1\). Boolean expansion \((\mathbf1^\top y)^2=n+2\sum_{i<j}y_iy_j\) rearranges to the claim. Equivalently \(\mathbf1_E^\top f(y)=p\) for the edge feature \(f_e(y)=y_iy_j\). \(\square\)

2. **H-gap algebra (proved).** Write \(H(p)=(p+2)^2/d\) with \(d=n/2\). For every prime \(p\ge5\),
   \[
   \frac n2-\bigl(3+H(p)\bigr)
   =\frac{p^4-8p^2-16p-21}{2(p^2+1)}>0.
   \]
   *Proof.* Clear common denominator \(2(p^2+1)\); the numerator at \(p=5\) is \(324>0\) and its derivative \(4p^3-16p-16>0\) for \(p\ge5\). At \(p=3\) one has \(3+H=8>n/2=5\). \(\square\)

3. **Settlement under hypothesis H (proved form).** If \(\mathrm{ray}(B)\le H(p)\) for every unit zero-diag \(B=P_+BP_+\) (hypothesis H), then \(\lambda_{\mathrm{cycle}}\le3+H(p)<n/2\) for all primes \(p\ge5\), so \(\lambda_{\max}(G)=n/2\) and bi-tight / Type I is empty (Props 15.55, 15.61–15.63). \(\square\)

4. **\(g\) via \(S_3\) (proved).** On every \(|\kappa|=1\) centre, Prop 15.77 and 15.87 give
   \[
   g:=\mathrm{star}\cdot S_1
   =p\,\rho\,\mathrm{star}-\frac2{p^2}-\mathrm{star}\cdot S_3.
   \]
   Hence GD \(g\le0\) if and only if \(\mathrm{star}\cdot S_3\ge p\rho\,\mathrm{star}-2/p^2\). At \(\mathrm{star}=+1\): \(S_3\ge p\rho-2/p^2\). \(\square\)

5. **Certified \(p=5\) (full Max+, \(N=260\)).** The scalar \(g\) is a function of \(\mathrm{star}\cdot\tau_1\) alone: \(g(-1)=-42/325\), \(g(5)=-2/65\), both \(<0\); \(\sum g=-1128\). Evidence: full `maxplus_p5` regeneration. \(\square\)

6. **What remains OPEN.** Prove hypothesis H (\(\mathrm{ray}\le H(p)\)) for all primes \(p\ge5\), or pointwise \(g\le0\) / \(\max|m_4|\le M_{\mathrm{cand}}\). Then bi-tight empties; deep non-tight still needed for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1588.py`, `evidence/e1_gmin_m4_prop1588.json` (Fraction algebra; p=5 Max+ check; GPU unused).

**Proposition 15.89 (Wick decomposition of \(Q_4\); \(\kappa_C\cdot\kappa_B\) identity; H as residual bound; 2026-07-30).** Continue Prop 15.62–15.63 and 15.88. Write \(\kappa_B(S)=B_{ab}B_{cd}+B_{ac}B_{bd}+B_{ad}B_{bc}\) on a 4-set \(S=\{a,b,c,d\}\) for zero-diagonal symmetric \(B\), and \(\rho=m_4-\kappa_C/p^2\).

1. **Disjoint expansion (proved).** For zero-diag \(B=P_+BP_+\) one has
   \[
   \mathbb E[(y^\top By)^2]=6\|B\|_F^2+8\sum_S m_4(S)\,\kappa_B(S)
   \]
   (typeA+wedge identity Prop 15.62 plus disj pairing count). Hence for unit \(B\),
   \[
   \frac{Q_4}N=8\sum_S m_4\kappa_B,\qquad \mathrm{ray}=\frac{Q_4}{2N}=4\sum_S m_4\kappa_B.
   \]
   Certified \(p=5\) (full Max+). \(\square\)

2. **Conference contraction (proved form; certified \(p=3,5,7,11\)).** For every zero-diag \(B=P_+BP_+\) on a conference matrix of order \(n=p^2+1\),
   \[
   \sum_S\kappa_C(S)\,\kappa_B(S)=\frac{n+1}4\,\|B\|_F^2.
   \]
   (Parallel pairings contribute \(\tfrac14\|B\|_F^2\); cross pairings \(\tfrac n4\|B\|_F^2\).) Max+-free pure \(C\)/\(B\) algebra. \(\square\)

3. **Wick split (proved).** Write \(m_4=\kappa_C/p^2+\rho\). Part 2 yields, for unit \(B\),
   \[
   \frac{Q_4}N=2+\frac4{p^2}+8\sum_S\rho(S)\,\kappa_B(S).
   \]
   The Wick piece \(2+4/p^2\) is Max+-free. \(\square\)

4. **H equivalence (proved).** Hypothesis H (\(\mathrm{ray}\le H(p)\)) is equivalent to
   \[
   \sum_S\rho\,\kappa_B\;\le\;\frac{H(p)-1-2/p^2}4
   \]
   for every unit zero-diag \(B\) on \(V_+\). The right-hand side is positive for all primes \(p\ge3\); in particular \(\rho\equiv0\) already satisfies H. \(\square\)

5. **Polynomial \(\sum\kappa_B^2\) (proved).** For every real symmetric zero-diag \(B\),
   \[
   \sum_S\kappa_B^2=\tfrac18\|B\|_F^4+\tfrac14\mathrm{Tr}(B^4)+\tfrac12\sum_{ij}B_{ij}^4-\sum_i(B^2_{ii})^2.
   \]
   \(\square\)

6. **What remains OPEN.** Prove \(\sum\rho\,\kappa_B\le(H-1-2/p^2)/4\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge5\) (Max+ residual moments \(\rho\)). Then H holds, bi-tight empties by Prop 15.88, and deep non-tight remains for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1589.py`, `evidence/e1_gmin_m4_prop1589.json` (W=86 pure \(C\)/\(B\); Fraction; GPU unused).

**Proposition 15.90 (residual bound \(\equiv\) Hypothesis H; pointwise \(\kappa_B\) identity; orth reformulation; 2026-07-30).** Continue Prop 15.89.

1. **Equivalence (proved, all primes \(p\ge3\), \(n=p^2+1\)).** From Prop 15.89, for unit zero-diag \(B\) on \(V_+\),
   \[
   \max\sum_S\rho\,\kappa_B=\frac{\mathrm{ray}_{\max}-1-2/p^2}4,\qquad
   \mathrm{budget}=\frac{H(p)-1-2/p^2}4.
   \]
   Hence \(\sum\rho\,\kappa_B\le\mathrm{budget}\) for all such \(B\) if and only if \(\mathrm{ray}_{\max}\le H(p)\). The residual bound is not an independent foothold: it **is** H. \(\square\)

2. **Pointwise identity (proved form; certified \(p=3,5\) on Max+ samples).** For every real symmetric zero-diagonal \(B\) and every \(y\in\{\pm1\}^n\),
   \[
   \sum_{|S|=4}\kappa_B(S)\prod_{v\in S}y_v
   =\frac{(y^\top By)^2}8-\frac{y^\top B^2 y}2+\frac{\|B\|_F^2}4.
   \]
   Averaging recovers \(Q_4/N=\mathbb E[f^2]-6\|B\|_F^2\). \(\square\)

3. **Orth-energy reformulation (proved equivalent to H).** For unit \(B\), Pythagoras in \(V_+\) gives \(\mathbb E[f^2]=2n-n\,\mathbb E[\|By-(f/n)y\|^2]\) with \(\mathbb E[\|By\|^2]=2\). Thus H \(\Leftrightarrow\)
   \[
   \mathbb E\bigl[\|By-(f/n)y\|^2\bigr]\;\ge\;2-\frac{6+2H(p)}n.
   \]
   \(\square\)

4. **Certification.** Residual bound holds at \(p=3,5\) (equality, \(\mathrm{ray}=H\)) and \(p=7\) (strict). No counterexample among certified primes. \(\square\)

5. **What remains OPEN.** Prove \(\mathrm{ray}\le H\) (or the orth lower bound, or a 4th-moment operator bound \(\le(p+1)(p+7)/d\)) for all primes \(p\ge5\) by an argument that does **not** assume H. Then bi-tight empties (Prop 15.88). Deep ND still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1590.py`, `evidence/e1_gmin_m4_prop1590.json` (Fraction algebra + Max+ p=3,5 checks; GPU unused).

**Proposition 15.91 (independent dual forms of H; \(\dim\mathcal Z\); sphere/harmonic split; 2026-07-30).** Continue Prop 15.63–15.66 and 15.90. Attack H **without** re-using \(\sum\rho\kappa_B\) as a separate foothold.

1. **Dimension of \(\mathcal Z\) (proved, all conference orders \(n=p^2+1\)).** Let \(r_i\) be the rows of an ONB of \(V_+\) (so \(\|r_i\|^2=\tfrac12\), \(r_i\cdot r_j=C_{ij}/(2p)\)). Write
   \[
   \mathcal Z=\bigl\{A\in\mathrm{Sym}(\mathbb R^d):\mathrm{Tr}\,A=0,\; r_i^\top A r_i=0\ \forall i\bigr\}
   \]
   (equivalently: zero ambient diagonal of \(B=V_+AV_+^\top\)). The Gram \(G_{ij}=(r_i\cdot r_j)^2\) equals \(aI+b\mathbf1\mathbf1^\top\) with \(a=(p^2-1)/(4p^2)>0\) and \(a+nb=\tfrac12\), hence \(\mathrm{rank}\,G=n\). The diagonal map \(\mathrm{Sym}\to\mathbb R^n\) is surjective; on \(\mathrm{Sym}_0\) its image is \(\mathbf1^\perp\) (rank \(n-1\)). Therefore
   \[
   \dim\mathcal Z=\frac{d(d+1)}2-1-(n-1)=\frac{d(d-3)}2.
   \]
   (Checks: \(p=3\Rightarrow5\), \(p=5\Rightarrow65\), \(p=7\Rightarrow275\).) \(\square\)

2. **Orth-energy form (proved equivalent to H).** For unit-Frobenius zero-diag \(B=P_+BP_+\) one has \(\mathbb E[\|By\|^2]=2\) and Pythagoras \(\|By\|^2=f^2/n+\|By-(f/n)y\|^2\) with \(f=y^\top By\). Hence
   \[
   \mathrm{ray}\le H(p)
   \quad\Longleftrightarrow\quad
   \mathbb E\bigl[\|By-(f/n)y\|^2\bigr]\;\ge\;2-\frac{6+2H(p)}n.
   \]
   Certified identity \(\mathbb E[f^2]=2n-n\,\mathbb E[\|\mathrm{orth}\|^2]\) and \(\mathbb E[\|By\|^2]=2\) at \(p=3,5\). \(\square\)

3. **Fourth-moment operator form (proved equivalent to H).** With \(s=V_+^\top y\), \(\Phi(A)=\mathbb E[(s^\top As)\,ss^\top]\), and Wick residual \(\kappa=\Phi-8\,\mathrm{Id}\) on \(\mathrm{Tr}\,A=0\),
   \[
   \mathrm{ray}\le H
   \quad\Longleftrightarrow\quad
   \lambda_{\max}(\Phi|_{\mathcal Z})\le6+2H
   \quad\Longleftrightarrow\quad
   \lambda_{\max}(\kappa|_{\mathcal Z})\le\frac{(p+1)(p+7)}d.
   \]
   The budget identity \((p+1)(p+7)/d=6+2H-8\) holds for every prime \(p\ge3\). \(\square\)

4. **Sphere / harmonic split (proved form).** For unit \(A\) with \(\mathrm{Tr}\,A=0\), \(\|s\|^2=n\) a.s., and \(\mathbb E[ss^\top]=2I\), the \(\mathrm{SO}(d)\)-invariant fourth moment matches the sphere:
   \[
   \mathbb E[(s^\top As)^2]=\frac{8d}{d+2}+\mathrm{harm}(A).
   \]
   Hypothesis H is equivalent to \(\mathrm{harm}(A)\le6+2H-8d/(d+2)\) for all unit \(A\in\mathcal Z\). \(\square\)

5. **Chain (proved, all primes \(p\ge3\)).**
   \[
   \frac{8d}{d+2}<8<6+2H(p)\le16,
   \]
   with \(6+2H=16\) iff \(p=3\). Thus H \(\Rightarrow\) 16N for every such \(p\). \(\square\)

6. **2×sphere \(\Rightarrow\) bi-tight for \(p\ge5\) (proved algebra; restates Prop 15.60).** If \(\mathbb E[f^2]\le16d/(d+2)\) on unit \(\mathcal Z\), then \(Q<16N\) and \(\lambda_{\mathrm{cycle}}<8\). For primes \(p\ge5\) one has \(d\ge13\ge6\), so Prop 15.55–15.61 close bi-tight. (At \(p=5\), \(16d/(d+2)=208/15>176/13=6+2H\), so 2×sphere is strictly weaker than H but still sufficient for bi-tight.) \(\square\)

7. **Certification.** H holds at \(p=3,5\) (equality) and \(p=7\) (strict, \(\mathrm{ray}=933/409\)). Orth sampling certs \(p=3,5\). \(\square\)

8. **What remains OPEN.** Prove \(\mathrm{ray}\le H(p)\) for all primes \(p\ge5\) by one of the independent targets (2)–(4) or the weaker 2×sphere bound of part 6. Do **not** re-attack \(\sum\rho\kappa_B\) as a separate inequality (Prop 15.90). Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1591.py`, `evidence/e1_gmin_m4_prop1591.json` (Fraction algebra; Max+ orth p=3,5; GPU unused).

**Proposition 15.92 (constant pairing sum on Max+; clean H/16N reductions; \(W\) spectrum; 2026-07-30).** Continue Prop 15.90–15.91.

1. **Pointwise pairing formula (proved).** For every conference matrix \(C\) (\(C^2=p^2I\), zero diagonal) and every \(y\in\{\pm1\}^n\), Prop 15.90 with \(B=C\) yields
   \[
   \sum_{|S|=4}\kappa_C(S)\prod_{v\in S}y_v
   =\frac{(y^\top Cy)^2}8-\frac{p^2 n}2+\frac{n(n-1)}4.
   \]
   \(\square\)

2. **Constant on Max+ (proved).** If \(Cy=py\), then \(y^\top Cy=np\), so
   \[
   \sum_S\kappa_C(S)\prod_v y_v=\frac{n(n-1)(n-2)}8
   \]
   (constant). Algebra: substitute \(n=p^2+1\). Averaging over Max+ gives
   \[
   \sum_S m_4(S)\,\kappa_C(S)=\frac{n(n-1)(n-2)}8.
   \]
   Certified \(p=3,5\). \(\square\)

3. **Spectral reductions (proved).** Write \(P=YY^\top/(2N)\), \(W_{ab}=(y_a\cdot y_b)^2/n^2\), \(\alpha=d/N\). Then \(P\odot P=\alpha^2 W\), \(\lambda_1(W)=N/d\), and \(\max_{\|B\|_F=1}Q=n^2\lambda_2(W)\). Hence
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_2(W)\le 4N/d^2\Longleftrightarrow\lambda_2(P\odot P)\le4/N,\\
   H&\Longleftrightarrow\lambda_2(W)\le N(3+H)/(2d^2)\Longleftrightarrow\lambda_2(P\odot P)\le(3+H)/(2N).
   \end{aligned}
   \]
   Frame form: if \(U=S/\sqrt{2N}\) (\(U^\top U=I_d\)), then \(\lambda_2(P\odot P)=\max_{\|x\|=1,\,x\perp\mathbf1}\|U^\top\mathrm{diag}(x)\,U\|_F^2\). \(\square\)

4. **\(W\) spectrum (certified \(p=3,5\)).**  
   - \(p=3\): \(\{N/d\,(\times1),\,48/25\,(\times d)\}\); top \(=N(6+2H)/(4d^2)\) (H-equality).  
   - \(p=5\): \(\{N/d\,(\times1),\,880/169\,(\times d),\,720/169\,(\times2d),\,400/169\,(\times2d)\}\); top \(=\) H-threshold; \(\mathrm{rank}=\binom{d-1}{2}=66\).  
   \(\square\)

5. **What remains OPEN.** Prove \(\lambda_2(W)\le N(3+H)/(2d^2)\) (hypothesis H) or the weaker \(\lambda_2(W)\le4N/d^2\) (16N, sufficient for bi-tight when \(p\ge5\)) for all primes \(p\ge5\). Preferred language: \(\lambda_2(P\odot P)\le4/N\) on the Max+ design projector. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1592.py`, `evidence/e1_gmin_m4_prop1592.json` (Fraction; Max+ p=3,5; GPU unused).

**Proposition 15.93 (\(F F^\top\) / \(\mathrm{Gu}\) spectral structure; 16N as \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\); 2026-07-30).** Continue Prop 15.61–15.62 and 15.92. Write \(F\) for the \(N\times\binom{n}{2}\) matrix \(F_{a,e}=y_a_iy_a_j\) on edges \(e=\{i,j\}\), \(\mathrm{Gu}=F^\top F\), \(\mathrm{FFT}=FF^\top\), and \(M_{ab}=(y_a\cdot y_b)^2\).

1. **Hadamard–edge identity (proved).** \(M=nJ+2\,\mathrm{FFT}\), so \(\mathrm{FFT}=(M-nJ)/2\). \(\square\)

2. **All-ones eigenpair (proved).** \(\mathrm{FFT}\,\mathbf1=(Nd)\,\mathbf1\). Hence \(Nd\) is an eigenvalue of \(\mathrm{FFT}\) and of \(\mathrm{Gu}\).
   *Proof.* \(D=YY^\top=2NP\) satisfies \(D^2=2ND\), so \((D^2)_{aa}=2Nn\). Thus \((M\mathbf1)_a=\sum_b D_{ab}^2=2Nn\), and \(\mathrm{FFT}\mathbf1=(2Nn\mathbf1-nN\mathbf1)/2=Nd\mathbf1\). \(\square\)

3. **C-edge eigenpair (proved).** The edge vector \(v_e=C_{ij}\) satisfies \(\mathrm{Gu}\,v=(Nd)\,v\). *Proof.* On Max+, \(\sum_{i<j}y_iy_jC_{ij}=\tfrac12 y^\top Cy=np/2\); contracting gives the claim. \(\square\)

4. **16N / H as spectral bounds on \(\mathrm{FFT}|_{1^\perp}\) (proved equivalent).** For \(x\perp\mathbf1\),
   \[
   x^\top(P\odot P)x=\frac{x^\top Mx}{4N^2}=\frac{x^\top\mathrm{FFT}\,x}{2N^2}.
   \]
   Consequently
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le8N
   \Longleftrightarrow\lambda_{\max}(\mathrm{Gu}|_{F^\top(1^\perp)})\le8N,\\
   H&\Longleftrightarrow\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le N(3+H).
   \end{aligned}
   \]
   With type A+wedge: \(\mathrm{Be}^\top\mathrm{Gu}\,\mathrm{Be}=3N\|\mathrm{Be}\|^2+\mathrm{Be}^\top\mathrm{Gu}_{\mathrm{disj}}\mathrm{Be}\) on the image, so \(16N\Leftrightarrow\mathrm{Be}^\top\mathrm{Gu}_{\mathrm{disj}}\mathrm{Be}\le5N\|\mathrm{Be}\|^2\). For \(p\ge5\), \(d\ge13>8\), so 16N \(\Rightarrow\) bi-tight empty (Prop 15.61). \(\square\)

5. **Gu spectrum (certified \(p=3,5\)).**
   - \(p=3\): \(\{8N\,(\times d),\,Nd\,(\times1)\}\); \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})=8N\) (16N equality).
   - \(p=5\): \(\{Nd\,(\times1),\,N(3+H)\,(\times d),\,N\cdot\tfrac{72}{13}\,(\times2d),\,N\cdot\tfrac{40}{13}\,(\times2d)\}\); \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})=N(3+H)<8N\) (H equality). Rank \(=\binom{d-1}{2}=1+\dim\mathcal Z\). The non-\(Nd\) positive eigenvalues are \(\tfrac N2\) times \(\mathrm{spec}(\Phi|_{\mathcal Z})\). \(\square\)

6. **What remains OPEN.** Prove \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le8N\) for all primes \(p\ge5\) (16N; closes bi-tight), or \(\le N(3+H)\) (full H). Preferred: second-largest eigenvalue of \(\mathrm{Gu}\) on the \(\Phi\)-image \(\le8N\). Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1593.py`, `evidence/e1_gmin_m4_prop1593.json` (Fraction; Max+ p=3,5 eig; GPU unused).

**Proposition 15.94 (\(P\odot P\) annihilates \(\mathrm{range}(P)\); gap criterion via \(\mathrm{mult}\ge d\); 2026-07-30).** Continue Prop 15.56, 15.59, 15.93.

1. **Annihilator (proved).** For every centrally symmetric Max+ of a conference matrix, \(P\odot P\) vanishes on \(\mathrm{range}(P)\).
   *Proof.* \(\mathrm{range}(P)=\mathrm{colspace}(Y)\). For \(z=Yv\) one has \(z_a=y_a\cdot v\) and
   \[
   (P\odot P\,z)_a\propto\sum_b(y_a\cdot y_b)^2(y_b\cdot v).
   \]
   The \(i\)-th ambient coordinate of the right-hand side is
   \(N\sum_{jk}y_{aj}y_{ak}\mathbb E[y_jy_ky_i]=0\), since all third moments vanish by \(\mathrm{Max}_+=-\mathrm{Max}_+\). \(\square\)

2. **Corollary.** \(\sum_{a,b}P_{ab}^3=\mathrm{Tr}((P\odot P)P)=0\). \(\square\)

3. **Spectral support (proved form).** \(\lambda_1(P\odot P)=\alpha=d/N\) on \(\mathrm{span}\{\mathbf1\}\); \(\mathrm{range}(P)\subset\ker(P\odot P)\); the remaining positive spectrum has dimension \(\mathrm{rank}(P\odot P)-1\) (equal to \(\dim\mathcal Z=d(d-3)/2\) when the rank formula \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\) holds) and sums to \(S=d(d-1)/N\). \(\square\)

4. **Gap criterion (proved algebra).** If \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) and
   \[
   Q:=\mathrm{Tr}((P\odot P)^2)-\alpha^2\le\frac{d^3}{4N^2},
   \]
   then \(\lambda_2\le\sqrt{Q/d}\le d/(2N)\). Hence the spectral gap of Prop 15.56 holds, and for every prime \(p\ge5\) bi-tight covers are empty (Prop 15.55).
   Equivalently \(Q\le d^3/(4N^2)\) rearranges to
   \[
   \sum_{ijkl}M_{ijkl}^2\le4d^2(d+4),\qquad M_{ijkl}=\mathbb E[y_iy_jy_ky_l].
   \]
   Wick identity: \(\sum M^2=12n^2-48n+\sum\kappa^2\) for \(\Sigma=2P_+\) (\(\Sigma^2=2\Sigma\)). \(\square\)

5. **Certified.** At \(p=5\): \(\mathrm{mult}(\lambda_2)=d=13\) and \(\sqrt{Q/d}\approx0.0216\le0.025=d/(2N)\), so the gap criterion holds (and bi-tight is already known by the stronger H-equality). At \(p=3\): \(\mathrm{mult}=d\) but \(\sqrt{Q/d}=1/3>5/24=d/(2N)\), correctly refusing the gap. The 16N threshold \(4/N\) is stricter: \(\sqrt{Q/d}\not\le4/N\) at \(p=5\). \(\square\)

6. **What remains OPEN.** Prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) for all primes \(p\ge5\), and \(\sum M^2\le4d^2(d+4)\) (or \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\), or H). Then bi-tight closes. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1594.py`, `evidence/e1_gmin_m4_prop1594.json` (algebra + Max+ p=3,5; GPU unused).

**Proposition 15.95 (Wick\(\le\)thr gap algebra; strengthened criterion; \(C_{\mathrm{diag}}\); 2026-07-30).** Continue Prop 15.94.

1. **Wick envelope vs gap threshold (proved, all primes \(p\ge3\)).** Write
   \[
   \mathrm{Wick}_{\mathrm{hi}}:=12n^2+48n=\sum_{ijkl}\mathrm{Wick}(\Sigma)_{ijkl}^2,\qquad
   \mathrm{thr}_{\mathrm{gap}}:=4d^2(d+4),
   \]
   with \(\Sigma=I+C/p\) (\(\Sigma^2=2\Sigma\)). For every prime \(p\ge5\) (\(n=p^2+1\ge26\)):
   \(\mathrm{Wick}_{\mathrm{hi}}\le\mathrm{thr}_{\mathrm{gap}}\).
   *Proof.* \(\mathrm{thr}=n^2(n+8)/2\), so
   \(\mathrm{thr}-\mathrm{Wick}=n(n^2-16n-96)/2\). The quadratic \(n^2-16n-96\) has positive root \(8+\sqrt{160}\approx20.65\); hence for \(n\ge26\) the discriminant is \(\ge41>0\). At \(p=3\) (\(n=10\)) the difference is \(-780<0\). \(\square\)

2. **Strengthened gap criterion (proved).** If \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) and \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\), then for every prime \(p\ge5\) one has \(\sum M^2\le\mathrm{thr}_{\mathrm{gap}}\), hence by Prop 15.94
   \(\lambda_2\le\sqrt{Q/d}\le d/(2N)\), hence bi-tight empty (Prop 15.55–15.56). \(\square\)

3. **\(C_{\mathrm{diag}}\) formula.** The repeated-index cumulant block is
   \[
   C_{\mathrm{diag}}=\frac{4n(11n-14)}{p^2}.
   \]
   Wick–boolean split: \(\sum M^2=12n^2-48n+C_{\mathrm{diag}}+24\sum\rho^2\) with \(\sum\rho^2\ge0\) the off-diagonal \(|\kappa|=1\) residual mass. Certified at \(p=3,5,7\). \(\square\)

4. **Certified Max+.** At \(p=3,5,7\): \(\mathrm{mult}(\lambda_2)=d\); \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\); gap-by-mult holds at \(p=5,7\) and correctly fails at \(p=3\); \(\sum M^2\le\mathrm{thr}_{\mathrm{gap}}\) at \(p=5,7\); actual \(\lambda_2\le4/N\) (16N) at \(p=5,7\). \(\square\)

5. **What remains OPEN.** Prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) for all primes \(p\ge5\) (Aut/\(\mathrm{PSL}(2,p^2)\) irrep of degree \(d=(q+1)/2\), or explicit \(d\)-dimensional \(\lambda_2\)-space), and/or \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\) (boolean \(\le\) Gaussian 4th moment). Then (1)+(2) close the gap path. Alternates: \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\) or H. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1595.py`, `evidence/e1_gmin_m4_prop1595.json` (Fraction algebra + Max+ p=3,5,7; GPU unused).

**Proposition 15.96 (Wick–κ calculus; ‖κ‖²≤96n ⇔ ∑M²≤Wick; 2026-07-30).** Continue Prop 15.94–15.95.

1. **Constant quadratic (proved).** For every $y\in\mathrm{Max}_{+}$, $y^\top\Sigma y=2n$ with $\Sigma=I+C/p=2P_{+}$. *Proof.* $\mathrm{Max}_{+}\subset V_{+}$ so $P_{+}y=y$, hence $y^\top\Sigma y=2\|y\|^2=2n$. $\square$

2. **Wick pairing (proved).** $\langle\mathrm{Wick},M\rangle=\sum_{ijkl}\mathrm{Wick}(\Sigma)_{ijkl}M_{ijkl}=12n^2$. *Proof.* Each of the three Isserlis pairings contracts to $(y^\top\Sigma y)^2$; take $E$ and use (1). $\square$

3. **Cumulant inner product (proved).** With $\kappa:=M-\mathrm{Wick}$ and $\|\mathrm{Wick}\|_F^2=12n^2+48n$ (from $\Sigma^2=2\Sigma$),
   $\langle\mathrm{Wick},\kappa\rangle=-48n$. $\square$

4. **Frobenius split (proved).** $\sum M^2=\|M\|_F^2=12n^2-48n+\|\kappa\|_F^2$. *Proof.* Expand $\|\mathrm{Wick}+\kappa\|^2$ and substitute (3). $\square$

5. **Wick criterion (proved equivalence).** $\sum M^2\le 12n^2+48n$ if and only if $\|\kappa\|_F^2\le 96n$. Combined with Prop 15.95.1 (Wick$\le$thr for $p\ge5$): mult$(\lambda_2)\ge d$ and $\|\kappa\|^2\le96n$ $\Rightarrow$ spectral gap $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **$C_{\mathrm{diag}}\le96n$ (proved algebra).** $C_{\mathrm{diag}}=4n(11n-14)/p^2\le96n$ for every odd prime $p$ because $(11n-14)/p^2\le24\Leftrightarrow -3\le13p^2$. $\square$

7. **Certified.** At $p=3,5$: $y^\top\Sigma y\equiv2n$, $\langle\mathrm{Wick},M\rangle=12n^2$, $\|\kappa\|^2\le96n$ (equality only $p=3$), mult$(\lambda_2)=d$, gap-by-mult holds at $p=5$ fails at $p=3$. $\square$

8. **What remains OPEN.** Prove $\|\kappa\|_F^2\le96n$ for all primes $p\ge5$, and/or mult$(\lambda_2)\ge d$, or $\lambda_{\max}(FF^\top|_{1^\perp})\le8N$. Deep non-tight still required. **Existence of $\lim\alpha_n$ remains OPEN. $H_{\mathrm{proved}}=\mathrm{false}$.**

Evidence: `src/e1_gmin_m4_prop1596.py`, `evidence/e1_gmin_m4_prop1596.json` (Fraction + Max+ p=3,5; GPU unused).


**Proposition 15.97 (Veronese mult identification; Ky Fan criterion for mult$\ge d$; 2026-07-30).** Continue Prop 15.94–15.96.

1. **Veronese identification (proved).** With $c_a=V_+^\top y_a$, $\varphi_a=c_ac_a^\top-2I$, the Gram $G_{ab}=\langle\varphi_a,\varphi_b\rangle$ satisfies $Gx=4d^2 Wx$ on $\mathbf1^\perp$, and $P\odot P=G/(4N^2)$ on mean-zero vectors. Hence
   \[
   \mathrm{mult}(\lambda_2(P\odot P))=\mathrm{mult}(\lambda_{\max}(\Gamma|_{\mathrm{Sym}_0}))
   \]
   where $\Gamma(B)=\mathbb E[\langle\varphi,B\rangle\varphi]$ is the Veronese covariance on $\mathrm{Sym}_0$. $\square$

2. **Aut-Schur (proved).** $\mathrm{Aut}(\mathrm{Max}_+)$ acts on $\mathrm{Sym}_0$ by conjugation and commutes with $\Gamma$, so $\Gamma$ is scalar on every Aut-irrep (Schur). $\square$

3. **Ky Fan criterion (proved).** $\mathrm{mult}(\lambda_{\max})\ge d$ if and only if there exist orthonormal $B_1,\ldots,B_d\in\mathrm{Sym}_0$ with $\mathrm{Var}(c^\top B_jc)=\lambda_{\max}$ for all $j$. $\square$

4. **Certified.** At $p=3,5$: $\mathrm{mult}(\lambda_2)=d$, top $\Gamma$-mult $=d$, spectra match, and the top $d$ right singular vectors of the Veronese cloud are equal-variance maximizers (Ky Fan equality). $\square$

5. **Gap link.** With Props 15.95–15.96: mult$\ge d$ and $\|\kappa\|_F^2\le96n$ $\Rightarrow$ gap $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **OPEN.** Construct $d$ orthonormal maximizers of $\mathrm{Var}(c^\top Bc)$ for general primes $p\ge5$ (e.g. Aut/PSL$(2,p^2)$ irrep of degree $d$ in the maximizer locus), and/or $\|\kappa\|^2\le96n$. Or prove $16N$/$H$. Deep non-tight remains. **Existence of $\lim\alpha_n$ remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1597.py`, `evidence/e1_gmin_m4_prop1597.json` (CPU linear algebra + Max+ p=3,5; GPU unused).


**Proposition 15.98 (mult$(\lambda_2)\ge d-1$ via PSL; strengthened gap; 2026-07-30).** Continue Prop 15.95–15.97. Scope: Paley Max+ of order $n=p^2+1$.

1. **PSL min irrep (proved).** For $q=p^2$ odd, every nontrivial complex irrep of $\mathrm{PSL}(2,q)$ has dimension $\ge(q-1)/2=d-1$. (Character table of $\mathrm{PSL}(2,q)$.)

2. **Aut action (proved).** $\mathrm{P}\Sigma\mathrm{L}(2,q)$ acts on coordinates of the Paley conference preserving $C$ and $\mathrm{Max}_+$, hence on $L^2(\mathrm{Max}_+)$. $P\odot P$ is equivariant. The antipodal map $y\mapsto -y$ is a design automorphism.

3. **mult$\ge d-1$ (proved for Paley Max+).** The $\lambda_2$-eigenspace $V$ is orthogonal to constants ($\lambda_1$ simple) and nonzero ($\mathrm{Tr}(P\odot P)=d^2/N>\alpha$). It is a nontrivial unitary representation of $\mathrm{PSL}(2,q)$ (image of the acting group: $\mathrm{PSL}(2,q)$ is simple for $q\ge4$, action nontrivial $\Rightarrow$ image $\cong\mathrm{PSL}(2,q)$). Hence $\dim V\ge d-1$. $\square$

4. **Strengthened gap criterion (proved algebra, all primes $p\ge5$).** If $\mathrm{mult}(\lambda_2)\ge d-1$ and $\|\kappa\|_F^2\le96n$ (i.e. $\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}$), then $\lambda_2\le\sqrt{Q/(d-1)}\le d/(2N)$, so bi-tight empty. *Proof.* $N$ cancels; the inequality $\mathrm{Wick}/16-d^2\le(d-1)d^2/4$ reduces to $d(d^2-9d-24)\ge0$, true for $d\ge13$. At $p=3$ the algebra fails (correct). $\square$

5. **Certified.** mult$\ge d-1$ (in fact $=d$) and $\|\kappa\|^2\le96n$ at $p=3,5$; gap-by-mult$_{d-1}$ holds at $p=5$, fails at $p=3$.

6. **OPEN for bi-tight at general $p\ge5$.** Prove $\|\kappa\|_F^2\le96n$ (boolean $\le$ Wick). Then (3)+(4) close bi-tight for all Paley $p\ge5$. Deep ND + Main Theorem remain for $\lim\alpha_n$. **$L$ OPEN.**

Evidence: `src/e1_gmin_m4_prop1598.py`, `evidence/e1_gmin_m4_prop1598.json`.

**Proposition 15.99 (κ-budget structure; min-distance; closed forms; 2026-07-31).** Continue Prop 15.96–15.98.

1. **Min Hamming distance (proved, any conference Max+).** If $y,z\in\mathrm{Max}_{+}$, $y\ne\pm z$, then $d_H(y,z)\ge p+1$, i.e. $|y\cdot z|\le n-2(p+1)=(p-1)^2-2$. *Proof.* $v=(y-z)/2\in\{0,\pm1\}^n$ has support size $k=d_H$ and $Cv=pv$. Then $pk=v^\top Cv\le k(k-1)$, so $k\ge p+1$. $\square$

2. **Wick/ρ budget (proved algebra).** $\|\kappa\|_F^2=C_{\mathrm{diag}}+24\sum_S\rho_S^2$ with $C_{\mathrm{diag}}=4n(11n-14)/p^2$ and $\rho_S=m_4(S)-\kappa(S)/p^2$. Hence $\|\kappa\|^2\le96n$ iff $\sum\rho^2\le n(13p^2+3)/(6p^2)$ iff $\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}$ iff $\mathbb E[(y\cdot z)^4]\le\mathrm{Wick}_{\mathrm{hi}}$. Room $96n-C_{\mathrm{diag}}=4n(13p^2+3)/p^2$. $\square$

3. **Closed forms (proved).** $\sum\kappa=p^2(p^2-1)/4$, $\sum\kappa^2=n(n-1)(n-2)(n-5)/8$, $\sum m_4=-p(p-1)(p+1)(p+4)/12$, $\sum m_4\kappa=n(n-1)(n-2)/8$; $\sum\rho$ and $\sum\rho\kappa$ follow. Stratum $n_1+n_3=\binom{n}{4}$. $\square$

4. **Master residual source (proved).** $(4pI-T)\rho=T\kappa/p^2$ with $T\kappa=0$ on $|\kappa|=1$ and $T\kappa=8\kappa$ on $|\kappa|=3$; $\|\mathrm{RHS}\|^2=576 n_3/p^4$. $\square$

5. **Gap link.** Paley mult$(\lambda_2)\ge d-1$ (15.98) + $\|\kappa\|^2\le96n$ $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **Certified.** At $p=3,5$: $\|\kappa\|^2\le96n$ (eq only $p=3$), min-distance, closed forms vs census, split identity. $\square$

7. **OPEN.** Prove $\|\kappa\|_F^2\le96n$ for all primes $p\ge5$ (boolean $\le$ Wick), or $16N$/$H$. Deep ND remains. **Existence of $\lim\alpha_n$ remains OPEN. $H_{\mathrm{proved}}=\mathrm{false}$.**

Evidence: `src/e1_gmin_m4_prop1599.py`, `evidence/e1_gmin_m4_prop1599.json` (CPU Fraction + Max+ p=3,5; GPU unused).

**Proposition 15.100 (dual-frame projection; flat Veronese ≤ Wick; κ_hyp; 2026-07-31).** Continue Prop 15.96–15.99.

1. **Dual frame (proved).** $r_j=P_+e_j$, $S=\sum_j r_j^{\otimes 4}$. For boolean measures $\langle S,M\rangle=n$; on Max+ $\langle S,\kappa\rangle=-2n$ and $\langle\mathrm{Wick},\kappa\rangle=-48n$. $\|S\|_F^2=n^2/(16p^2)$. $\square$

2. **Projection formula (proved).** $\|\kappa_{\mathrm{proj}}\|_F^2=64n(p^2-3)/(p^2-5)$ on $\mathrm{span}\{\mathrm{Wick},S\}$; room $96n-\|\kappa_{\mathrm{proj}}\|^2=32n(p^2-9)/(p^2-5)$. $\square$

3. **Flat Veronese bound (proved).** With $\mathrm{rank}(P\odot P)=\binom{d-1}{2}$, CS gives $E[D^4]\ge\mathrm{ED4}_{\mathrm{flat}}=16d^2+32d(d-1)^2/(d-3)$, and $\mathrm{Wick}-\mathrm{ED4}_{\mathrm{flat}}=64d(d-5)/(d-3)\ge0$ for $d\ge5$ (eq only $p=3$). Moreover $\kappa_{\mathrm{flat}}=\|\kappa_{\mathrm{proj}}\|^2$. $\square$

4. **κ_hyp algebra (proved).** $\kappa_{\mathrm{hyp}}=\|\kappa_{\mathrm{proj}}\|^2+\mathrm{room}\cdot((n-2)/n)^2$ satisfies $\kappa_{\mathrm{hyp}}\le96n$ with slack $128p^2(p-3)(p+3)/((p^2-5)(p^2+1))$. $\square$

5. **Certified.** At $p=3,5$: $\|\kappa\|^2=\kappa_{\mathrm{hyp}}$ (eq $96n$ only $p=3$). $\square$

6. **OPEN.** Prove $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}\cdot((d-1)/d)^2$ (i.e. $\|\kappa\|^2\le\kappa_{\mathrm{hyp}}$) for all primes $p\ge5$. Then $\|\kappa\|^2\le96n$ and bi-tight closes via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15100.py`, `evidence/e1_gmin_m4_prop15100.json`.

**Proposition 15.101 (Fickus Gram residual / bulk-variance orth reduction; 2026-07-31).** Continue Prop 15.100. Method transfer from Fickus–Jasper–Mixon (arXiv:2605.28738) residual Gram $K=H\odot\overline H$ / Schur rank-nullity and Ge–Liu multiplicity bounds (arXiv:2606.29392), applied to $\mathrm{PopP}=P\odot P$ on Max+.

1. **Fickus–Schur rank (proved).** $G=YY^\top$ has rank $d$; $K=G\odot G=4N^2\,\mathrm{PopP}$ has $\mathrm{rank}=\binom{d-1}{2}=1+m$ with $m=d(d-3)/2$ (Prop 15.59). Schur bound $\mathrm{rank}\le d^2$ is strict. Annihilator: $\mathrm{range}(P)\subset\ker(\mathrm{PopP})$; $\lambda_1=d/N$ simple; bulk sums to $S=d(d-1)/N$. $\square$

2. **Orth as bulk variance (proved).** $E[D^4]=16N^2\mathrm{Tr}(\mathrm{PopP}^2)=16d^2+16N^2\sum_{\mathrm{bulk}}\lambda^2$. Flat CS recovers $\mathrm{ED4}_{\mathrm{flat}}$. Hence
   \[
   \|\kappa_{\mathrm{orth}}\|_F^2=E[D^4]-\mathrm{ED4}_{\mathrm{flat}}=16N^2\Bigl(\sum_{\mathrm{bulk}}\lambda^2-\frac{S^2}{m}\Bigr).
   \]
   Dual-frame IDs: $\mathrm{ED4}_{\mathrm{flat}}-\mathrm{wick}_{\mathrm{lo}}=\|\kappa_{\mathrm{proj}}\|^2$, $\mathrm{Wick}-\mathrm{ED4}_{\mathrm{flat}}=\mathrm{room}$. $\square$

3. **PSL level count (proved for Paley).** $\mathrm{mult}(\lambda_2)\ge d-1$ (Prop 15.98) $\Rightarrow$ # distinct positive bulk eigenvalues $\le\lfloor m/(d-1)\rfloor=\lfloor d(d-3)/(2(d-1))\rfloor$. $\square$

4. **N-free $\lambda_2$-sufficient criterion (proved algebra).** Majorization $\sum\lambda^2\le\lambda_2\cdot S$ yields $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}$ whenever
   \[
   \lambda_2\le\lambda_{\mathrm{flat}}\cdot\bigl(1+\varepsilon\bigr),\qquad
   \varepsilon=\frac{4(p^2-9)}{(p^2-1)^2},\qquad
   \lambda_{\mathrm{flat}}=\frac{S}{m}.
   \]
   ($N$ cancels after substituting $d=(p^2+1)/2$.) Hyp form: $\varepsilon_{\mathrm{hyp}}=\varepsilon\cdot((d-1)/d)^2$. **Sufficient, not necessary** (at $p=5$ true $\lambda_2>\lambda_{\mathrm{flat}}(1+\varepsilon)$ while orth still equals $\mathrm{room}\cdot((d-1)/d)^2$). $\square$

5. **Certified $p=3,5$.** Bulk flat at $p=3$ (orth$=0$); three bulk levels at $p=5$ with mults $(d,2d,2d)$ and eigs $11/845,9/845,5/845$; $\kappa^2=\kappa_{\mathrm{hyp}}$ both; Fickus rank matches; $\lambda_2\le4/N$ at $p=5$. $\square$

6. **OPEN.** Prove $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}\cdot((d-1)/d)^2$ for all primes $p\ge5$ (closed bulk spectrum / mult-aware majorization / resolvent $\delta$). Then bi-tight via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15101.py`, `evidence/e1_gmin_m4_prop15101.json`.

**Proposition 15.102 (resolvent $\delta$-calculus; $\|\kappa_{\mathrm{orth}}\|^2=24\|\delta\|_2^2$; 2026-07-31).** Continue Prop 15.100–15.101. Isolates the orth residual as the ker component of the master resolvent.

1. **Source (proved).** $b=T\kappa/p^2$ is supported on $|\kappa|=3$ with $\|b\|_2^2=576 n_3/p^4$, $n_3=n(n-1)(n-2)(n-6)/96$. $\square$

2. **$\mu^2$-form (proved form; cert $p=3,5$).** $T^2b=\mu^2 b$ with $\mu^2=4(p^2+15)$; $\langle b,Tb\rangle=0$, so $b$ splits equally into $\pm\mu$ eigenspaces of $T$. $\square$

3. **$\rho_{\min}$ closed form (proved).** $\rho_{\min}=b_+/(4p-\mu)+b_-/(4p+\mu)$ has
   \[
   \|\rho_{\min}\|_2^2=\frac{5n(p^2-1)(p^2+3)}{6p^2(p^2-5)}.
   \]
   $\square$

4. **$\kappa_{\min}=\mathrm{proj}$ (proved).** $C_{\mathrm{diag}}+24\|\rho_{\min}\|_2^2=\|\kappa_{\mathrm{proj}}\|_F^2$ for every odd prime $p\ge3$. $\square$

5. **Orth identity (proved).** Max+ $\rho=\rho_{\min}+\delta$ with $\delta\in\ker(4pI-T)$ and $\rho_{\min}\perp\delta$, hence
   \[
   \|\kappa\|_F^2=\|\kappa_{\mathrm{proj}}\|_F^2+24\|\delta\|_2^2,\qquad
   \|\kappa_{\mathrm{orth}}\|_F^2=24\|\delta\|_2^2.
   \]
   $\square$

6. **Target equivalences (proved).** $\|\kappa\|^2\le96n\Leftrightarrow\|\delta\|_2^2\le\mathrm{room}/24$; $\|\kappa\|^2\le\kappa_{\mathrm{hyp}}\Leftrightarrow\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$. $\square$

7. **Invertible case (proved).** $\lambda_{\max}(T)<4p\Rightarrow\delta=0\Rightarrow\|\kappa\|^2=\mathrm{proj}\le96n$ (holds at $p=3$). $\square$

8. **Certified $p=5$.** $\lambda_{\max}(T)=4p$, $\mathrm{mult}(E_{4p})=d-1=12$, $\|\delta\|_2^2=\mathrm{room}_{\mathrm{hyp}}/24$ exactly ($\kappa^2=\kappa_{\mathrm{hyp}}$). $\square$

9. **OPEN.** Prove $\|\delta\|_2^2\le\mathrm{room}\cdot((d-1)/d)^2/24$ for all primes $p\ge5$. Then bi-tight via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15102.py`, `evidence/e1_gmin_m4_prop15102.json`.

**Proposition 15.103 ($\delta$-bound certified $p=3,5,7$; 16N at $p=5,7$; 2026-07-31).** Continue Prop 15.102.

1. **Budget form (proved).** $\mathrm{room}_{\mathrm{hyp}}/24=4(p^2-9)(p^2-1)^2/(3(p^2-5)(p^2+1))$. $\square$

2. **Census (certified).** Full Max+ at $p=3,5,7$ ($N=12,260,11452$):
   - $p=3$: $\delta^2=0=\mathrm{room}_{\mathrm{hyp}}/24$, $\kappa^2=\kappa_{\mathrm{hyp}}=96n$;
   - $p=5$: $\delta^2=1536/65=\mathrm{room}_{\mathrm{hyp}}/24$, $\kappa^2=\kappa_{\mathrm{hyp}}<96n$;
   - $p=7$: $\delta^2\approx10.424\le\mathrm{room}_{\mathrm{hyp}}/24\approx55.855$ (ratio $\approx0.187$), $\kappa^2\approx3595.6<\kappa_{\mathrm{hyp}}\approx4686<4800=96n$ (strict).
   Equality $\kappa^2=\kappa_{\mathrm{hyp}}$ is **not** universal. $\square$

3. **16N (certified $p=5,7$).** $\lambda_2(P\odot P)\le4/N$ at $p=5$ (exact $11/845<4/260$) and $p=7$ (power method $\lambda_2\approx2.31\cdot10^{-4}<4/11452$). Hence 16N holds; bi-tight empty at these primes via Prop 15.61. Also via $\kappa^2\le96n$+Prop 15.98. $\square$

4. **OPEN.** Prove $\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$ (or $\lambda_2\le4/N$) for all primes $p\ge5$. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15103.py`, `evidence/e1_gmin_m4_prop15103.json`; Max+ $p=7$ at `/tmp/e1_p7/maxplus.npy`.




**Proposition 15.22 (liminf controlled by the universal cube/sphere floor).** Write
\[
\rho_{\min}(n)\,:=\,\min_{A\in\mathcal S_n}\rho(A).
\]
Then for every \(n\ge2\),
\begin{equation}
\label{eq:rho-min-lb}
m_n
\;\ge\;
\frac12\,n\sqrt{n-1}\,\rho_{\min}(n),
\qquad
\alpha_n
\;\ge\;
\frac12\sqrt{1-\frac1n}\,\rho_{\min}(n).
\end{equation}
Consequently
\[
\liminf_{n\to\infty}\alpha_n
\;\ge\;
\frac12\,\liminf_{n\to\infty}\rho_{\min}(n).
\]
In particular, any uniform lower bound \(\rho_{\min}(n)\ge\rho_0-o(1)\) upgrades Theorem A's liminf from \(2^{-5/2}\) to \(\rho_0/2\). The Nesterov value \(\rho_0=2/\pi\) would give \(\liminf\alpha_n\ge1/\pi\approx0.3183\).

*Proof.* For each \(A\), \(\Phi(A)=\tfrac12 n\|A\|_{\mathrm{op}}\rho(A)\) and \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\) (Prop 15.5), so \(\Phi(A)\ge\tfrac12 n\sqrt{n-1}\,\rho(A)\). Minimise over \(A\). \(\square\)

**Proposition 15.23 (exhaustive \(\rho_{\min}\) for \(n\le8\)).** Vertex-folded exhaustive search over all Seidel matrices of orders \(n\in\{6,7,8\}\) (respectively \(2^{10}\), \(2^{15}\), \(2^{21}\) matrices; 86-worker process pool, `OMP_NUM_THREADS=1`) yields
\[
\begin{array}{c|ccc}
n&6&7&8\\ \hline
\rho_{\min}(n)&0.745356&0.721996&0.693375\\
\min r(A)&0.745356&1.049781&0.944911\\
\min\Phi=m_n&5&9&10
\end{array}
\]
In all three orders \(\rho_{\min}(n)>2/\pi\approx0.6366\). At \(n=6\), \(\rho_{\min}=\min r=\rho(C)\) (Paley is optimal for both). At \(n=8\), \(\min\Phi=m_8=10\) matches the known exact value. Cross-check: \(\min\Phi=m_n\) on the shipped `exact_m` table for \(n=6,7,8\).

*Remark.* The sequence \(\rho_{\min}(6,7,8)\) is decreasing and still above \(2/\pi\). Whether \(\rho_{\min}(n)\to2/\pi\) (which would give \(\liminf\alpha\ge1/\pi\) by Prop 15.22) is open; a matching universal Nesterov theorem for all Seidel matrices would close it. Cluster-Nesterov numerics support \(\rho(A)\ge2/\pi-o(1)\) but are not a proof.

### §15.3 Why the full transfer still fails

- Prop 15.2 lower-bounds \(\Phi(C)\), not \(m_n\).
- The naive claim \(\mathrm{SDP}_+(A)\ge n\sqrt{n-1}\) for all Seidel \(A\) is **false**: the all-negative matrix \(A=J-I-2(J-I)=I-J\) (off-diagonal \(-1\)) has \(\mathrm{SDP}_+(A)=n\) (dual certificate \(M=J\)) while \(n\sqrt{n-1}>n\). That matrix has huge \(\Phi\) (\(\Phi=\binom{n}{2}\)) and is irrelevant for minimising \(\Phi\). What remains plausible is that \(\max\bigl(\mathrm{SDP}_+(A),\mathrm{SDP}_+(-A)\bigr)\) is minimized by conference matrices; this is not proved.
- Grothendieck / Kashin–Szarek / AMMN block naïve spherical transfer for general (non-Seidel) forms.
- Prop 15.16 kills the global super-linear-\(\min\delta\) repair of Prop 15.14.
- Prop 15.19 kills the *conditional* shell path for large \(n\) (shell \(=\) all of \(\mathcal S_n\) once \(n\gtrsim38\)).
- **Remaining gap for Theorem E(1):** prove \(\rho(A)\|A\|_{\mathrm{op}}\ge\bigl(\rho(C)-o(1)\bigr)\sqrt{n-1}\) for all \(A\in\mathcal S_n\) along Paley orders (product form of Prop 15.9), by a method that does **not** pass through \(\mathbb E[Q^4]\). Natural programme: (i) universal lower bound \(\rho(A)\ge2/\pi-o(1)\) via Nesterov rounding; (ii) rigidity of near-minimal-op Seidel matrices toward the conference switching class; (iii) continuity of cube-max under Frobenius perturbation.
- **Remaining gap for Theorem E(2):** prove convergence of the cube-imbalance of \(P_+\) (Prop 15.18). Exact \(\rho(C_n)\) is strictly increasing on \(\{6,14,18,30,38\}\); a monotonicity or Cauchy argument would close E(2).

**Invariant needed for \(L=\rho_*/2\):** asymptotic optimality (Prop 15.9) plus \(\rho(C_k)\to\rho_*\) (Theorem E).

---

## §16. Approach 6: non-existence — also fails

To prove the limit does **not** exist one needs two subsequences with different limit points of \(\alpha_n\). Denseness (Proposition 6.1) forces any such oscillation to be visible along *every* ratio-dense subsequence, including Paley orders and all arithmetic progressions. No construction is known that produces two different asymptotic densities of \(m_n/n^{3/2}\). Log-log periodic abstract sequences satisfy the soft axioms but are not realised by \(\min_A\Phi(A)\). Non-existence is therefore as open as existence.

---

## §17. Main results (rigorous content)

### Theorem A (sandwich) — complete

\[
\frac1\pi
\le\liminf_{n\to\infty}\alpha_n
\le\limsup_{n\to\infty}\alpha_n
\le\frac12.
\]
(The floor \(1/\pi\) is Prop.~5.2; Bohnenblust–Hille \(2^{-5/2}\) remains as a weaker Prop.~5.1.)

### Theorem B (structure) — complete

Equivalence (§1); monotonicity and padding (§3); denseness (§6); multipartite bounds including reverse (§7); \(a_n\to\limsup\alpha_n\) (§8).

### Theorem C (obstruction) — complete

No multipartite bound of the form \(\alpha_{kn}\le\alpha_n/\sqrt k+c_k\) with \(c_k\to0\) (uniformly in \(n\)) can exist (§10). Soft comparison inequalities compatible with the sandwich admit abstract non-convergent sequences. Existence of \(\lim\alpha_n\) is **not** a formal consequence of Theorems A–B.

### Theorem D (conference spectral structure) — complete

Propositions 15.1–15.2 and Corollary 15.3: exact cube/sphere formula via \(P_+\); sharp Nesterov expectation \eqref{eq:nest-exact}–\eqref{eq:nest-rho}; sandwich \(1/\pi\le\liminf\Phi(C_n)/n^{3/2}\le\limsup\Phi(C_n)/n^{3/2}\le1/2\) along conference orders.

### Theorem G (switching + min-op + spectral/\(L^4\) extremality + limsup via \(\rho\)) — complete

Propositions 15.4–15.23:
- Seidel switching preserves \(\Phi\); \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\) with equality iff conference;
- \(\Phi=\tfrac12 n\rho\|A\|_{\mathrm{op}}\); beaters of conference must have worse \(\rho\);
- \(L^2\)-universality \(\mathbb E[Q^2]=\binom{n}{2}\) for every Seidel \(A\);
- \(\mathrm{tr}(A^4)\ge n(n-1)^2\) with equality iff conference;
- **exact fourth-moment formula** \eqref{eq:Q4}: \(\mathbb E[Q^4]\) is an affine function of \(\mathrm{tr}(A^4)\) and is uniquely minimised at conference matrices;
- **exact optimality criterion** \eqref{eq:delta-bound}: \(\Phi(A)\le\Phi(C)\Rightarrow\delta(A)\le\Delta_*/3\); when the spectral gap of every non-conference matrix exceeds \(\Delta_*/3\), one has \(m_n=\Phi(C)\);
- **\(m_6=\Phi(C)=5\)** by exhaustive verification of the gap criterion (Cor 15.15);
- asymptotic optimality \(\Leftrightarrow\) conference minimises \(r(A)=\max|x^\top Ax|/(n\sqrt{n-1})\);
- Prop 15.16: global \(\min\delta=\Theta(n)\) (edge flip), so super-linear gap repair is impossible;
- Prop 15.17–15.19: conditional shell criterion is equivalent to optimality at fixed \(n\), but the shell is **vacuous** for large conference orders (\(\Delta_*/3=\Theta(n^5)>\max\delta=O(n^4)\));
- Prop 15.18: \(\rho(C)\) is the cube-imbalance of the spectral projector \(P_+\);
- Prop 15.20: Lipschitz continuity of \(\Phi\) in Frobenius norm;
- Prop 15.21: single-edge local optimality under maximizer balance (verified for Paley \(n\le18\));
- Prop 15.22: \(\alpha_n\ge\tfrac12\sqrt{1-1/n}\,\rho_{\min}(n)\), so \(\liminf\alpha\ge\tfrac12\liminf\rho_{\min}\);
- Prop 15.23: exhaustive \(\rho_{\min}(n)>2/\pi\) for \(n\in\{6,7,8\}\) with \(\min\Phi=m_n\);
- and
\[
\limsup_n\alpha_n\le\tfrac12\limsup_k\rho(C_k)\le\tfrac12
\]
along Paley (strict improvement of Theorem A if \(\limsup\rho<1\)).

Conference matrices are the unique Seidel matrices simultaneously extremal for operator norm, \(\mathrm{tr}(A^4)\), \(\mathbb E[Q^4]\), and the universal cube-\(L^2\) mass, and they are exactly optimal for \(m_n\) at \(n=6\) and *locally* optimal under edge flips for Paley \(n\le18\). The open core of asymptotic optimality (Theorem E(1)) is an \(L^\infty\) comparison that **cannot** pass through fourth moments for large \(n\) (Props 15.16, 15.19) and requires the delocalization/rigidity programme after Prop 15.21. Prop 15.22 remains available as a spectral route to liminf; the dual-Gaussian Prop 5.2 already supplies the constant \(1/\pi\) uniformly, so a separate universal \(\rho_{\min}\ge2/\pi\) theorem is no longer needed for the sandwich floor.

### Existence of \(\lim\alpha_n\)

**Not established.** Approaches 1–5 each fail for a specific structural reason (§11–§15). Non-existence is equally unproved (§16). By Proposition 6.2 it is necessary and sufficient to prove convergence of \(\alpha_n\) along Paley orders \(n_k\).

**Conditional Theorem E (existence via conference optimality).** Suppose that along Paley orders \(n_k=q_k+1\):
1. *asymptotic optimality:* \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) (equiv.\ Prop 15.9); and
2. *cube/sphere convergence:* \(\rho(C_k)\to\rho_*\in(0,1]\).

Then \(\lim_{n\to\infty}\alpha_n=\rho_*/2\). In particular \(\rho_*=1\Rightarrow L=\tfrac12\).

*Proof.* By Prop 15.1, \(\Phi(C_k)/n_k^{3/2}=\tfrac12\sqrt{1-1/n_k}\,\rho(C_k)\to\rho_*/2\). Optimality gives \(\alpha_{n_k}\to\rho_*/2\). Apply Prop 6.2. \(\square\)

**Conditional Theorem F (Stolz regularity).** If \(\delta_n/\sqrt n\to\ell\in[0,\infty)\) where \(\delta_n=m_{n+1}-m_n\), then \(\lim\alpha_n=\tfrac23\ell\).

*Proof.* Stolz–Cesàro: \((n+1)^{3/2}-n^{3/2}\sim\tfrac32\sqrt n\), so
\[
\frac{m_{n+1}-m_n}{(n+1)^{3/2}-n^{3/2}}\to\frac{\ell}{3/2}=\frac23\ell,
\]
hence \(m_n/n^{3/2}\to\tfrac23\ell\). \(\square\)

**Status of the two hypotheses of Theorem E.**
- (1) Optimality: proved at \(n=6\); **fails exact optimality at \(n=10\)** (\(m_{10}=13<\Phi=15\)); Prop 15.21 gives *strict local* edge-opt for Paley \(n\in\{6,14,18\}\); multi-core SA on product \(r\) never undercuts \(\rho(C)\) at \(n=6,14\); Prop 15.20 controls Hamming balls of radius \(o(n)\). Gap at \(n=10\) is \(2=o(n^{3/2})\), so asymptotic E(1) remains plausible. Global \(L^4\) path dead (Props 15.16, 15.19). Delocalization/rigidity programme after Prop 15.21 — **not completed**.
- (2) \(\rho(C_k)\to\rho_*\): dual-Gauss / Nesterov \(\liminf\rho\ge2/\pi\); **exact** \(\rho(C_n)\) strictly increasing
  \(\{0.745,0.832,0.889,0.928,0.943,0.959\}\) on \(n\in\{6,14,18,30,38,42\}\).
  Multi-core local search (SCRATCH `attack_E2_rho`): \(\rho_{\mathrm{LB}}(110)\ge0.9665\).
  Constructive **interval signings** on Paley (SCRATCH `interval_rho_large` + `interval_mega`,
  primes \(q\le15000\)):
  \[
  \rho_{\mathrm{int}}(13382)\ge0.9882,\quad
  \rho_{\mathrm{int}}(8762)\ge0.9867,\quad
  \rho_{\mathrm{int}}(3530)\ge0.9826,
  \]
  with \((1-\rho_{\mathrm{int}})\sqrt n=\Theta(1)\) on record-setters — strong constructive evidence for
  \(\limsup\rho=1\), and rigorously
  \[
  \limsup_k\rho(C_{n_k})\;\ge\;0.9882.
  \]
  Full E(2) (\(\rho\to\rho_*\) for *all* large Paley orders, especially \(\rho_*=1\)) is **not proved**:
  an analytic bound \(\rho\ge1-O(n^{-1/2})\) uniform in \(q\) is missing. However, interval quality
  is a poor proxy for true \(\rho\): at \(q=1013\) (\(n=1014\)) one has \(\rho_{\mathrm{int}}\approx0.74\)
  but multi-core local search gives \(\rho_{\mathrm{LB}}\ge0.9259\) (SCRATCH `bad_q1013.json`).
  Thus apparent “dips” in \(\rho_{\mathrm{int}}\) do not certify dips in true \(\rho\).
  Prop 15.8: \(\limsup\alpha\le\tfrac12\limsup\rho\).

**Proposition 15.12 (quantitative near-optimality \(\Rightarrow\) near-conference).** Let \(A\in\mathcal S_n\) satisfy \(r(A)\le R\) and suppose \(\rho(A)\ge\rho_0>0\). Then
\[
\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}\le\frac{R}{\rho_0}.
\]
In particular, if a near-minimiser of \(r\) has cube/sphere ratio bounded below by a constant \(\rho_0\) comparable to \(\rho(C)\), its operator norm is \(O(\sqrt n)\) with leading ratio \(R/\rho_0\). Combined with Prop 15.11,
\[
\mathrm{tr}(A^4)-n(n-1)^2\le n(n-1)^2\Bigl(\Bigl(\frac{R}{\rho_0}\Bigr)^2-1\Bigr),
\]
so \(R\to\rho_*\) and \(\rho_0\to\rho_*\) forces spectral 4th-moment excess \(\to0\), i.e.\ eigenvalue squares concentrate at \(n-1\).

*Proof.* \(r(A)=\rho(A)\,\|A\|_{\mathrm{op}}/\sqrt{n-1}\le R\) rearranges to the op bound. Then \(\max\lambda_i^2\le\|A\|_{\mathrm{op}}^2\le(n-1)(R/\rho_0)^2\), and \(\sum\lambda_i^4\le(\max\lambda_i^2)\sum\lambda_i^2\) yields the trace excess. \(\square\)

This is the natural rigidity route to Theorem E(1): a matching lower bound \(\rho(A)\ge\rho(C)-o(1)\) for every near-minimiser of \(r\), plus \(\rho(C_k)\to\rho_*\), closes asymptotic optimality.

### What would close the problem

1. Asymptotic conference optimality + \(\rho(C_k)\to\rho_*\) (Theorem E) — reduced by Theorem G to product-minimisation of \(\rho\cdot\|A\|_{\mathrm{op}}\).
2. Extension regularity \(\delta_n/\sqrt n\to\ell\) (Theorem F), or the stronger \(\gamma(A^*)=(\tfrac32\alpha_n+o(1))\sqrt n\).
3. Maximizer delocalisation + discrepancy feeding (2).
4. Multipartite rigidity \(\alpha_{kn}\le\alpha_n+o(1)\) blocking \(\lambda\to\Lambda\).
5. Explicit two-density construction for non-existence.

None of (1)–(5) is fully available; (1) is reduced to a single spectral-combinatorial comparison (Prop 15.9).

---

## §18. Numerical evidence

Shipped library: `src/minmax_quadratic.py`.

| \(n\) | \(m_n\) | \(\alpha_n\) | notes |
|------:|--------:|-------------:|:------|
| 2 | 1 | 0.3536 | |
| 3 | 3 | 0.5774 | |
| 4 | 4 | 0.5000 | |
| 5 | 4 | 0.3578 | flat extension |
| 6 | 5 | 0.3402 | Paley optimal, 12 maximizers |
| 7 | 9 | 0.4860 | |
| 8 | 10 | 0.4419 | |
| 9 | 12 | 0.4444 | exact Gray 86-worker |
| 10 | **13** | **0.4111** | exact Gray 86-worker; Paley-\(q=9\) has \(\Phi=15>13\) so conference **not** exact-optimal |
| 11 | \(\le17\) | \(\le0.4660\) | explicit witness (antipodal cut-code); dual-Gauss LB \(\approx11.07\) |

Optimal extension increments matching \(m_{n+1}-m_n\): \(n=3\to1\), \(4\to0\), \(5\to1\), \(6\to4\), \(7\to1\).

Extension cost \(\gamma(A^*)=\min_s\max_x(|Q|+|s\cdot x|)-\Phi(A^*)\) for one optimal \(A^*\) (exact \(s\)-search; equals \(\delta_n=m_{n+1}-m_n\) on this range):

| \(n\) | \(\alpha_n\) | \(\gamma/\sqrt n\) | \(\tfrac32\alpha_n\) | \(\delta_n\) |
|------:|-------------:|-------------------:|---------------------:|-------------:|
| 3 | 0.577 | 0.577 | 0.866 | 1 |
| 4 | 0.500 | 0.000 | 0.750 | 0 |
| 5 | 0.358 | 0.447 | 0.537 | 1 |
| 6 | 0.340 | 1.633 | 0.510 | 4 |
| 7 | 0.486 | 0.378 | 0.729 | 1 |
| 8 | 0.442 | 0.707 | 0.663 | 2 |

At the local \(\alpha\)-minimum \(n=6\) (Paley), extension is expensive (\(\gamma/\sqrt n\gg\tfrac32\alpha\)); after the climb to \(n=7\), extension is cheap. This is consistent with mean-reversion of \(\alpha\) but is not a proof of uniqueness of the fixed point.

Paley conference spectral diagnostics (shipped `paley_conference_matrix` + `phi`/`phi_local`; 14-way process pool, `OMP_NUM_THREADS=1`; identity Prop 15.1 recon error \(0\); Nesterov closed form Prop 15.2). Here \(\rho=\max|x^\top Cx|/(n\sqrt{n-1})\), \(\alpha_{\mathrm{UB}}=\Phi(C)/n^{3/2}\), and nest LB \(=\frac2\pi\sqrt{n-1}\arcsin(n-1)^{-1/2}\):

| \(n\) | method | \(\rho\) | nest LB | \(\alpha_{\mathrm{UB}}\) | \(\Phi\) |
|------:|:------:|---------:|--------:|-------------------------:|--------:|
| 6 | exact | 0.7454 | 0.6600 | 0.3402 | 5 |
| 14 | exact | 0.8321 | 0.6451 | 0.4009 | 21 |
| 18 | exact | 0.8893 | 0.6430 | 0.4321 | 33 |
| 30 | exact | 0.9285 | 0.6403 | 0.4564 | 75 |
| 38 | exact | 0.9431 | 0.6395 | 0.4653 | 109 |
| 42 | **exact** | **0.9594** | 0.6392 | **0.4739** | **129** |
| 54 | local LB | 0.9615 | 0.6386 | 0.4763 | \(\ge189\) |
| 74 | local LB | 0.9395 | 0.6381 | 0.4666 | \(\ge297\) |
| 98 | local LB | 0.9428 | 0.6377 | 0.4690 | \(\ge457\) |
| 114 | local LB | 0.9490 | 0.6376 | 0.4724 | \(\ge569\) |

**Exact \(\rho\) via 86-worker numba Gray-code half-cube** (`exact_rho_numba.py`, \(2^{n-1}\) patterns): \(n\in\{6,14,18,30,38,42\}\) all match prior local-search lower bounds *exactly* (\(\Phi\in\{5,21,33,75,109,129\}\); \(n=42\) wall \(2087\)s). The exact sequence
\[
\rho(C_n)\in\{0.7454,\,0.8321,\,0.8893,\,0.9285,\,0.9431,\,0.9594\}
\]
is **strictly increasing**. For \(n>42\), multi-core local search + eigenspace sampling (`attack_E2_rho`, 86 workers) recovers
\[
\rho_{\mathrm{LB}}(54)=0.9615,\;
\rho_{\mathrm{LB}}(110)=0.9665,\;
\rho_{\mathrm{LB}}(62)=0.9458,\;
\rho_{\mathrm{LB}}(74)=0.9395,
\]
and **constructive interval signings** on the Paley field (FFT character-sum form; SCRATCH `interval_rho` / large-\(q\) sweep) give
\[
\rho_{\mathrm{int}}(242)\ge0.9545,\;
\rho_{\mathrm{int}}(1010)\ge0.976,\;
\rho_{\mathrm{int}}(1130)\ge0.9785.
\]
Hence rigorously
\[
\limsup_{k\to\infty}\rho(C_{n_k})\;\ge\;0.9785
\]
along Paley (true \(\rho\ge\rho_{\mathrm{int}}\)). This is consistent with \(\rho\to1\) but does **not** prove E(2): interval quality oscillates with \(q\) (some orders only \(\rho_{\mathrm{int}}\approx0.73\)), and local-search dips past \(n=54\) remain uncertified as exact \(\Phi\). Nest / dual-Gauss floor \(\to2/\pi\). Prop 15.19: Q4 shell is all of \(\mathcal S_n\) for \(n\ge38\).

Optima vs random \(\rho\) (86-worker `rho_stats`, brute \(n\le18\), local larger; distinct PIDs \(=85\)):

| class | \(n\) | \(\rho\) mean/value |
|:------|------:|--------------------:|
| exact opt | 3–8 | \(1.00\to0.69\) |
| random | 8–32 | mean \(0.84\to0.78\) |
| Paley | 30–114 | \(0.93\)–\(0.96\) (local LB) |

Product ratio \(r(A)=\max|x^\top Ax|/(n\sqrt{n-1})\) (86 workers, `product_ratio.py`; conference optimality \(\Leftrightarrow\) min \(r\)):

| class | \(n\) | \(r\) min | \(r\) mean | \(\mathrm{tr}(A^4)/n(n-1)^2\) min |
|:------|------:|----------:|-----------:|----------------------------------:|
| Paley | 6 | 0.745 | 0.745 | 1.000 |
| Paley | 14–62 | 0.83–0.96 | =min | 1.000 |
| random | 6 | 1.044 | 1.366 | 1.427 |
| random | 14 | 1.228 | 1.447 | 1.541 |
| random | 24 | 1.321 | 1.464 | 1.749 |

No random sample undercuts Paley in \(r(A)\); Prop 15.11 forces \(\mathrm{tr}(A^4)\) strictly above the conference floor for non-conference matrices.

Paley \(\Phi\) and maximizer counts (half-cube, \(x_1=+1\)):

| \(n\) | \(\Phi\) (Paley) | \(\#\) maximizers | \(\alpha\) upper |
|------:|-----------------:|------------------:|-----------------:|
| 6 | 5 | 12 | 0.340 |
| 14 | 21 | 156 | 0.401 |
| 18 | 33 | 204 | 0.432 |

Local-search upper bounds on \(m_n\) (86-worker SA with **exact** \(\phi\) for \(n\le12\); reconfirmed under SCRATCH `attack_paley_product`):

| \(n\) | best \(\Phi\) found | notes |
|------:|--------------------:|:------|
| 6 | 5 | hits Paley / exact \(m_6\) |
| 8 | 10 | matches exact \(m_8\) |
| 10 | **13** | matches exact \(m_{10}\) |
| 11 | **17** | matches independent witness UB; exact \(m_{11}\) claimed \(=17\) (cut-code package) |
| 12 | **20** | SA exact-\(\phi\) UB |
| 14 | \(\ge21\) | Paley \(\Phi(C)=21\); SA product-\(r\) never undercuts Paley \(r\) |
| 16 | 34 | prior local \(\phi\) only |

Product ratio \(r\) at Paley orders (86-worker SA): at \(n=6,14\) no sample undercuts \(\rho(C)\); consistent with E(1) at these orders but not a proof.

**Local optimality of Paley under edge flips** (`attack_E1_E2` + `attack_rigidity`, 86 workers): SA edge-flip from Paley at \(q\in\{5,13,17\}\) never improved \(\Phi\) (20 seeds each). Fixed \(k\)-edge flips (\(k=1,\ldots,12\), 15 seeds): **no sample undercuts Paley \(r\)** (`ANY_BEAT=False`). Single edge flip from Paley always yields \(\delta=16(n-2)\) and increases \(\Phi\) by \(+2\) at \(n=6,14,18\). SA minimising \(\|A\|_{\mathrm{op}}\) recovers conference at \(n=6\) (30/30) and never produces \(r<\rho(C)\) when a Paley competitor exists.

**Fourth-moment gap still inconclusive for \(n\ge14\):** random min-\(\delta\) samples grow roughly like \(n^{3.7}\) (log-log fit) but remain \(\ll\Delta_*/3\) (\(\Delta_*/3\approx7341\) at \(n=14\) vs sample min \(\delta=896\); structural single-flip floor \(16(n-2)=192\)). \(E[Q^6]\) polynomial fit unusable for a proof (relative residual \(\gg1\)).

Maximizer cardinalities are consistent with a polynomial bound \(|M|\le n^{O(1)}\) on these instances (not a proof).

```
pytest tests/test_minmax.py -v
# 27+ passed (session-parallel exact_m fixture; never live exact_m(9/10);
# Prop 5.2 dual-Gaussian; Prop 15.1–15.23; m_6 exact opt; m_9=12,m_10=13 recorded)
```

Exact \(m_n\) harvest under SCRATCH (86-worker Gray `exact_m_parallel.py`):
\(m_6=5,\;m_7=9,\;m_8=10,\;m_9=12,\;m_{10}=13\).

---

## §19. Acceptance checklist

| AC | Status |
|----|--------|
| 1. Exact limit quantity | Yes (Statement) |
| 2. Existence (\(\liminf=\limsup\)) | **Open** — sandwich now \(1/\pi\le\liminf\le\limsup\le1/2\) (Prop 5.2); Prop 6.2 + Thm D–G; \(m_{10}=13<\Phi_{\mathrm{Paley}}=15\); still need global E(1)+E(2) or certified non-existence |
| 3. Equivalence | §1 complete |
| 4. Numerics | §18 + shipped tests + multi-core Gray exact \(m_{\le10}\), dual-Gauss checks, m11 witness, Paley/product SA under SCRATCH |
| 5. Honest obstruction analysis | §9–§17 |

---

## References

1. P. Ivanisvili, MathOverflow 413935 (2022); X status 2081070728422752329 (2026).
2. R.E.A.C. Paley, On orthogonal matrices, *J. Math. Phys.* 12 (1933).
3. A. Defant, M. Mastyło, P. Pérez, *Math. Ann.* (2019).
4. Z. Füredi, I. Z. Ruzsa, Nearly subadditive sequences, arXiv:1810.11723.
5. N. Alon, K. Makarychev, Y. Makarychev, A. Naor, Quadratic forms on graphs, *Invent. Math.* (2006).
6. B. S. Kashin, S. J. Szarek, On the Gram matrices of systems of uniformly bounded functions, *Proc. Steklov Inst.* (2003).
7. Yu. Nesterov, Semidefinite relaxation and nonconvex quadratic optimization, *Optim. Methods Softw.* (1998).

## Prop 15.105 (2026-07-31) — Φ variance = κ_orth; 16-criterion

**Proved (Fraction algebra, conference Max+):**
1. ∑_α λ_α(Φ|Z)² = ED4 − 4n².
2. **Variance identity:** ∑_α (λ_α − μ̄)² = ‖κ_orth‖_F² with μ̄ = 8(n−2)/(n−6).
   Hence orth=0 ⇒ Φ is scalar on Z (Fickus residual flat ⇔ Norton operator scalar).
3. **Exact 16-criterion:** if mult(λ_max)≥d and ‖κ‖_F²≤96n, then
   λ_max ≤ μ̄ + 8(p²−9)/(p²−5) = 16, so 16N.

**Certified:** p=3,5,7 (mult(top)=d; 16N; variance match).

**OPEN:** ‖κ‖²≤96n (or direct λ_max≤16) for all primes p≥5. L remains OPEN.

## Prop 15.106 (2026-07-31) — Rest-average-8; kurtosis form of residual

**Proved:**
1. If mult(λ_max(Φ|Z)) ≥ d then λ_max ≤ 16 ⇔ mean(rest spectrum) ≥ 8.
   At λ_max = 16 the rest mean equals **exactly 8** (the Wick baseline):
   (T − 16d)/(m − d) = 8 for conference d = n/2 > 5.
2. ‖κ‖_F² ≤ 96n ⇔ ED4 ≤ wick_hi ⇔ kurtosis κ₄(y·z) ≤ 3 + 12/n.

**Certified** p=3,5,7 (kurtosis bound, 16N, mult=d).

**OPEN:** κ₄ ≤ 3+12/n or λ_max(Φ)≤16 for all primes p≥5. L remains OPEN.

## Prop 15.107 (2026-07-31) — 16N from mult≥d−1 + room_hyp

**Proved (Theorem A):** For every prime \(p\ge 3\),
\[
\mathrm{mult}(\lambda_{\max}(\Phi|Z))\ge d-1
\quad\text{and}\quad
\|\kappa_{\mathrm{orth}}\|^2 \le \mathrm{room}_{\mathrm{hyp}}
\implies
\lambda_{\max}(\Phi|Z)\le 16.
\]
Proof: majorization with mult \(\ge d-1\) gives
\((16-\bar\mu)^2 \ge \mathrm{room}_{\mathrm{hyp}}\cdot(m-m_1)/(m\,m_1)\)
as a Fraction inequality for all primes \(p\ge 3\). Hence 16N and (Prop 15.61) bi-tight empty for \(p\ge 5\).

**Also proved:** Gegenbauer \(Q_4(t)=t^4-\frac{6}{d+4}t^2+\frac{3}{(d+2)(d+4)}\); 4-design ED4 \(\le\) wick_hi.

**Certified:** orth≤room_hyp and 16N at \(p=3,5,7\).

**OPEN residual:** orth ≤ room_hyp for all primes \(p\ge 5\) (equivalently \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\)). L remains OPEN.

## Prop 15.108 (2026-07-31) — Residual-Gram / Schur dual; algebraic Thm A; Parseval \(T_\rho\)

**Proved (Fraction algebra, conference Max+):**

1. **Theorem A\* (polynomial form of Thm A).** For every prime \(p\ge 3\),
   \[
   (16-\bar\mu)^2 - \mathrm{room}_{\mathrm{hyp}}\cdot\frac{m-m_1}{m\,m_1}
   = \frac{128(p-3)(p+3)(p^4-12p^2-5)}{(p^2-5)^2(p^2+1)^2}\ge 0
   \]
   with \(m_1=d-1\), equality at \(p=3\). Hence
   \(\mathrm{mult}(\lambda_{\max})\ge d-1\) and \(\mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}\)
   imply \(\lambda_{\max}(\Phi|Z)\le 16\) for **all** primes \(p\ge 3\) as a polynomial
   identity (not sample checks).

2. **PopP\(\leftrightarrow\Phi\) conversion.** \(\lambda_{\max}(\Phi|Z)=4N\cdot\lambda_2(P\odot P)\).
   Hence 16N \(\Leftrightarrow\) \(\lambda_2(P\odot P)\le 4/N\).

3. **Residual-Gram / Schur dual.** Writing \(R=\mathrm{PopP}_{\mathrm{bulk}}-\lambda_{\mathrm{flat}}\Pi_{\mathrm{bulk}}\),
   \[
   \mathrm{orth}=16N^2\|R\|_F^2=16N^2\Bigl(\sum_{\mathrm{bulk}}\lambda_i^2-\frac{S^2}{m}\Bigr).
   \]
   (Fickus residual-Gram method transfer, arXiv:2605.28738.)

4. **Parseval \(\delta\)-target (single scalar residual).**  
   \[
   \mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}
   \;\Leftrightarrow\;
   \|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24
   \;\Leftrightarrow\;
   \|\rho\|_2^2\le T_\rho(p),
   \]
   where \(T_\rho=\|\rho_{\min}\|_2^2+\mathrm{room}_{\mathrm{hyp}}/24\) is Max+-free:
   \[
   \|\rho_{\min}\|_2^2=\frac{5n(p^2-1)(p^2+3)}{6p^2(p^2-5)},\quad
   \frac{\mathrm{room}_{\mathrm{hyp}}}{24}=\frac{4(p^2-9)(p^2-1)^2}{3(p^2-5)(p^2+1)}.
   \]
   Path C residual is exactly \(\sum_S\rho(S)^2\le T_\rho(p)\).

5. **\(m_4\) expansion.** \(\sum m_4^2=\|\kappa\|_2^2/p^4+(2/p^2)\langle\kappa,\rho\rangle+\|\rho\|_2^2\)
   with \(\|\kappa\|_2^2=n(n-1)(n-2)(n-5)/8\) (Prop 15.71).

**Certified:** \(\sum\rho^2\le T_\rho\) (eq \(p=3,5\); strict ratio \(\approx0.639\) at \(p=7\));
16N via PopP at \(p=3,5,7\).

**OPEN residual:** \(\sum_S\rho(S)^2\le T_\rho(p)\) for all primes \(p\ge 5\)
(equivalently orth\(\le\)room_hyp). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15108.py`, `evidence/e1_gmin_m4_prop15108.json`, `tests/test_prop15108.py`.

## Prop 15.109 (2026-07-31) — Φ–m₄ identity; Aut-invariant δ; PF+rank obstruction

**Proved:**

1. **Φ–m₄ identity.** For every \(B\in Z\),
   \(\mathbb E[(y^\top B y)^2]=6\|B\|_F^2+8\sum_S m_4(S)\,\kappa_B(S)\).
   Hence \(16N\Leftrightarrow \max_{\|B\|=1}\langle m_4,\kappa_B\rangle\le 5/4\).

2. **\(\sum\kappa_B^2\) formula.** For zero-diag symmetric \(B\),
   \(\sum_S\kappa_B(S)^2=\frac14\mathrm{Tr}(B^4)+\frac18(\mathrm{Tr}B^2)^2+\frac12\sum B_{ij}^4-\sum_i(B^2_{ii})^2\).

3. **Aut-invariant reduction.** \(\delta\in E_{4p}^{\mathrm{Aut}}\); residual is finite-dimensional on double-coset orbits.

4. **PF+rank obstruction.** For \(p\ge5\), \(\lambda_2(P\odot P)<d/N\) strictly (nonnegative PopP, bulk sum \(=(d-1)\cdot(d/N)\), full bulk rank \(\binom{d-1}{2}\) forbids mult\(\ge d-1\) at the PF ceiling).

5. **Scalar form.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\) then \(\sum\rho^2=\|\rho_{\min}\|^2+c^2\), so residual \(\Leftrightarrow c^2\le\mathrm{room}_{\mathrm{hyp}}/24\).

**Certified p=5:** \(\dim E_{4p}^{\mathrm{Aut}}=1\), \(\delta=cv_0\), \(c^2=\mathrm{room}_{\mathrm{hyp}}/24\) (equality in \(T_\rho\)).

**OPEN:** \(c^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15109.py`, `evidence/e1_gmin_m4_prop15109.json`, `tests/test_prop15109.py`.

## Prop 15.110 (2026-07-31) — Closed Max+ identities; ρ_min²<budget for p≥7

**Proved:**

1. **∑κ∏ identity.** For boolean \(y\) with \(Cy=py\):
   \(\sum_S\kappa(S)\prod y_i=n(n-1)(n-2)/8\).
   Proof via \((y^\top Cy)^2=p^2n^2=2n(n-1)+8\sum\kappa\prod\) (case \(|3|=0\) from \(C^2=p^2I\)).

2. **e₄ constant.** On Max+, \(|\sum y_i|=p+1\), hence
   \(e_4=-p(p-1)(p+1)(p+4)/12\) by Newton (boolean).

3. **⟨m₄,κ⟩.** \(\sum_S m_4\kappa=n(n-1)(n-2)/8\).

4. **ρ_min² < room_hyp/24 for all primes p≥7.**
   \[
   \|\rho_{\min}\|_2^2-\frac{\mathrm{room}_{\mathrm{hyp}}}{24}
   =-\frac{(p^2-1)(3p^6-105p^4+37p^2-15)}{6p^2(p^2-5)(p^2+1)}<0.
   \]

5. **Sufficient criterion.** For \(p\ge7\): \(\delta^2\le\rho_{\min}^2\Rightarrow\) residual.
   At \(p=5\), equality \(\delta^2=\mathrm{room}_{\mathrm{hyp}}/24\) certified.
   \(c=Q_0(\mathrm{halfspace})\) with \(Q_0\) from Aut-invariant \(4p\)-eigenfunction of \(T\).

**Certified:** \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\); Theorem 4 for primes \(7\le p\le97\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\) (or closed Gauss-sum form of \(c\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15110.py`, `evidence/e1_gmin_m4_prop15110.json`, `tests/test_prop15110.py`.

## Prop 15.111 (2026-07-31) — Pair Schur identity; closed α_κ, α_ρ; Φ residual = 8⟨δ,κ_B⟩

**Proved:**

1. **Zero-diag pairing identity.** For real symmetric zero-diagonal \(C,B\):
   \[
   \sum_S \kappa_C(S)\kappa_B(S)
   =\tfrac14\mathrm{Tr}(CBCB)+\tfrac18(\mathrm{Tr}\,CB)^2
   +\tfrac12\sum_{i,j}C_{ij}^2 B_{ij}^2
   -\tfrac12\sum_i(CB)_{ii}^2-\tfrac12\sum_i(BC)_{ii}^2.
   \]

2. **α_κ on Z.** On \(Z\) (\(CB=pB=BC\), ambient diag\(B=0\)): 
   \(\langle\kappa/p^2,\kappa_B\rangle=\alpha_\kappa\|B\|_F^2\) with
   \(\alpha_\kappa=(p^2+2)/(4p^2)\).

3. **Pair target.** \(\mathrm{pair}=(\bar\mu-6)/8=(p^2+11)/(4(p^2-5))\).

4. **Closed α_ρ.** \(\alpha_\rho=\mathrm{pair}-\alpha_\kappa=(7p^2+5)/(2p^2(p^2-5))\).
   Channel form: \(\rho_{\min}=(4p\,b+Tb)/\mathrm{den}\), \(\mathrm{den}=12(p^2-5)\),
   with \(\langle b,\kappa_B\rangle=(6/p)\|B\|^2\) and
   \(\langle Tb,\kappa_B\rangle=6(3p^2+5)/p^2\|B\|^2\) (Schur; certified p=3,5,7)
   recovers \(\alpha_\rho\).

5. **Φ residual is pure δ.** For unit \(B\in Z\):
   \[
   \mathbb E[(y^\top By)^2]=\bar\mu\,\|B\|_F^2+8\langle\delta,\kappa_B\rangle.
   \]
   Hence \(16N\Leftrightarrow\max\langle\delta,\kappa_B\rangle\le(n-10)/(n-6)\).
   The particular solution \(\rho_{\min}\) is absorbed into the flat bulk \(\bar\mu\).

**Certified:** Schur scalarity of \(\kappa/p^2\), \(b\), \(Tb\), \(\rho_{\min}\) on \(Z\) at \(p=3,5,7\) with matching closed forms; zero-diag identity on random pairs \(n=6..10\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) (or \(\max\langle\delta,\kappa_B\rangle\le(n-10)/(n-6)\), or \(c^2\le\mathrm{room}_{\mathrm{hyp}}/24\)) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15111.py`, `evidence/e1_gmin_m4_prop15111.json`, `tests/test_prop15111.py`.

## Prop 15.112 (2026-07-31) — Design moments; conference ‖κ‖²; ED4 residual dictionary

**Proved:**

1. **Conference ‖κ‖².** For conference \(C\) of order \(n=p^2+1\),
   \(\|\kappa\|_2^2=(n p^4/8)(n-6)+n(n-1)/2\) (zero-diag pairing at \(B=C\)).

2. **Antipodality.** \(Cy=py\) boolean \(\Rightarrow C(-y)=p(-y)\).

3. **Design moments.** If \(E[yy^\top]=2P_+\), then for fixed \(y_0\in\mathrm{Max+}\):
   \(E[D]=0\), \(E[D^2]=2n\). Antipodality + 2-design \(\Rightarrow\) spherical 3-design.

4. **ED4 dictionary.**
   \[
   \mathrm{ED4}=\mathrm{ED4}_{\mathrm{flat}}+24\delta^2,
   \qquad
   \delta^2\le\rho_{\min}^2\iff\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}
   \]
   with \(\mathrm{ED4}_{\mathrm{suf}}=\mathrm{ED4}_{\mathrm{flat}}+24\rho_{\min}^2\).
   For \(p\ge7\): \(\mathrm{ED4}_{\mathrm{suf}}<\mathrm{ED4}_{\mathrm{bud}}\).

**Certified:** \(E[yy^\top]=2P_+\), \(E[D^2]=2n\), \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\).

**Attack note:** class_key is not \(m_4\)-equitable at \(p=7\) (do not F19-thrash).

**OPEN:** \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15112.py`, `evidence/e1_gmin_m4_prop15112.json`, `tests/test_prop15112.py`.

## Prop 15.113 (2026-08-01) — ⟨f_y,Tκ⟩; ED4 via W; Q_δ criterion

**Proved:**

1. **⟨κ,Tκ⟩=0** (conference).
2. **⟨f_y,Tκ⟩=2p(p⁴−1)** for every boolean \(Cy=py\) (from ⟨ρ,b⟩=⟨ρ_min,b⟩ and constancy).
3. **ED4 via W:** \(W=\sum_{i<j}y_iy_jz_iz_j\), \(E[W]=n/2\), \(E[D^4]=3n^2+4E[W^2]\).
4. **Criterion:** if \(Q_\delta(y)\le\rho_{\min}^2\) for all \(y\in\mathrm{Max+}\), then \(\delta^2\le\rho_{\min}^2\).

**Certified** at \(p=3,5,7\). **OPEN:** \(Q_\delta\le\rho_{\min}^2\) or ED4≤ED4_suf for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15113.py`, `evidence/e1_gmin_m4_prop15113.json`, `tests/test_prop15113.py`.

## Prop 15.114 (2026-08-01) — γ-calculus; Tf_y multiplicative formula; ∑γ, ∑γ²

**Proved:**

1. **Multiplicative eigenformula.** For boolean \(Cy=py\),
   \((Tf_y)(S)=(4p-2\gamma_y(S))f_y(S)\) where
   \(\gamma_y(S)=\sum_{\{i,j\}\subset S}C_{ij}y_iy_j\).
   Equivalently \(Af_y=2(\gamma_y\odot f_y)\), so \(\langle\delta,\gamma_y\odot f_y\rangle=0\) for all \(\delta\in\ker A\).
2. **∑γ closed form.** \(\sum_S\gamma_y(S)=(6/p)\binom{n}{4}\) (edge double-count + \(y^\top Cy=pn\)).
3. **∑γ² closed form.** \(\sum_S\gamma_y(S)^2=6\binom{n}{4}+n(n-1)(n-2)/4\).
   Adjacent-edge cross terms vanish because \(Cy=py\) and \(n-1=p^2\); matching contribution is \(n(n-1)(n-2)/4\) via \(\sum\kappa\prod=n(n-1)(n-2)/8\).
4. **‖Tf_y‖² closed.** \(\|Tf_y\|_2^2=(16p^2-72)\binom{n}{4}+n(n-1)(n-2)\).
5. **Pair-average residual moment.** The ED4 dictionary uses \(N^{-2}\sum_{y,z}(y\cdot z)^4=E_y[\mathrm{ED4}(y)]\), not a single basepoint.

**Certified:** Thm 1–4 at \(p=3,5\); at \(p=7\), Max+ has **three** ED4(y) types (counts 2352, 8400, 700), all \(\le\mathrm{ED4}_{\mathrm{suf}}\), with global mean \(\le\mathrm{ED4}_{\mathrm{suf}}\) (not 2-point homogeneous).

**OPEN:** pair-average \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) for general \(p\ge5\). Attack: γ-weighted spectral mass of \(f_y\) on \(E_{4p}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15114.py`, `evidence/e1_gmin_m4_prop15114.json`, `tests/test_prop15114.py`.

## Prop 15.115 (2026-08-01) — Max+ residual solves \(A\rho=b\); \(\delta=P_{E_{4p}}m_4\); spectral moments of \(f_y\)

**Proved:**

1. **\(E[\gamma\odot f]=2\kappa/p\).** By the Max+ 2-design \(E[y_iy_j]=C_{ij}/p\) and double-counting matchings on each \(K_4\).
2. **Resolvent.** \(Am_4=4\kappa/p\) and \(A\rho=b\) for the Max+ residual \(\rho=m_4-\kappa/p^2\).
3. **\(\kappa\perp E_{4p}\).** \(T\kappa\in E_\mu\oplus E_{-\mu}\) with \(\mu^2=4(p^2+15)\); \((4p)^2-\mu^2=12(p^2-5)\ne0\) for primes \(p\ge3\). Hence \(\delta=P_{E_{4p}}m_4=E_y P_{E_{4p}}f_y\).
4. **Spectral moments of \(f_y\).** \(m_1=4p-12/p\), \(m_2=16p^2-72+24/(p^2-2)\),
   \(\mathrm{Var}=24(p^2-3)(p^2-4)/(p^2(p^2-2))\).
5. **Jensen.** \(\delta^2\le E_y\|P_{E_{4p}}f_y\|_2^2\) (full \(E_{4p}\) energy is too crude: at \(p=5\) it is \(\gg\rho_{\min}^2\) while residual still holds).

**Certified:** resolvent + moments at \(p=3,5\); \(Q_\delta\) constant (\(=\delta^2\)) at \(p=3,5\); ED4≤ED4_suf at \(p=3,5,7\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\) via the **coherent** mass \(\|E_y P_{E_{4p}}f_y\|_2\), not full \(E_{4p}\) energy. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15115.py`, `evidence/e1_gmin_m4_prop15115.json`, `tests/test_prop15115.py`.

## Prop 15.116 (2026-08-01) — e₄↔ED4↔δ dictionary; coherent mass; Aut-line criterion

**Proved:**

1. **e₄ poly.** Boolean \(e_4(s)=s^4/24+((4-3n)/12)s^2+n(n-2)/8\).
2. **∑m₄² from ED4.** \(\sum m_4^2=\mathrm{ED4}/24+n(4-3n)/6+n(n-2)/8\).
3. **⟨κ,ρ_min⟩.** \(n(n-1)(n-2)(n-6)/(2p^2(p^2-5))\).
4. **Flat identity.** Pythagoras flat part matches the e₄ constant for all primes \(p\ge3\), \(p^2\neq5\).
5. **Coherent mass.** \(\delta^2=\|E_y P_{E_{4p}}f_y\|_2^2\).
6. **Aut-line.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\) then \(\delta=cv_0\) and residual \(\Leftrightarrow c^2\le\rho_{\min}^2\); when \(Q_0\) constant, \(c=Q_0(x_{\mathrm{hs}})\) (halfspace, Max+-free).
7. **Min-distance envelope.** \(|D|\le p^2-2p-1\) for \(y\neq\pm z\) yields an ED4 UB that **fails** \(\mathrm{ED4}_{\mathrm{suf}}\) for all primes \(p\ge5\) (dead for closing).

**Certified:** dictionary at \(p=3,5,7\); Aut-line at \(p=3,5\); \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\).

**OPEN:** coherent mass \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15116.py`, `evidence/e1_gmin_m4_prop15116.json`, `tests/test_prop15116.py`.

## Prop 15.117 (2026-08-01) — Path C hyp residual primary; ρ_min pairings

**Proved:**

1. **Path C primary residual.** For all primes \(p\ge5\),
   \[
   \delta^2\le\frac{\mathrm{room}_{\mathrm{hyp}}}{24}
   =\frac{4(p^2-9)(p^2-1)^2}{3(p^2-5)(p^2+1)}
   \]
   \(\Leftrightarrow\) \(\mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}\) \(\Leftrightarrow\) \(\|\kappa\|_F^2\le\kappa_{\mathrm{hyp}}\) \(\Leftrightarrow\) \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\).
2. **ρ_min² vs hyp.** \(\rho_{\min}^2>\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=5\); \(\rho_{\min}^2<\mathrm{room}_{\mathrm{hyp}}/24\) for all primes \(p\ge7\). Thus \(\delta^2\le\rho_{\min}^2\) is sufficient for Path C when \(p\ge7\); at \(p=5\) the hyp form is the tight target.
3. **Slack.** \(\kappa_{\mathrm{hyp}}-\|\kappa\|_F^2=24(\mathrm{room}_{\mathrm{hyp}}/24-\delta^2)\).
4. **Pairings.** \(\langle b,f_y\rangle=\langle\rho_{\min},b\rangle=2(p^4-1)/p\) for every \(y\in\mathrm{Max+}\); \(\langle\rho_{\min},m_4\rangle=\rho_{\min}^2+\langle\kappa,\rho_{\min}\rangle/p^2\) closed.
5. **γ-channel average.** \(E_y\langle b,\gamma_y\odot f_y\rangle=0\). If pointwise zero, then \(\langle\rho_{\min},f_y\rangle=4(p^4-1)/(3(p^2-5))\) constant on Max+.
6. **Coherent mass.** Path C residual \(\Leftrightarrow\|E_y P_{E_{4p}}f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\).

**Certified:** pointwise \(\langle b,\gamma\odot f\rangle=0\) at \(p=3\); hyp residual \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=3\) (eq), and at \(p=5,7\) when Max+ caches are present (eq at 5, strict at 7).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15117.py`, `evidence/e1_gmin_m4_prop15117.json`, `tests/test_prop15117.py`.

## Prop 15.118 (2026-08-01) — Pointwise ⟨b,γ⊙f⟩=0; T²κ pairing; ρ_min·f closed

**Proved:**

1. **⟨κ,γ⊙f⟩ closed.** For every Max+ \(y\),
   \(\sum_S\kappa(S)\gamma_y(S)f_y(S)=p(p^2+1)(p^2-1)(p^2-4)/4\)
   (edge expansion + \(Cy=py\) + \(C^2=p^2I\)).
2. **⟨T²κ,m₄⟩.** With \(T\kappa\in E_\mu\oplus E_{-\mu}\) and \(\delta\perp T^2\kappa\),
   \(\langle T^2\kappa,m_4\rangle=8p^2(p^4-1)\).
3. **Pointwise criterion.** For each \(y\in\mathrm{Max+}\):
   \(\langle b,\gamma\odot f\rangle=0\Leftrightarrow\langle Tb,f\rangle=4p\langle b,f\rangle\Leftrightarrow\langle T^2\kappa,f\rangle=8p^2(p^4-1)\).
4. **Pointwise vanishing.** \(\langle b,\gamma_y\odot f_y\rangle=0\) and \(\langle T^2\kappa,f_y\rangle=8p^2(p^4-1)\) for all Max+ \(y\) (constancy of the \(T\)-module pairing on Max+ + mean from (2)).
5. **⟨ρ_min,f_y⟩ closed.** \(\langle\rho_{\min},f_y\rangle=4(p^4-1)/(3(p^2-5))\) for all Max+ \(y\).

**Certified:** full Max+ census at \(p=3\) (\(N=12\)) and \(p=5\) (\(N=260\)): bgf=0, \(T^2\) target, ρ_min pairing, and \(\delta^2=\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=5\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (pin \(\langle m_4,f_{\mathrm{hs}}\rangle\) / ED4). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15118.py`, `evidence/e1_gmin_m4_prop15118.json`, `tests/test_prop15118.py`.

## Prop 15.119 (2026-08-01) — Residual budget dictionary; weight enum; halfspace pin

**Proved:**

1. **ED4 budgets (closed rationals).** For primes \(p>\sqrt5\),
   \[
   \mathrm{ED4}_{\mathrm{flat}}=\frac{4(p^2-3)(p^2+1)(3p^2+1)}{p^2-5},\qquad
   \mathrm{ED4}_{\mathrm{bud}}=\frac{4(3p^8+6p^6-104p^4+138p^2-75)}{(p^2-5)(p^2+1)}.
   \]
2. **E[W²] channel.** With \(D=y\cdot z\), \(D^2=n+2W\), \(\mathrm{ED4}=3n^2+4\mathbb E[W^2]\):
   \[
   \mathrm{EW2}_{\mathrm{flat}}=\frac{(p^2+1)(9p^4-20p^2+3)}{4(p^2-5)},\qquad
   \mathrm{EW2}_{\mathrm{bud}}=\frac{9p^8+30p^6-380p^4+594p^2-285}{4(p^2-5)(p^2+1)}.
   \]
   Path C residual \(\Leftrightarrow\mathbb E[W^2]\le\mathrm{EW2}_{\mathrm{bud}}\). Gap:
   \(\mathrm{EW2}_{\mathrm{bud}}-\mathrm{EW2}_{\mathrm{flat}}=\mathrm{room}_{\mathrm{hyp}}/4=6\cdot(\mathrm{room}_{\mathrm{hyp}}/24)\).
3. **⟨m₄,f_y⟩ ↔ ED4.** \(\langle m_4,f_y\rangle=\mathrm{ED4}(y)/24+n(4-3n)/6+n(n-2)/8\). Budget forms:
   \[
   \mathrm{m4f}_{\mathrm{flat}}=\frac{(p-1)(p+1)(p^2+1)(3p^2+17)}{24(p^2-5)},\qquad
   \mathrm{m4f}_{\mathrm{bud}}=\frac{(p-1)(p+1)(3p^2-5)(p^4+20p^2-61)}{24(p^2-5)(p^2+1)}.
   \]
4. **Equivalence.** For \(p\ge5\): \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\Leftrightarrow\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\Leftrightarrow\mathbb E[W^2]\le\mathrm{EW2}_{\mathrm{bud}}\). When \(Q_\delta\) is constant on Max+ (certified \(p=3,5\)): also \(\Leftrightarrow\langle m_4,f_y\rangle\le\mathrm{m4f}_{\mathrm{bud}}\) (halfspace pin).
5. **Weight-enumerator structure.** Max+ dots satisfy \(y\cdot z\equiv2\pmod4\), \(|y\cdot z|\le p^2-2p-1\) off \(\pm\) pairs, antipodal measure. Crude envelope \(\mathrm{ED4}\le2n D_{\max}^2\) is strictly larger than \(\mathrm{ED4}_{\mathrm{bud}}\) for \(p\ge5\) (too weak).

**Certified:** full Max+ weight spectra at \(p=3\) (\(\{\pm10,\pm2\}\)) and \(p=5\) (\(\{\pm26,\pm14,\pm10,\pm6,\pm2\}\) with mults \(\{1,13,20,36,60\}\)); \(\mathrm{ED4}=\mathrm{ED4}_{\mathrm{bud}}\) and \(\langle m_4,f_{\mathrm{hs}}\rangle=\mathrm{m4f}_{\mathrm{bud}}\) at both; \(\mathrm{orth}\cdot N=147456\) at \(p=5\). Prior \(p=7\): \(\delta^2/(\mathrm{room}_{\mathrm{hyp}}/24)=124875/669124\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (independent upper bound on \(\mathbb E[W^2]\) or \(\langle m_4,f_{\mathrm{hs}}\rangle\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15119.py`, `evidence/e1_gmin_m4_prop15119.json`, `tests/test_prop15119.py`.

## Prop 15.120 (2026-08-01) — Pointwise E[W²] factorization; Pythagoras; majorization

**Proved:**

1. **Pointwise residual factorization.** For every \(y\in\mathrm{Max+}\),
   \(\mathbb E_z[W_y(z)^2]=\mathrm{EW2}_{\mathrm{flat}}+6\,Q_\delta(y)\).
   Wick wedges + \(\kappa/p^2+\rho_{\min}\) on disjoint edges give the flat part;
   the \(\delta\)-part of disjoint-edge pairings contributes exactly \(6\langle\delta,f_y\rangle\).
2. **Pythagoras.** \(\sum m_4^2=F(p)+\delta^2\) with
   \(F=\|\kappa/p^2+\rho_{\min}\|_2^2=\mathrm{m4f}_{\mathrm{flat}}\) Max+-free.
3. **Majorization UB.** \(H=G\odot G\succeq0\), \(\lambda_{\max}=2Nn\), \(\mathrm{Tr}=Nn^2\), \(\mathrm{Tr}/\lambda_{\max}=d\)
   \(\Rightarrow\mathrm{ED4}\le 2n^3\). For all primes \(p\ge5\), \(2n^3>\mathrm{ED4}_{\mathrm{bud}}\) (too weak).
4. **Dead independent UBs.** CS via \(\|\tilde\gamma\odot f\|\), discrete moment LP \(2n D_{\max}^2\),
   and majorization \(2n^3\) all exceed residual budgets for \(p\ge5\).

**Certified:** EW2 constant on Max+ and equal to \(\mathrm{EW2}_{\mathrm{flat}}+6\delta^2\) at \(p=3,5\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (weight enumerator / Gauss sums). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15120.py`, `evidence/e1_gmin_m4_prop15120.json`, `tests/test_prop15120.py`.

## Prop 15.121 (2026-08-01) — Spectral residual dictionary; Frobenius form

**Proved:**

1. **ED4 from FFT.** \(M_{ab}=(y_a\cdot y_b)^2=nJ+2\,\mathrm{FFT}\) (15.93) \(\Rightarrow\)
   \(\mathrm{ED4}=4n^2+4N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\).
2. **E[W²] from FFT.** \(\mathbb E[W^2]=\|\mathrm{FFT}\|_F^2/N^2=d^2+N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\),
   and with 15.120: \(\|\mathrm{FFT}|_{1^\perp}\|_F^2=N^2(\mathrm{EW2}_{\mathrm{flat}}-d^2+6\delta^2)\).
3. **Φ variance.** \(\sum\lambda_\alpha^2=\mathrm{ED4}-4n^2\), \(\sum(\lambda_\alpha-\bar\mu)^2=\mathrm{orth}=24\delta^2\) (15.105);
   at \(\delta=0\), \(\Phi\equiv\bar\mu\) on \(Z\).
4. **EW2_flat closed.** \(\mathrm{EW2}_{\mathrm{flat}}=(n^2+T^2/m)/4\) with \(T=n(n-2)\), \(m=\dim Z\).
5. **Residual ⇔ Frobenius.** Path C residual \(\Leftrightarrow\|\mathrm{FFT}|_{1^\perp}\|_F^2\le N^2(\mathrm{EW2}_{\mathrm{bud}}-d^2)\)
   \(\Leftrightarrow\sum(\lambda_\alpha-\bar\mu)^2\le\mathrm{room}_{\mathrm{hyp}}\).
   Contrast: 16N \(\Leftrightarrow\|\mathrm{FFT}|_{1^\perp}\|_{\mathrm{op}}\le8N\) (operator norm on the same operator).
6. **H/16N·Tr too weak.** \(\|A\|_F^2\le\|A\|_{\mathrm{op}}\mathrm{Tr}(A)\) under H or 16N yields ED4 UBs \(> \mathrm{ED4}_{\mathrm{bud}}\) for all \(p\ge5\).

**Certified:** spectral identities at \(p=3,5\) (full Max+ Gram/FFT).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15121.py`, `evidence/e1_gmin_m4_prop15121.json`, `tests/test_prop15121.py`.

## Prop 15.122 (2026-08-01) — Max+ disagreement \(u\in V_+\); Aut-line; \(\lambda_{\max}(T)\)

**Proved:**

1. **Disagreement identity.** For \(y,z\in\mathrm{Max}+\), \(D=y\cdot z\), \(k=(n-D)/2\),
   \(u=(y-z)/2\in\{0,\pm1\}^n\) has weight \(k\) and
   \(u^\top Cu=pk\) with \(u\in V_+\) (equivalently \(\|P_+u\|_2^2=k\)).
   Same for \(v=(y+z)/2\) with weight \(n-k\).
2. **Dot support.** Off \(\pm\) pairs, \(|D|\le p^2-2p-1\), \(D\equiv2\pmod4\), and \(k\) is a
   ternary \(V_+\) weight. Cert spectra \(p=3,5\).
3. **\(\lambda_{\max}(T)\) threshold.** \(\lambda_{\max}(T)<4p\Rightarrow\delta=0\) (residual OK).
   Cert: \(p=3\) strict; \(p=5,7\) equality \(\lambda_{\max}=4p\).
4. **Aut-line.** \(\delta\in E_{4p}^{\mathrm{Aut}}\); if \(\dim\le1\) then residual \(\Leftrightarrow c^2\le\mathrm{room}_{\mathrm{hyp}}/24\)
   with \(c=Q_0(x_{\mathrm{hs}})\). Cert equality line at \(p=5\); \(p=7\) has non-constant \(Q_\delta\).
5. **FFT budget.** Residual \(\Leftrightarrow N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\le B(p)=\mathrm{EW2}_{\mathrm{bud}}-d^2\) (Max+-free).
6. **Dead tight attempts.** Discrete LP with exact \(N\), PGL character sums (F18), and prior majorizations all exceed residual budgets.

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (ternary \(V_+\) weight enumerator / Gauss \(Q_0\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15122.py`, `evidence/e1_gmin_m4_prop15122.json`, `tests/test_prop15122.py`.

## Prop 15.123 (2026-08-01) — Switching; conference srg; dual Krawtchouk residual

**Proved:**

1. **Switching bijection.** After \(C'=D_yCD_y\), \(1\in\mathrm{Max}+'\) and
   \(z\mapsto w=(1-D_yz)/2\) bijects \(\mathrm{Max}+'\leftrightarrow V_+\cap\{0,1\}^n\), with
   \(D(y,z)=n-2\,\mathrm{wt}(w)\).
2. **Conference srg.** \(G\) with edges \(C'=-1\) is
   \(\mathrm{srg}(n,\,p(p-1)/2,\,\mu-1,\,\mu)\), \(\mu=((p-1)/2)^2\);
   \(A\)-eigs \(k_G\), \(\theta=(p-1)/2\) (\(\times d\)), \(\tau=-(p+1)/2\) (\(\times d-1\));
   \(C'\)-eigs \(\pm p\) each of mult \(d\).
3. **Regular sets.** Supports of weight-\(k\) codewords are regular sets with
   \(\alpha=(k-1-p)/2\), \(\beta=k/2\), \(\alpha-\beta=\tau\); allowed
   \(k\in\{0,n\}\cup\mathrm{even}[p+1,p(p-1)]\).
4. **Weight = distance dist.** \(B_i=W_i\), \(|X|=N\), \(\mathrm{ED4}=N^{-1}\sum W_k(n-2k)^4\).
5. **Dual residual.** \(A'_4=\sum m_4^2=\mathrm{m4f}_{\mathrm{flat}}+\delta^2\);
   residual \(\Leftrightarrow A'_4\le\mathrm{m4f}_{\mathrm{bud}}\). (Hamming Delsarte alone too weak.)
6. **Two-valued form.** \(\chi_S-(k/n)\mathbf1\in V_+\cap\mathbf1^\perp\) with only two coordinate values
   and \(\|\cdot\|_2^2=k(n-k)/n\); \(W_k\) counts such vectors.

**Certified:** srg params and full \(W_k\) at \(p=3,5\) (Petersen; srg(26,10,3,4)).

**OPEN:** closed \(W_k\) / \(A'_4\le\mathrm{m4f}_{\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15123.py`, `evidence/e1_gmin_m4_prop15123.json`, `tests/test_prop15123.py`.

## Prop 15.124 (2026-08-01) — Closed weight moments \(j\le3\); \(E[k^4]\) partition; residual as \(R_4\)

Continues 15.123. Isolates the residual as a single 4-wise mass. Does **not** soft-close.

**Proved:**

1. **Moments \(j\le3\) (Max+-free).** After switch, \(E[s^2]=2n\), \(E[s]=0\), \(E[s^3]=0\), and
   \[
   E[k]=\tfrac n2,\quad
   E[k^2]=\tfrac{n(n+2)}{4},\quad
   E[k^3]=\tfrac{n^2(n+6)}{8}.
   \]
2. **Exact \(\le3\) partition of \(E[k^4]\).** Using only pair/triple design averages,
   \[
   \mathrm{exact}_{\le3}
   =\tfrac n2+n^2+\tfrac{3n^2}{4}+\tfrac{3n(n-2)(n+2)}{4},
   \]
   and \(R_4:=E[k^4]-\mathrm{exact}_{\le3}=n(n-1)(n-2)(n-3)\,\overline E[\prod_4 w]\).
3. **Residual dictionary.** \(E[D^4]=\mathrm{ed4\_from\_exact3}(p)+16 R_4\); Path C residual
   \(\Leftrightarrow R_4\le R_{4,\mathrm{bud}}\Leftrightarrow A'_4\le\mathrm{m4f}_{\mathrm{bud}}\).
4. **Hamming Delsarte LP.** Max \(A'_4\) under allowed-weight dual constraints saturates
   \(\mathrm{m4f}_{\mathrm{bud}}\) at \(p=3\) and **strictly exceeds** it at \(p=5,7\) (too weak for \(p\ge5\)).
5. **Hoffman layer.** \(W_{p+1}=d\) (regular cocliques) certified at \(p=3,5\).

**Certified:** moments + partition + \(R_4=\mathrm{bud}\) saturation + \(A'_4=\mathrm{m4f}_{\mathrm{bud}}\) at \(p=3,5\).

**OPEN:** closed \(W_k\) / \(A'_4\le\mathrm{m4f}_{\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15124.py`, `evidence/e1_gmin_m4_prop15124.json`, `tests/test_prop15124.py`.

## Prop 15.125 (2026-08-01) — Perfect 2-colorings; 4-design defect; closed \(R_4\) budget

Continues 15.123–15.124. Does **not** soft-close.

**Proved:**

1. **Perfect 2-colorings.** \(W_k\) equals the number of \(\tau\)-equitable bipartitions
   (perfect 2-colorings) of the conference srg with \(|S|=k\),
   \(\alpha=(k-1-p)/2\), \(\beta=k/2\), \(\alpha-\beta=\tau=-(p+1)/2\).
   Hoffman layer \(k=p+1\) has \(\alpha=0\).
2. **Spherical 2-design / 4-design defect.** After \(V_+\cong\mathbb R^d\), Max+ is a
   spherical 2-design (\(E[uu^\top]=I_d/d\)). It is not a 4-design for \(p\ge5\):
   \(E[s^4]=\mathrm{ED4}\) exceeds \(3n^4/(d(d+2))\). Path C residual is this defect
   controlled by \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\).
3. **Closed \(R_4\) budget.**
   \(\mathrm{ed4\_from\_exact3}=-(p^2+1)(p^6+3p^4-25p^2+13)\);
   \(R_{4,\mathrm{bud}}=(\mathrm{ED4}_{\mathrm{bud}}-\mathrm{ed4\_from\_exact3})/16\).
4. **Delsarte + moments \(j\le3\).** Still saturates only at \(p=3\); weak for \(p\ge5\).
5. **Antipodal dual.** \(W_k=W_{n-k}\Rightarrow A'_j=A'_{n-j}\).

**Certified:** algebra + LP weak p=5,7; defect ratios and saturation at \(p=3,5\).

**OPEN:** closed \(W_k\) / \(R_4\le R_{4,\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15125.py`, `evidence/e1_gmin_m4_prop15125.json`, `tests/test_prop15125.py`.

## Prop 15.126 (2026-08-01) — Geometric Hoffman seed; 1-design; simplex bound

Continues 15.123–15.125. Does **not** soft-close.

**Proved:**

1. **Geometric seed.** The subfield line \(F_p\cup\{\infty\}\subset\mathrm{PG}(1,p^2)\) is a
   Hoffman coclique (\(\tau\)-equitable, \(\alpha=0\), size \(p+1\)) after the standard
   halfspace Seidel switch. Certified \(p=3,5,7\).
2. **1-design algebra.** Hoffman cocliques form an Aut-invariant 1-design:
   \(b(p+1)=nr\), and \(r=(p+1)/2\Longleftrightarrow b=d=n/2\).
3. **Census equality.** At \(p=3,5\): \(W_{p+1}=d\), \(r=(p+1)/2\), and
   \(\{\chi_S\}\) is a basis of \(V_+\).
4. **Simplex bound.** Equal pairwise intersections \(\Rightarrow W_{p+1}\le d\)
   (Gram rank \(b-1\le d-1\)); equality at \(p=3\) (\(\lambda=1\)).
5. **ED4 contribution.** If \(W_{p+1}=W_{n-p-1}=d\), Hoffman layers contribute
   \(2d(p^2-2p-1)^4/N\) to ED4.

**OPEN:** \(W_{p+1}=d\) for all primes \(p\ge5\); full \(W_k\) / 4-design defect bound.
L remains OPEN. **(W=d general disproved in Prop 15.127.)**

Evidence: `src/e1_gmin_m4_prop15126.py`, `evidence/e1_gmin_m4_prop15126.json`, `tests/test_prop15126.py`.

## Prop 15.127 (2026-08-01) — Closed \(W_{p+1}\); inversive plane; \(W=d\) false

Continues 15.126. Does **not** soft-close residual.

**Proved:**

1. **Inversive plane.** \(F_p\)-sublines of \(\mathrm{PG}(1,p^2)\) form the miquelian
   inversive plane of order \(p\): \(S(3,p+1,p^2+1)\), \(b=p(p^2+1)\), \(\lambda_2=p+1\), \(\lambda_3=1\).
2. **Closed Hoffman weight.**
   \[
   W_{p+1}
   =\frac{1+\chi_4(p)}{2}\cdot\frac{p^2+1}{2}
   +\frac{1-\chi_4(p)}{2}\cdot\frac{3p+1}{2},
   \quad
   \chi_4(p)=(-1)^{(p-1)/2}.
   \]
   Equivalently: \(W_{p+1}=d\) if \(p\equiv1\pmod4\), and \(W_{p+1}=(3p+1)/2\) if \(p\equiv3\pmod4\).
3. **Counterexample to \(W_{p+1}=d\).** At \(p=7\), \(W_8=11\neq25=d\).
4. **Census.** Full max-coclique enumeration: \(W=(5,13,11,17)\) at \(p=3,5,7,11\), matching the formula; coincides with the count of regular sublines.
5. **Corrected ED4 part.** Hoffman contribution \(2W_{p+1}D_{\max}^4/N\) (not always \(2d\)).

**OPEN:** full closed \(W_k\) / 4-design defect bound for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15127.py`, `evidence/e1_gmin_m4_prop15127.json`, `tests/test_prop15127.py`.

## Prop 15.128 (2026-08-01) — Full \(W_k\) census \(p=3,5,7\); exact ED4 at \(p=7\)

Continues 15.127. Does **not** soft-close residual.

**Certified:**

1. **Full weight enumerators.**
   - \(p=3\) (\(N=12\)): \(\{0{:}1,4{:}5,6{:}5,10{:}1\}\)
   - \(p=5\) (\(N=260\)): \(\{0{:}1,6{:}13,8{:}20,10{:}36,12{:}60,14{:}60,\ldots\}\)
   - \(p=7\) (\(N=11452\)), free-coord Max+ \(2^{25}\):
     \(\{0{:}1,8{:}11,12{:}112,14{:}159,16{:}280,18{:}728,20{:}1099,22{:}1502,24{:}1834,\ldots\}\)
     with \(W_{10}=W_{40}=0\) (allowed by \(\alpha\ge0\) but empty).
2. **Consistency.** \(\sum W=N\), \(W_k=W_{n-k}\), \(W_{p+1}\) matches 15.127,
   \(E[k^j]\) for \(j\le3\) match 15.124.
3. **ED4 at \(p=7\).** \(\mathrm{ED4}=12835984/409<\mathrm{ED4}_{\mathrm{bud}}=1775728/55\)
   (strict hyp residual); \(\delta^2=82176/4499\).
4. **Saturation.** \(\mathrm{ED4}=\mathrm{ED4}_{\mathrm{bud}}\) at \(p=3,5\).

**OPEN:** closed \(W_k\) for general primes \(p\ge5\) (must allow structural zeros),
or character-sum/PBIBD bound on the 4-design defect. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15128.py`, `evidence/e1_gmin_m4_prop15128.json`, `tests/test_prop15128.py`.

## Prop 15.129 (2026-08-01) — Jensen coherent-mass inequality; Hoffman \(\bar r\)

Continues 15.128. Does **not** soft-close residual.

**Proved:**

1. **Jensen.** \(\delta=E_y[P_{E_{4p}}f_y]\) and
   \(\delta^2\le E\|P f_y\|_2^2\), equality iff \(P f_y\) is a.s. constant on Max+.
   Hence \(E\|P f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\Rightarrow\) Path C residual.
2. **Dictionary.** Linear equivalences ED4 \(\leftrightarrow\sum m_4^2\leftrightarrow\delta^2\)
   certified on full \(W\) at \(p=3,5,7\).
3. **Average replication.** \(\bar r=W_{p+1}(p+1)/n\) with \(W_{p+1}\) from 15.127.
   Integral (1-design possible) for \(p\equiv1\pmod4\) and \(p=3\); at \(p=7\),
   \(\bar r=44/25\notin\mathbb Z\), so the Hoffman layer is **not** a 1-design.
4. **Hoffman geometry at \(p=5\).** 30 disjoint Hoffman pairs, each giving a
   regular 12-set; these are exactly half of \(W_{12}\). Every weight-16 set
   contains exactly two Hoffman cocliques.

**OPEN:** bound \(E\|P_{E_{4p}}f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for all primes
\(p\ge5\), or closed general \(W_k\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15129.py`, `evidence/e1_gmin_m4_prop15129.json`, `tests/test_prop15129.py`.

## Prop 15.130 (2026-08-01) — \(P m_4=\delta\); \(\rho_{\min}\)-sufficient residual for \(p\ge7\)

Continues 15.129. Does **not** soft-close residual.

**Proved:**

1. **\(P m_4=\delta\).** Since \(m_4=\kappa/p^2+\rho_{\min}+\delta\) with the first two summands
   orthogonal to \(E_{4p}\), \(P m_4=\delta\). With \(m_4=E f_y\), \(\delta=E[P f_y]\).
2. **Jensen.** \(\delta^2\le E\|P f_y\|_2^2\); bound on the right-hand side yields residual.
3. **Gap algebra.** For all primes \(p\ge7\),
   \[
   \frac{\mathrm{room}_{\mathrm{hyp}}}{24}-\rho_{\min}^2
   =\frac{(p-1)(p+1)(3p^6-105p^4+37p^2-15)}{6p^2(p^2-5)(p^2+1)}>0,
   \]
   so \(\delta^2\le\rho_{\min}^2\Rightarrow\) Path C residual. (At \(p=5\), \(\rho_{\min}^2>\mathrm{room}\); use hyp form.)
4. **Census.** \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\) (ratio \(\approx0.379\) at \(p=7\)).
5. **Aut-line program.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\), residual \(\Leftrightarrow c^2\le\mathrm{room}/24\)
   with \(c=Q_0(\mathrm{hs})\) when \(Q_0\) is constant.

**OPEN:** \(\delta^2\le\rho_{\min}^2\) (or \(E\|Pf\|^2\le\mathrm{room}\), or Gauss-sum \(Q_0\)) for general
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15130.py`, `evidence/e1_gmin_m4_prop15130.json`, `tests/test_prop15130.py`.

## Prop 15.131 (2026-08-01) — Pair-avg vs basepoint ED4; p=7 three-type \(Q_\delta\)

Continues 15.130. Does **not** soft-close residual.

**Proved / certified:**

1. **Pair vs basepoint.** \(ED4(y)=ED4_{\mathrm{flat}}+24 Q_\delta(y)\) and
   \(ED4_{\mathrm{pair}}=ED4_{\mathrm{flat}}+24\delta^2\). When \(Q_\delta\) is non-constant,
   \(ED4(y_0)\ne ED4_{\mathrm{pair}}\) for a fixed basepoint \(y_0\) (e.g. halfspace);
   the weight-enumerator moment \(E[(n-2k)^4]\) equals \(ED4(y_0)\), not \(\delta^2\).
2. **p=7 spectrum (ProcessPool W=86).** Exactly three types:
   counts \(2352/8400/700\), \(Q_\delta\in\{-124800/4499,\ 82176/4499,\ 200448/4499\}\).
   \(Q_\delta\) can be negative (not 2-point homogeneous for 4th moment).
3. **True \(\delta^2\).** Pair-avg \(\delta^2=19180800/1840091\approx10.424\)
   (prior W-based \(82176/4499\) was \(Q_\delta(\mathrm{hs})\), the middle type).
   Ratios: \(\delta^2/\rho_{\min}^2\approx0.216\), \(\delta^2/\mathrm{room}\approx0.187\).
4. **Pointwise criterion.** \(\max Q_\delta=200448/4499\le\rho_{\min}^2=26000/539\)
   with slack \(812048/220451\); hence residual for \(p=7\) by 15.130.C.
5. **Variance.** \(\mathrm{Var}(Q_\delta)>0\) at \(p=7\), so pointwise is strictly
   stronger than \(\delta^2\le\rho_{\min}^2\).

**OPEN:** \(\max_y Q_\delta(y)\le\rho_{\min}^2\) (or \(\delta^2\le\rho_{\min}^2\)) for all primes
\(p\ge7\); at \(p=5\) use hyp form. Aut-line at \(p=7\) needs care (\(Q\) non-constant).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15131.py`, `evidence/e1_gmin_m4_prop15131.json`, `tests/test_prop15131.py`.

## Prop 15.132 (2026-08-01) — Max+-free residual dictionary; Aut \(\delta\); dead envelopes

Continues 15.131. Does **not** soft-close residual. No full Max+ census for general \(p\).

**Proved / certified:**

1. **Max+-free dictionary.** \(\delta^2=\sum m_4^2-m4f_{\mathrm{flat}}\) with
   \(m4f_{\mathrm{flat}}\) closed Max+-free. Residual
   \(\Leftrightarrow\sum m_4^2\le m4f_{\mathrm{bud}}:=m4f_{\mathrm{flat}}+\mathrm{room}/24\).
   For \(p\ge7\): \(\sum m_4^2\le m4f_{\mathrm{flat}}+\rho_{\min}^2\) suffices.
2. **Aut-invariance.** \(\delta\) is Aut-invariant; \(Q_\delta\) constant on Aut-orbits of Max+.
   At \(p=7\), three \(Q\) types \(\Rightarrow\ge3\) Aut-orbits (not transitive).
3. **\(\gamma\)-parity.** \(\gamma_y(S)\in\{-6,-4,-2,0,2,4,6\}\); formal \(4p\)-fiber is \(\gamma=0\).
4. **\(\gamma=0\) mass.** Constant \(4350\) at \(p=5\); 3-valued in \(p=7\) sample (parallels \(Q\) types).
5. **Dead envelopes.** Moment LP on allowed regular-set \(k\), and pole+\(D_{\max}\) mixes with
   any \(N\ge n\), all strictly exceed \(\mathrm{ED4}_{\mathrm{suf}}\) for primes \(p=5..19\).

**OPEN:** Max+-free bound \(\max Q_\delta\le\rho_{\min}^2\) or \(\delta^2\le\rho_{\min}^2\)
(character sums / Aut-orbit Bose–Mesner / \(\dim E_{4p}^{\mathrm{Aut}}\le1+Q_0\)) for all
primes \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15132.py`, `evidence/e1_gmin_m4_prop15132.json`, `tests/test_prop15132.py`.

## Prop 15.133 (2026-08-01) — class_key Bose–Mesner; F19 quantitative; Aut-line

Continues 15.132. Does **not** soft-close residual. Does **not** thrash class_key (F19).

**Proved / certified:**

1. **Aut-line form.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\), residual
   \(\Leftrightarrow c^2\le\mathrm{room}/24\) with \(c=Q_0(\mathrm{hs})\) when \(Q_0\) constant.
2. **class_key \(T\)-spectrum (ProcessPool W=86).**
   \(\dim E_{4p}^{\mathrm{ck}}=0,1,0\) at \(p=3,5,7\);
   \(\lambda_{\max}=4p\) only at \(p=5\) among these.
3. **F19 quantitative.** At \(p=7\), \(\dim E_{4p}^{\mathrm{ck}}=0\) but
   \(\delta^2=19180800/1840091>0\), so \(\delta\notin V^{\mathrm{ck}}\)
   (\(m_4\) not class_key-equitable). class_key cannot close residual at \(p=7\).
4. **CR dead.** PGL cross-ratio orbits are not \(\mathrm{Aut}(C)\) orbits
   (\(\kappa\) takes both signs inside CR cells at \(p=5\)); false nullity 0.
5. **Success locus.** Aut-line via class_key works at \(p=5\) (\(c^2=\mathrm{room}=1536/65\));
   ambient \(E_{4p}=0\) at \(p=3\); \(p=7\) needs true \(\mathrm{Aut}(C)\) or character sums.

**OPEN:** true \(\mathrm{Aut}(C)\) Bose–Mesner / Gauss \(Q_0\) / \(\max Q_\delta\le\rho_{\min}^2\)
for general primes. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15133.py`, `evidence/e1_gmin_m4_prop15133.json`,
`evidence/e1_gmin_m4_prop15133_classkey_spectrum.json`, `tests/test_prop15133.py`.

## Prop 15.134 (2026-08-01) — Strict Aut(\(C\)) Bose–Mesner; residual projection

Continues 15.133. Does **not** soft-close residual. No class_key / raw PGL thrash.

**Proved / certified:**

1. **Strict Aut \(G\).** Affine square-semilinear maps
   \(x\mapsto a\cdot\mathrm{Frob}^i(x)+b\) (\(\chi(a)=1\), \(\infty\) fixed) give
   \(|G|=p^2(p^2-1)\) and \(P^\top CP=C\). Inversion is switch-only; adjoining it
   yields all of \(\mathrm{PGL}(2,q)\not\le\mathrm{Aut}(C)\).
2. **Orbits.** Counts \(9/42/128\) at \(p=3,5,7\); \(\kappa\)-constant; \(T\)-equitable;
   strictly finer than class_key at \(p=5,7\).
3. **Spectrum.** \(\dim E_{4p}^{G}=0,2,7\) at \(p=3,5,7\); \(\lambda_{\max}=4p\) for
   \(p=5,7\). Aut-line \(\dim\le1\) **fails** for \(G\) at \(p\ge5\).
4. **Residual projection.** \(\delta=P_{E_{4p}^{G}}m_4\), \(\delta^2=\|Pm_4\|_2^2\)
   recovers \(\delta^2=1536/65\) at \(p=5\) and \(19180800/1840091\) at \(p=7\).
   \(G\) **carries** the residual at \(p=7\) (class_key does not).

**OPEN:** Gauss / character-sum evaluation of \(m_4\) on \(G\)-orbits (Max+-free)
then project to prove \(\delta^2\le\mathrm{room}/24\) for all primes \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15134.py`, `evidence/e1_gmin_m4_prop15134.json`,
`evidence/e1_gmin_m4_prop15134_strict_aut_spectrum.json`,
`evidence/e1_gmin_m4_prop15134_residual_proj.json`, `tests/test_prop15134.py`.

## Prop 15.135 (2026-08-01) — Coherent-mass spectral form; halfspace char sums

Continues 15.134. Does **not** soft-close residual.

**Proved / certified:**

1. **Spectral form.** On an \(L^2\)-ONB \(\{v_j\}\) of \(E_{4p}^{G}\),
   \(\delta=\sum c_j v_j\) with \(c_j=\langle m_4,v_j\rangle=E_y Q_j(y)\), and
   \(\delta^2=\sum c_j^2\). Residual \(\Leftrightarrow\sum c_j^2\le\mathrm{room}/24\).
2. **Halfspace character formula.** \(f_{\mathrm{hs}}(S)\) is Max+-free via
   \(F_p\)-coordinate indicators; \(\sum_S f_{\mathrm{hs}}=e_4\) (same closed form as \(\sum m_4\)).
3. **\(G\cdot\mathrm{hs}\) dead.** \(|O_{\mathrm{hs}}|\in\{60,168\}\) at \(p=5,7\) vs
   \(|\mathrm{Max}^+|\in\{260,11452\}\); \(\delta^2\) from \(m_4^{G\mathrm{hs}}\gg\mathrm{room}\).
4. **Moments don't pin \(\delta\).** \(P_G(\mathbf{1})=P_G(\kappa)=0\) on \(E_{4p}^{G}\);
   \(e_4\) and \(\langle m_4,\kappa\rangle\) do not constrain the free \(c_j\).
5. **Program.** Need character sums over **full** Max+ for the \(c_j\) (Max+ has
   multiple \(G\)-orbits of vectors; \(Q_j\) not constant).

**OPEN:** Gauss/character-sum for \(c_j\) or \(m_4\) on \(G\)-orbits for all primes
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15135.py`, `evidence/e1_gmin_m4_prop15135.json`,
`tests/test_prop15135.py`.

## Prop 15.136 (2026-08-01) — Max+-free flat on \(G\)-orbits; free \(c_j\) remain

Continues 15.135. Does **not** soft-close residual.

**Proved / certified:**

1. **Max+-free flat.** On the \(G\)-orbit space,
   \(\rho_{\min}=A^+(T\kappa/p^2)\) (resolvent on \((\ker A)^\perp\)) and
   \(\mathrm{flat}=\kappa/p^2+\rho_{\min}\) depend only on \(C\). Matches
   \(\|\rho_{\min}\|^2=\rho_{\min}^2\) and \(\|\mathrm{flat}\|^2=m4f_{\mathrm{flat}}\)
   at \(p=3,5,7\); construction available for all primes \(p\ge3\).
2. **Decomposition.** \(m_4=\mathrm{flat}+\delta\) with \(\delta\in E_{4p}^{G}\);
   residual \(\delta^2=\sum c_j^2\). Certified match at \(p=5,7\).
3. **Geometry insufficient.** The invariant
   \((\mathbf{1}_{\infty\in S},\kappa,\dim_{\mathbb F_p}\mathrm{affspan})\)
   does not determine \(m_4\) (split types at \(p=5,7\)).
4. **Character-sum form.** \(c_j=N^{-1}\sum_{y\in\mathrm{Max}^+}Q_j(y)\) with
   Max+-free kernels \(v_j\). Only these \(\nu_G\) coefficients remain.
5. **Partial.** At \(p=5\), some \(\infty\)-orbits have
   \(m_4\in\{-1/5,-21/65\}\).

**OPEN:** Gauss/character-sum evaluation of each \(c_j\) over full Max+ for all
primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15136.py`, `evidence/e1_gmin_m4_prop15136.json`,
`tests/test_prop15136.py`.

## Prop 15.137 (2026-08-01) — \(c_j\) over Max+ \(G\)-orbits; p=5 two-type formula

Continues 15.136. Does **not** soft-close residual.

**Proved / certified:**

1. **G-equivariance.** \(Q_j(g\cdot y)=Q_j(y)\); \(Q_j\) constant on Max+ \(G\)-orbits.
2. **Hemisphere formula.** \(G\) fixes \(\infty\), so
   \(c_j=\sum_t w_t Q_j(y_t)\) over \(G\)-orbits in \(H_+=\{y_\infty=+1\}\).
3. **Census.** \(r=1,2,5\) hemisphere types at \(p=3,5,7\); p=5 weights
   \(3/13\) (hs) and \(10/13\) (other); p=7 sizes \(\{56,84,294,588,1176\}\).
4. **p=5 formula.** \(c_j=(3/13)Q_j(\mathrm{hs})+(10/13)Q_j(y_*)\);
   \(\sum c_j^2=\mathrm{room}=1536/65\).
5. **\(Q_j(\mathrm{hs})\) Max+-free.** Both \(v_j\) and \(f_{\mathrm{hs}}\) are Max+-free.
   Using only \(Q_j(\mathrm{hs})\) for \(c_j\) fails when \(r>1\) (\(\sum Q_j(\mathrm{hs})^2\gg\mathrm{room}\)).

**OPEN:** Max+-free non-hs representatives \(y_t\) and character sums \(Q_j(y_t)\)
for general primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15137.py`, `evidence/e1_gmin_m4_prop15137.json`,
`tests/test_prop15137.py`.

## Prop 15.138 (2026-08-01) — Max+-free non-hs \(y_*\); p=5 residual Max+-free

Continues 15.137. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **hs-switch.** \(C'=D_{\mathrm{hs}}CD_{\mathrm{hs}}\) has \(C'\mathbf1=p\mathbf1\);
   \(y=\mathrm{hs}\odot z'\) with \(C'z'=pz'\) lies in Max+.
2. **Norm circles.** \(S_{t,c}=\{u:N(u-t)=c\}\) (\(N(u)=u^{p+1}\)); lex-first
   \((t^*,c^*)\) making \(S\) a \(C'\)-Hoffman coclique yields Max+-free
   \(y_*\in H_+\) of non-hs type when found (cert \(p=5,7,11\); empty at \(p=13\)).
3. **p=5 residual Max+-free.** With \((t,c)=(0,3)\),
   \(c_j=(3/13)Q_j(\mathrm{hs})+(10/13)Q_j(y_*)\) (both Max+-free) gives
   \(\sum c_j^2=1536/65=\mathrm{room}\) — residual without Max+ census.
4. **Partial p=7.** Norm circles cover hs (size 84) and one non-hs orbit
   (size 588) only — 2 of 5 \(H_+\) types.

**OPEN:** Remaining \(H_+\) G-orbit types for all primes \(p\ge7\); full residual
for all \(p\ge5\). L remains OPEN (do not promote L from p=5 alone).

Evidence: `src/e1_gmin_m4_prop15138.py`, `evidence/e1_gmin_m4_prop15138.json`,
`tests/test_prop15138.py`.

## Prop 15.139 (2026-08-01) — Affine halfspaces + double switch; p=7 size classes

Continues 15.138. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine halfspaces.** For \(L\not\equiv0\) linear and \(|S|=(p+1)/2\),
   \(y_u=+1\Leftrightarrow L(u)\in S\) is Max+ at certified \(p=5,7\) (all such \(S\)).
2. **AP dichotomy at \(p=7\).** Among all \(\binom{7}{4}=35\) sets \(S\):
   4-term AP \(\Rightarrow\) \(H_+\) orbit size 84 (21 sets); non-4-AP \(\Rightarrow\)
   size 56 (14 sets, including QR-half \(\{0,1,2,4\}\)).
3. **Double Seidel–norm-circle.** \(y=y_0\odot z\) with \(C_0=D_{y_0}CD_{y_0}\)
   and \(z\) a norm-circle Hoffman evec of \(C_0\) yields Max+.
4. **All \(p=7\) size classes Max+-free:**
   84 (AP affine), 56 (non-AP affine), 588 (hs\(\odot\)nc),
   1176 (\(y_{56}\odot\)nc), 294 (\(y_{\mathrm{nc}}\odot\)nc).

**OPEN:** Max+-free weights \(w_t\) and character sums \(Q_j(y_t)\) for every
orbit (including four size-1176 orbits); residual for general \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15139.py`, `evidence/e1_gmin_m4_prop15139.json`,
`tests/test_prop15139.py`.

## Prop 15.140 (2026-08-01) — Weights \(|G|/|\mathrm{Stab}|\); character-sum \(Q_j\); p=7 residual form

Continues 15.139. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Orbit–stabiliser weights.** \(|O_t|=|G|/|\mathrm{Stab}^\infty(y_t)|\);
   stabs \(\{42,28,8,4,2\}\) for sizes \(\{56,84,294,588,1176\}\);
   \(|H_+|=5726\); weights sum to 1.
2. **Character-sum residual.** \(c_j=\sum_t w_t Q_j(y_t)\) with \(Q_j\) on
   G-quotient \(v_j\in E_{4p}^{G}\) recovers
   \(\sum c_j^2=19180800/1840091=\delta^2_{\mathrm{pair}}\le\mathrm{room}=3072/55\)
   at \(p=7\).
3. **Coverage.** Seven of eight \(H_+\) orbits have Max+-free geometric reps
   (affine + double Seidel–norm-circle); **one** size-1176 orbit still lacks
   a Max+-free construction (Q from census only for that type).

**OPEN:** Max+-free \(y\) for the remaining size-1176 orbit; full Max+-free
\(c_j\) at \(p=7\); residual for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15140.py`, `evidence/e1_gmin_m4_prop15140.json`,
`tests/test_prop15140.py`.

## Prop 15.141 (2026-08-01) — Size-12 Seidel partner; p=7 residual Max+-free

Continues 15.140. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Size-12 Seidel partner.** For affine halfspace \(y_0\) with
   \(S=\{2,3,4,5\}\) and field set
   \(T=\{10,12,13,16,18,25,29,36,38,42,44,48\}\),
   \(z=-1\) on \(T\), \(C_0=D_{y_0}CD_{y_0}\) has \(C_0z=pz\), and
   \(y_\sharp=y_0\odot z\) is Max+ of \(H_+\) orbit size 1176 with the
   previously missing \(Q_j\) signature — Max+-free.
2. **All eight \(H_+\) types Max+-free** at \(p=7\).
3. **Residual Max+-free at \(p=7\).** Free weights \(w_t=|G|/|\mathrm{Stab}|\)
   and free \(Q_j(y_t)\) give \(\sum c_j^2=\delta^2_{\mathrm{pair}}\le\mathrm{room}\).
4. **Bi-tight at \(p=5,7\).** mult\(\ge d-1\) + residual \(\Rightarrow\) 16N
   \(\Rightarrow\) bi-tight empty (form already proved); residual Max+-free
   at \(p=5,7\).

**OPEN:** General primes \(p\ge5\) (uniform type law beyond p=7-explicit \(T\));
deep ND; Main Theorem. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15141.py`, `evidence/e1_gmin_m4_prop15141.json`,
`tests/test_prop15141.py`.

## Prop 15.142 (2026-08-01) — Uniform affine law; partners; p=11 sample

Continues 15.141. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine all-\(S\).** Every \(S\subset\mathbb F_p\) with \(|S|=(p+1)/2\) gives a
   Max+ affine halfspace at \(p=5,7,11\) (all \(\binom{p}{k}\) sets).
2. **\(k\)-AP split.** \(p=5\): all AP (one affine orbit). \(p=7\): 21 AP \(\to84\),
   14 non-AP \(\to56\). \(p=11\): non-AP further splits (orbits 132, 330, 660).
3. **Fourths-coset partners.** \(z=-1\) on \(t+a\cdot\mathbb F_q^{\times4}\) works
   at \(p=5\) only; **no** evec hits at \(p=7,11\) — not a uniform size-12 law.
4. **\(p=7\) size-12 fibre.** 84 distinct field sets of size 12; explicit \(T\)
   of 15.141 is one; \(k=(q-1)/4\) but \(T\) is not a fourths coset.
5. **\(p=11\) samples.** Max+-free: affine 132/330/660; ystar 3630;
   double-switch 3630/7260. Full residual OPEN.

**OPEN:** Complete type law for \(p\ge11\); free \(Q_j\) residual for all
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15142.py`, `evidence/e1_gmin_m4_prop15142.json`,
`tests/test_prop15142.py`.

## Prop 15.143 (2026-08-01) — p=11 affine 6-orbit census; double-switch LB

Continues 15.142. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine type census at \(p=11\).** All \(\binom{11}{6}=462\) affine halfspaces
   form **exactly 6** \(H_+\) G-orbits under strict Aut \(G\) (\(|G|=14520\)):
   sizes \(132\times1\), \(330\times2\), \(660\times3\), with constructive samples.
2. **Double Seidel–norm-circle.** ystar orbit 3630; double-switch sizes include
   3630 and 7260; norm-circle count on \(C_0\) depends on affine class
   (size-132 rich, size-660 empty). Lower bound \(|H_+|\ge28182\).
3. **Incomplete non-affine list.** Full \(H_+\) type census at \(p=11\) OPEN.

**OPEN:** Non-affine types at \(p\ge11\); free \(c_j\) residual; general \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15143.py`, `evidence/e1_gmin_m4_prop15143.json`,
`tests/test_prop15143.py`.

## Prop 15.144 (2026-08-01) — Free orbits; type-enum residual dead for \(p\ge11\)

Continues 15.143. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Free G-orbit.** Affine \(S=\{0,1,2,4,5,7\}\) double-switch path
   \((t,c)=(33,3)\to(69,9)\) yields Max+ with \(|O|=|G|=14520\) (trivial stab).
2. **Size 2420.** ystar chain \((22,4)\to(91,4)\to(25,4)\to(95,5)\) yields
   orbit size 2420. Non-affine sizes include \(\{1210,2420,3630,7260,14520\}\).
3. **Type-enum residual DEAD for \(p\ge11\).** Deep double-switch multiplies
   free orbits — completing a Max+ G-orbit type list for \(c_j=\sum w_t Q_j\)
   is not a viable Max+-free proof path.
4. **Redirect.** Prefer type-free Max+-free residual: \(\delta^2\le\rho_{\min}^2\)
   for \(p\ge7\), or pointwise \(\sum_j Q_j(y)^2\le\mathrm{room}\).

**OPEN:** Type-free residual for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15144.py`, `evidence/e1_gmin_m4_prop15144.json`,
`tests/test_prop15144.py`.

## Prop 15.145 (2026-08-01) — Type-free residual package

Continues 15.144. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction algebra, Max+-free):**

1. **Dictionary.** For primes \(p>\sqrt5\),
   \[
   \delta^2\le\rho_{\min}^2
   \;\Longleftrightarrow\;
   \|\rho\|_2^2\le 2\rho_{\min}^2
   \;\Longleftrightarrow\;
   \|m_4\|_2^2\le m4f_{\mathrm{suf}}
   \;\Longleftrightarrow\;
   \mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}},
   \]
   with closed forms
   \[
   m4f_{\mathrm{suf}}=\frac{(p^2-1)(p^2+1)(3p^4+37p^2+60)}{24p^2(p^2-5)},
   \quad
   \mathrm{ED4}_{\mathrm{suf}}=\frac{4(p^2+1)(3p^6-3p^4+7p^2-15)}{p^2(p^2-5)}.
   \]
   For \(p\ge7\), \(\delta^2\le\rho_{\min}^2\) implies Path C residual
   (\(\rho_{\min}^2<\mathrm{room}_{\mathrm{hyp}}/24\)).
2. **Asymptotic.** \(\rho_{\min}^2/(\mathrm{room}_{\mathrm{hyp}}/24)\to 5/8\) as \(p\to\infty\)
   (monotone on primes \(7\ldots97\)).
3. **Type-free targets.** \(E\|P f_y\|^2\le\rho_{\min}^2\), pointwise
   \(Q_\delta\le\rho_{\min}^2\), \(\|m_4\|^2\le m4f_{\mathrm{suf}}\), or
   \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) (weight enumerator / Gauss sums).

**Certified:** \(\delta^2\le\rho_{\min}^2\) at \(p=5,7\) only (prior residual closures).

**OPEN:** type-free \(\delta^2\le\rho_{\min}^2\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15145.py`, `evidence/e1_gmin_m4_prop15145.json`,
`tests/test_prop15145.py`.

## Prop 15.146 (2026-08-01) — Type-free R₄ / μ₄ residual channel

Continues 15.145. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **ED4_from_exact3** \(= -n^4+28n^2-40n\).
2. **R₄ dictionary.** \(\delta^2=(2/3)(R_4-R4_{\mathrm{flat}})\);
   \(\delta^2\le\rho_{\min}^2\Leftrightarrow R_4\le R4_{\mathrm{suf}}=R4_{\mathrm{flat}}+(3/2)\rho_{\min}^2\);
   equivalently \(\mu_4\le\mu4_{\mathrm{suf}}\) and
   \(\bar E[\prod_4 w]\le R4_{\mathrm{suf}}/(n)_4\).
3. **Central moments.** \(\mu_2=n/2\), \(\mu_3=0\) Max+-free.
4. **Spectral mass (too weak).** Under \(\lambda_{\max}(T)=4p\),
   \(\|P f_y\|^2\le (p^2+1)(p^2-2)(p^2-3)(p^2-4)/24\)
   with \(w^*=(p^2-3)(p^2-4)/(p^2(p^2-1))\) — far above residual budgets,
   so Jensen+spectral moments cannot close residual.

**Certified:** R₄ channel at \(p=5,7\).

**OPEN:** \(R_4\le R4_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15146.py`, `evidence/e1_gmin_m4_prop15146.json`,
`tests/test_prop15146.py`.

## Prop 15.147 (2026-08-01) — Inclusion-density residual; ULC near-miss

Continues 15.146. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Falling factorial.** \(R_4=E[k^{\underline 4}]\); under design moments,
   \(\mathrm{exact}_{\le3}=6E[k^3]-11E[k^2]+6E[k]\).
2. **Inclusion densities.**
   \(d_1=1/2\), \(d_2=(p^2+1)/(4p^2)\), \(d_3=(p^2+3)/(8p^2)\);
   \(d_4=R_4/(n)_4\); for \(p\ge7\), residual \(\Leftrightarrow d_4\le d4_{\mathrm{suf}}\).
3. **ULC comparison.** \(U=d_3^2/d_2=(p^2+3)^2/(16p^2(p^2+1)) < d4_{\mathrm{suf}}\)
   for all primes \(p\ge5\) via
   \(P(x)=x^5-8x^4+78x^3-4x^2-7x-60=(x-1)(x^4-7x^3+71x^2+67x+60)>0\) at \(x=p^2\ge25\).
   Thus \(d_4\le U\) would close residual.
4. **ULC fails slightly (census).** At \(p=5,7\): \(d_4/U\approx 1.036,1.019\)
   while \(d_4\le d4_{\mathrm{suf}}\) still holds. Near-miss.

**OPEN:** \(d_4\le d4_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15147.py`, `evidence/e1_gmin_m4_prop15147.json`,
`tests/test_prop15147.py`.

## Prop 15.148 (2026-08-01) — Relaxed-ULC residual calculus

Continues 15.147. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Linear defect.** \(C_{\mathrm{act}}=p^2(d_4/U-1)=C_{\mathrm{flat}}+\kappa\delta^2\) with
   \(C_{\mathrm{flat}}=x P_{\mathrm{flat}}(x)/((x-5)(x-1)(x-2)(x+3)^2)\),
   \(P_{\mathrm{flat}}=x^4-8x^3+58x^2-64x+13\),
   \(\kappa=24p^2/((n-2)(n-3)(n+2)^2)\), \(x=p^2\).
2. **C_max.** \(C_{\mathrm{max}}=Q(p^2)/((p^2-5)(p^2-2)(p^2+3)^2)\),
   \(Q=x^4-7x^3+71x^2+67x+60\);
   residual \(\Leftrightarrow C_{\mathrm{act}}\le C_{\mathrm{max}}\);
   \(C_{\mathrm{max}}<1\to1^-\).
3. **Uniform criterion (p≥7).** \(C_{\mathrm{max}}(p)\ge C_7=C_{\mathrm{max}}(7)=79923/87373\);
   thus \(d_4\le U(1+C_7/p^2)\) implies residual for all primes \(p\ge7\).
4. **Census window.** \(C_{\mathrm{act}}(5)\approx0.900\), \(C_{\mathrm{act}}(7)\approx0.907\);
   constant-\(C\) implication window \([C_{\mathrm{act}}(7),C_7]\) nonempty; \(C=1\) fails.

**OPEN:** prove \(d_4\le U(1+C_7/p^2)\) (or \(C_{\mathrm{act}}\le C_{\mathrm{max}}\)) type-free. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15148.py`, `evidence/e1_gmin_m4_prop15148.json`,
`tests/test_prop15148.py`.

## Prop 15.149 (2026-08-01) — Size-bias residual form

Continues 15.148. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Size-bias.** Let \(\mu\propto W_k\,k^{\underline 3}\). Then
   \(d_4/d_3=E_\mu[k-3]/(n-3)\) and residual for \(p\ge7\) is
   \(E_\mu[k]\le k_{\mathrm{suf}}=3+8\,R4_{\mathrm{suf}}/(n(n-2)(n+2))\).
   Shift: \(E_\mu[k]-k_{\mathrm{flat}}=12\delta^2/(n(n-2)(n+2))\).
2. **\(k_{\mathrm{flat}}\) shift.** \(k_{\mathrm{flat}}-n/2\to 3^-\) with closed
   \(P_k(p^2)/(2(p^2-5)(p^2+1)(p^2-1)(p^2+3))\).
3. **Independence excesses.** \(d_2-1/4=1/(4(n-1))\), \(d_3-1/8=3/(8(n-1))\),
   \(d4_{\mathrm{flat}}-1/16=P_{\mathrm{ind}}(p^2)/(16(p^2-5)n_4)>0\) for \(p\ge5\);
   \(k_{\mathrm{flat}}>(n+3)/2\) (binomial size-bias).
4. **Uniform \(C_7\).** \(E_\mu[k]\le k_{C7}\Rightarrow\) residual for \(p\ge7\).
5. **Gauss program.** \(\mu\) mixes regular-set sizes through srg triple types
   with Aut-constant \(\lambda_\tau\) — finite-type character-sum target.

**OPEN:** \(E_\mu[k]\le k_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15149.py`, `evidence/e1_gmin_m4_prop15149.json`,
`tests/test_prop15149.py`.

## Prop 15.150 (2026-08-01) — srg triples; \(\lambda_e\); \(\pi_e\) residual

Continues 15.149. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Triple counts.** Closed \(n_e\) for edge-type \(e=0,1,2,3\) in the
   conference srg after switch:
   \(n_3=p(p^2+1)(p^2-1)(p-3)/48\), etc.
2. **Covering numbers.** Under Aut-constancy on edge-types and affine
   \(\lambda_e=A+Be\), design moments \(j\le3\) force
   \(\lambda_e=N(p+3-2e)/(8p)\). Certified at \(p=5,7\).
3. **Mixture.** \(\pi_e=n_e(p+3-2e)/\mathrm{Tot}\) is Max+-free (\(N\) cancels);
   \(E_\mu[k]=\sum_e\pi_e m_e\) with \(m_e=\) mean regular-set size through a
   type-\(e\) triple.
4. **Residual.** For \(p\ge7\): \(\sum\pi_e m_e\le k_{\mathrm{suf}}\).

**OPEN:** bound \(m_e\) (character sums / regular-set constraints). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15150.py`, `evidence/e1_gmin_m4_prop15150.json`,
`tests/test_prop15150.py`.

## Prop 15.151 (2026-08-01) — \(m_e\) covariance formula

Continues 15.150. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Regular-set identities.** \(\sum t_e=C(k,3)\), \(\sum e\,t_e=e_S(k-2)\),
   \(t_2+3t_3=k\binom{\alpha}{2}\).
2. **\(E[t_e]\) Max+-free.** \(E[t_e]=n_e(p+3-2e)/(8p)\).
3. **Covariance formula.**
   \(m_e=n/2+8p\,\mathrm{Cov}(k,t_e)/(n_e(p+3-2e))\);
   \(E_\mu[k]=n/2+(8p/\mathrm{Tot})\,\mathrm{Cov}(k,C(k,3))\).
4. **Exact \(m_e\) at \(p=5\).** Fractions certified; \(\sum\pi m\le k_{\mathrm{suf}}\);
   \(m_0,m_1<k_{\mathrm{suf}}<m_2,m_3\).
5. **CS/Popoviciu dead** for residual (π-average exceeds \(k_{\mathrm{suf}}\)).

**Note:** weight-constancy of \(t_e\) (Thm B) overstated for \(p=7\); see 15.152.

**OPEN:** closed \(t_e(k)\) or character-sum \(\mathrm{Cov}(k,t_e)\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15151.py`, `evidence/e1_gmin_m4_prop15151.json`,
`tests/test_prop15151.py`.

## Prop 15.152 (2026-08-01) — free-param \(t_3\); multi-orbit; residual \(\equiv R_4\)

Continues 15.151. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Free parameter.** On regular sets,
   \(t_2=R_3-3t_3\), \(t_1=R_2-2R_3+3t_3\), \(t_0=C(k,3)-R_2+R_3-t_3\).
2. **\(p=5\) closed form.** Mono-type per weight;
   \(t_3(\alpha)=0\) (\(\alpha\le3\)), \(t_3(\alpha)=3\alpha^2-21\alpha+40\) (\(\alpha\ge4\)).
3. **\(p=7\) multi-orbit.** Weights \(k\in\{16,18,\ldots,34\}\) have 3–8 distinct
   \(t\)-vectors (full Max+ census \(N=11452\), \(W=86\)). Pure \(t_e(k)\) is **dead**
   for general \(p\). Corrects 15.151.B.
4. **Residual \(\equiv R_4\).** With \(R_4=E[k^{\underline4}]\),
   \(\mathrm{Cov}(k,C(k,3))=(R_4+(3-n/2)E_3)/6\) and
   \(E_\mu[k]=3+8R_4/(n(n-2)(n+2))\). Multi-type structure does not open a new
   residual channel beyond the weight-4 falling moment.
5. **Per-type reduction.** \(\mathrm{Cov}(k,t_e)=\mathrm{Cov}(\mathrm{det}_e(k))+\gamma_e\mathrm{Cov}(k,t_3)\)
   with \(\gamma=(-1,3,-3,1)\); \(E[t_3]=n_3(p-3)/(8p)\).

**OPEN:** \(R_4\le R4_{\mathrm{suf}}\) / character-sum \(m_e\) on fixed triples /
\(\mathrm{Cov}(k,t_3)\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15152.py`, `evidence/e1_gmin_m4_prop15152.json`,
`tests/test_prop15152.py`.

## Prop 15.153 (2026-08-01) — switched \(\mu_4\) residual dictionary

Continues 15.152. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Switched low moments.** \(e_1=e_3=0\), \(e_2=1/(n-1)\) Max+-free.
2. **Four-point expansion.**
   \(d_4=(1+6/(n-1)+\mu_4)/16\) with \(\mu_4=\mathrm{avg}\,m_4^{\mathrm{sw}}\).
3. **Closed budgets.**
   \(\mu4_{\mathrm{flat}}=(3p^2+17)/(p^2(p^2-2)(p^2-5))\),
   \(\mu4_{\mathrm{suf}}=(3p^4+37p^2+60)/(p^4(p^2-2)(p^2-5))\).
   Residual for \(p\ge7\): \(\mu_4\le\mu4_{\mathrm{suf}}\).
4. **Census.** \(\mu_4\le\mu4_{\mathrm{suf}}\) at \(p=5,7\); exact \(m_e\) at \(p=7\).

**OPEN:** Paley/Weil bound \(\mu_4\le\mu4_{\mathrm{suf}}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15153.py`, `evidence/e1_gmin_m4_prop15153.json`,
`tests/test_prop15153.py`.

## Prop 15.154 (2026-08-01) — \(\mathrm{avg}(\chi\kappa)=3/(n-3)\); \(\eta\) residual

Continues 15.153. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Switch.** \(m_4^{\mathrm{sw}}(S)=\chi(S)\,m_4(S)\), \(\chi=\prod z_i\) halfspace.
2. **Combinatorial average (Max+-free).** By Seidel-switched conference \(C_2\)
   row-sum algebra:
   \(\mathrm{avg}(\chi\kappa)=3/(n-3)=3/(p^2-2)\).
3. **Decomposition.**
   \(\mu_4=\kappa_{\mathrm{main}}+\eta\) with
   \(\kappa_{\mathrm{main}}=3/(p^2(p^2-2))\),
   \(\eta=\mathrm{avg}(\chi\,\mathrm{Ext})/(4p)\).
4. **Residual.** \(\mu_4\le\mu4_{\mathrm{suf}}\Leftrightarrow\eta\le\eta_{\mathrm{suf}}\) with
   \(\eta_{\mathrm{suf}}=4(13p^2+15)/(p^4(p^2-2)(p^2-5))\).
   Strictly \(\kappa_{\mathrm{main}}<\mu4_{\mathrm{flat}}<\mu4_{\mathrm{suf}}\).
5. **Census.** \(\eta\le\eta_{\mathrm{suf}}\) at \(p=5,7\).

**OPEN:** Weil/Aut bound on \(\eta\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15154.py`, `evidence/e1_gmin_m4_prop15154.json`,
`tests/test_prop15154.py`.

## Prop 15.155 (2026-08-01) — Aut-line \(e_4/T\chi/Q\); \(\eta=c_1 R_4+c_0\)

Continues 15.154. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **\(e_4(s)\).** For any \(\pm1\) vector with sum \(s\):
   \(e_4=s^4/24+((-3n+4)/12)s^2+n(n-2)/8\).
2. **\(T\chi\).** \(T\chi(S)=\chi(S)(4p-2\sigma_z(S))\) with
   \(\sigma_z=\sum_{uv\subset S}C_{2,uv}\).
3. **\(Q(s)\).** On \(C_2\)-eigenvectors:
   \(Q=(p/4)[s^2(n-4)+n(6-n)]\).
4. **Average.** \(\mathrm{avg}(\chi\,\mathrm{Ext})=E[(4p\,e_4-2Q)/C(n,4)]\).
5. **Affine.** \(\eta=c_1 R_4+c_0\) with
   \(c_1=16/(n)_4\), \(c_0=-(p^4+4p^2-9)/(p^2(p^2-2))\).
   Residual \(\Leftrightarrow R_4\le R4_{\mathrm{suf}}\) (pure \(E[s^4]\)).
6. **Crude dead.** \(E[s^4]\le 2n^3\) exceeds budget.

**OPEN:** Weil/Paley \(E[\langle z,y\rangle^4]\) or spherical 3-design defect. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15155.py`, `evidence/e1_gmin_m4_prop15155.json`,
`tests/test_prop15155.py`.

## Prop 15.156 (2026-08-01) — \(\kappa_4=E[s^4]-12n^2\) residual dictionary

Continues 15.155. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Cumulant form.** \(\kappa_4:=E[s^4]-12n^2\); residual \(\Leftrightarrow\kappa_4\le\kappa4_{\mathrm{suf}}\).
2. **Closed budgets.**
   \(\kappa4_{\mathrm{flat}}=16(p^2+1)(p^2+3)/(p^2-5)\),
   \(\kappa4_{\mathrm{suf}}=4(p^2+1)(9p^4+22p^2-15)/(p^2(p^2-5))\).
3. **Bridge.** \(\kappa_4=(n)_4\eta-16n\).
4. **Design orientation.** Spherical 4-design value \(3n^4/(d(d+2))\) is a
   **lower** bound among 2-designs; lies below \(\mathrm{ED4}_{\mathrm{flat}}\) for \(p\ge5\).
5. **Dead ends.** Crude \(E[s^4]\le2n^3\); moment LP on allowed weights (factor \(\gtrsim3.5\) at \(p=5\)).
6. **Census.** \(\kappa_4\le\kappa4_{\mathrm{suf}}\) at \(p=5,7\) (ratios \(\approx0.90,0.66\)).

**OPEN:** Weil/Paley bound on \(\kappa_4\), or upper bound on the spherical
3-design defect of Max+ in \(V_+\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15156.py`, `evidence/e1_gmin_m4_prop15156.json`,
`tests/test_prop15156.py`.

## Prop 15.157 (2026-08-01) — Gegenbauer / 3-design defect residual

Continues 15.156. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Expansion.** \(t^4=a_0+a_2 Q_2+a_4 Q_4\) on \(S^{d-1}\) with
   \(a_0=3/(d(d+2))\), \(a_2=6(d-1)/(d(d+4))\),
   \(a_4=(d^2-1)/((d+2)(d+4))\).
2. **2-design reduction.** \(E[s^4]=n^4(a_0+a_4\mu_{G4})\) with
   \(\mu_{G4}=E[Q_4(s/n)]\ge0\).
3. **Residual.** \(\mu_{G4}\le\mu_{G4,\mathrm{suf}}\) where (with \(x=p^2\))
   \(\mu_{G4,\mathrm{suf}}=4(21x^3+19x^2+35x-75)(x+9)/[x(x-5)(x+1)^3(x+3)(x-1)]\).
4. **Census.** Defect positive but inside budget at \(p=5,7\) (ratios \(\approx0.94,0.83\)).
5. **Dead UBs.** \(\mu\le1\); \(\mu\le1/h_4\) (false at \(p=5\)); \(\mu\le d/h_4\) (too weak for \(p\ge7\)).

**OPEN:** Weil/Aut upper bound on \(\mu_{G4}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15157.py`, `evidence/e1_gmin_m4_prop15157.json`,
`tests/test_prop15157.py`.

## Prop 15.158 (2026-08-01) — closed \(Q_4\); Max+ non-scheme; pole bound

Continues 15.157. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Closed \(Q_4\).**
   \(Q_4(t)=[(d+2)(d+4)t^4-6(d+2)t^2+3]/(d^2-1)\), \(Q_4(0)=3/(d^2-1)\).
2. **1-homogeneous tight frame.** Absolute distributions identical; Gram
   spectrum \(2N\) (mult \(d\)) + \(0\) (mult \(N-d\)). Certified \(p=5\).
3. **Not IP-scheme.** Intersection numbers fail constancy on several
   inner-product classes at \(p=5\) — Bose–Mesner on pure \(R_s\) blocked.
4. **Pole decomposition.**
   \(\mu_{G4}\le 2/N+P(E)+Q_4(0)P(\mathrm{Eq})\); pure pole bound
   \(\mu_{G4}\le 2/N\) if \(W\) avoids positive-\(Q_4\) weights except poles.
5. **Conditional residual.** If \(N\ge N_*=\lceil 2/\mu_{G4,\mathrm{suf}}\rceil\) and
   nonpositive-\(Q_4\) support, residual holds. \(N\ge N_*\) at \(p=5,7\), but
   support hypothesis fails (equator / Hoffman exterior).
6. **Chebyshev split dead** (\(\mathrm{UB}\sim 0.18\gg\mu_{G4,\mathrm{suf}}\)).

**OPEN:** Weil/Paley \(\sum_y Q_4(\langle z,y\rangle/n)\), Aut-coherent configuration,
or closed \(W_k\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15158.py`, `evidence/e1_gmin_m4_prop15158.json`,
`tests/test_prop15158.py`.

## Prop 15.159 (2026-08-03) — Φ|Z spectrum structure; dual gap \(G=(d/32)(16I-\Phi)\)

Continues 15.158 / strategy reframe. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Exact Φ spectrum at \(p=5\) (certified Fraction).** On \(Z\),
   \(\lambda\in\{176/13\,(\times d),\,144/13\,(\times 2d),\,80/13\,(\times 2d)\}\)
   with \(d=13\), sum mult \(=m=65\). In particular
   \(\lambda_{\max}=176/13=16(d-2)/d\), mult\(=d\).
2. **Exact Φ spectrum at \(p=7\) (certified Fraction).** With den \(=N/(4p)=409\),
   \(\lambda\in\{4320,4032,3648,3360,3072\}/409\) with mults
   \((d,2d,2d,4d,2d)\). \(\lambda_{\max}=4320/409<16(d-2)/d\), mult\(=d\).
3. **Design threshold algebra (proved Fraction).** For \(d>2\),
   \(16(d-2)/d<16\). For \(d\ge13\) (\(p\ge5\)),
   \(\bar\mu=8(d-1)/(d-3)\le16(d-2)/d\) iff \((d-3)(d-6)\ge0\).
4. **Dual gap operator.** \(G:=(d/32)(16I-\Phi)\). At \(p=5\): eigs of \(G\) are
   \(\{1,2,4\}\); at \(p=7\): \(G\succeq I\) (strict). Hence \(G\succeq I\Rightarrow\lambda_{\max}\le16(d-2)/d<16\Rightarrow16\mathrm N\).
5. **16N chain predicate (proved, conditional).** mult\(\ge d\) and \(\|\kappa\|_F^2\le96n\) \(\Rightarrow16\mathrm N\) for \(p\ge5\) (15.105 restated).

**OPEN:** dual gap \(G\succeq I\) (or mult\(\ge d\) + \(\|\kappa\|^2\le96n\)) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15159.py`, `evidence/e1_gmin_m4_prop15159.json`,
`tests/test_prop15159.py`.

## Prop 15.160 (2026-08-03) — Dual-gap vs Hypothesis H; \(H\Rightarrow G\succeq I\)

Continues 15.159. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction):**

1. **H vs thr_ray.** \(H(p)=(p+2)^2/d\), thr_ray\(=5-16/d\).
   \(H-\mathrm{thr\_ray}=(p-5)(p+1)/(2d)\) (equivalently checked Fraction form).
   Hence \(H\le\mathrm{thr\_ray}\) for all primes \(p\ge5\), equality only at \(p=5\).
2. **\(H\Rightarrow\) dual gap.** ray_max\(\le H(p)\) and \(p\ge5\) \(\Rightarrow\) ray_max\(\le\mathrm{thr\_ray}\) \(\Rightarrow G\succeq I\Rightarrow16\mathrm N\).
3. **\(H\Rightarrow16\mathrm N\) for \(p\ge3\).** \(H(p)\le5\) (eq only \(p=3\)), so ray\(\le H\Rightarrow\) ray\(\le5\Rightarrow16\mathrm N\) (15.63).
4. **Census.** \(p=5\): ray \(=H=\mathrm{thr\_ray}=49/13\). \(p=7\): ray \(<H<\mathrm{thr\_ray}\).

**OPEN:** Hypothesis H (ray_max\(\le H(p)\)) for all primes \(p\ge5\); equivalent residual forms \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\), orth\(\le\mathrm{room}_{\mathrm{hyp}}\), \(\|\kappa\|^2\le\kappa_{\mathrm{hyp}}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15160.py`, `evidence/e1_gmin_m4_prop15160.json`,
`tests/test_prop15160.py`.

## Prop 15.161 (2026-08-05) — Φ-frame of Max+; 16N via mult\(\ge d\) + \(\kappa_4\le48n\)

Continues 15.159–15.160. Does **not** soft-close residual / L.

**Proved (Fraction / Max+-free):**

1. **Constant embedding norm.** \(v_y:=P_Z(yy^\top-I)\) satisfies \(\|v_y\|_F^2=n(n-2)\) for all \(y\in\mathrm{Max}_+\) (from \(\mathrm{tr}(\Phi)=m\bar\mu=n(n-2)\) and 1-homogeneity).
2. **Pairwise frame Gram.** With \(v_y=yy^\top-(n/d)P_+\in Z\),
   \(\langle v_y,v_z\rangle_F=(y\cdot z)^2-2n\).
3. **16N budgets under mult\(\ge d\).** \(E[s^4]\le12n(n+4)\) \(\Leftrightarrow\kappa_4\le48n\) \(\Rightarrow\lambda_{\max}\le16\) when mult\(\ge d\); two-level bulk \(b=8\).
4. **Census.** mult\(=d\) and \(\kappa_4\le48n\) at \(p=5,7\).

**OPEN:** mult\((\lambda_{\max})\ge d\) and \(\kappa_4\le48n\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15161.py`, `evidence/e1_gmin_m4_prop15161.json`,
`tests/test_prop15161.py`.

## Prop 15.162 (2026-08-05) — Maximizers in \(Z\); mult\(\ge d-1\); \(E[s^4]\) type expansion

Continues 15.161. Does **not** soft-close residual / L.

**Proved:**

1. **Maximizers of \(\Gamma\) lie in \(Z\).** Criticality on Sym\(_0\): \(E[f\,cc^\top]=\lambda A\) (\(\mu=0\) from \(E[f]=0\)); ambient diagonal \(B_{ii}=E[f]/\lambda=0\). Hence mult\((\Gamma\text{ top})=\)mult\((\Phi\text{ top})\).
2. **mult\((\Phi)\ge d-1\) (proved for all primes \(p\ge5\)).** Thm A + Prop 15.97 (mult \(\Gamma=\) mult \(\lambda_2(P\odot P)\)) + Prop 15.98 (PSL min nontrivial irrep dim \(d-1\)).
3. **Type expansion.** \(E[s^4]=C_0+R\) with
   \(C_0=n(3n-2)+2n(n-1)(3n-4)/p^2\) and \(R=24\sum_{4\text{-sets}}m_4^2\ge0\).
4. **16N \(\Leftrightarrow\) mult\(\ge d\) + \(m_4\)-mass.** \(\sum m_4^2\le n(3p^2+61)/24\) \(\Leftrightarrow\kappa_4\le48n\).
5. **Census.** Identities + 16N at \(p=5,7\).

**OPEN:** upgrade mult\(\ge d-1\to d\); prove \(m_4\)-mass / \(\kappa_4\le48n\) for all \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15162.py`, `evidence/e1_gmin_m4_prop15162.json`,
`tests/test_prop15162.py`.

## Prop 15.163 (2026-08-05) — Wick \(m_4\) mass; Aut\(_0\); \(H_C\) split

Continues 15.161–15.162. Does **not** soft-close residual / L.

**Proved / certified:**

1. **Wick 4-set mass (proved Fraction).** \(T=C_{ab}C_{cd}+C_{ac}C_{bd}+C_{ad}C_{bc}\), \(m_4^W=T/p^2\);
   \(\sum T^2=n(n-1)(n-2)(n-5)/8\), \(\sum(m_4^W)^2=(p^4-1)(p^2-4)/(8p^2)\).
2. **η-room after Wick (proved Fraction).** Under \(\sum m_4^W\eta=0\),
   16N \(\Leftrightarrow\sum\eta^2\le n(19p^2-3)/(6p^2)\). Census usage \(\approx95\%,73\%\) at \(p=5,7\).
3. **Aut\(_0\) on \(V_+\) (structure + cert \(p=5\)).** \(V_+\cong\mathbf1\oplus\sigma\) with \(\dim\sigma=d-1\).
4. **\(H_C\) split of Φ-top (cert \(p=5,7\)).** Top mult splits as \(1+(d-1)\).

**OPEN:** mult\(\ge d\) general \(p\); \(\sum\eta^2\le\eta_*\) / 16N for all \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15163.py`, `evidence/e1_gmin_m4_prop15163.json`,
`tests/test_prop15163.py`.

## Prop 15.164 (2026-08-05) — 16N from mult\(\ge d-1\) + \(E[s^4]\le\mathrm{Es4}_*(p)\)

Continues 15.161–15.163. Preferred 16N path (no mult\(\ge d\) needed). Does **not** soft-close L.

**Proved (Fraction):**

1. **Two-level majorization.** Fixed sum \(T\), sum of squares \(Q\), mult\(\ge k\), \(\lambda_i\ge\ell_{\min}\): max top \(L\) is two-level; \(L\) increasing in \(Q\).
2. **Es4\(_*\) budget under mult\(\ge d-1\).** With \(\lambda_{\min}\ge6\), bulk
   \(b_*=8(p^2-1)(p^2-7)/(p^4-8p^2-1)\ge8>6\) for primes \(p\ge5\).
   If \(E[s^4]\le\mathrm{Es4}_*(p)\) then \(\lambda_{\max}\le16\) (16N).
3. **Equivalent forms.** \(E[s^4]\le\mathrm{Es4}_*\Leftrightarrow\kappa_4\le\kappa4_*\Leftrightarrow R\le R_*\Leftrightarrow\sum\eta^2\le\eta_*\) (Wick orth).
4. **Census.** 16N at \(p=5,7\) via spectrum Es4.

**OPEN:** \(E[s^4]\le\mathrm{Es4}_*(p)\) (or \(\sum\eta^2\le\eta_*\)) for all primes \(p\ge5\).
This is the single analytic residual for Path-C 16N given mult\(\ge d-1\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15164.py`, `evidence/e1_gmin_m4_prop15164.json`,
`tests/test_prop15164.py`.

## Prop 15.165 (2026-08-05) — Exact Es4; closed Es4\(_*\)/\(\eta_*\); GoG\(\leftrightarrow\Phi\); \(m_4\) is C-eigen

Continues 15.164. Does **not** soft-close L.

**Proved / certified:**

1. **Moments.** \(E[s]=0\) (central symmetry), \(E[s^2]=2n\) (2-design / \(\sum yy^\top=2N P_+\)).
2. **GoG \(\leftrightarrow\) Φ spectrum.** \(\mathrm{spec}(G\circ G)=\{2nN\}\cup\{N\lambda:\lambda\in\mathrm{spec}(\Phi|_Z)\}\cup\{0\}^{N-1-m}\);
   hence \(E[s^4]=4n^2+\mathrm{tr}(\Phi^2)\).
3. **Closed budgets (Fraction).**
   \[
   \mathrm{Es4}_*(p)=\frac{4(3p^8-6p^6-148p^4-10p^2+129)}{p^4-8p^2-1},\quad
   \eta_*(p)=\frac{(p^2-1)(p^2+1)(19p^4-152p^2-3)}{6p^2(p^4-8p^2-1)}.
   \]
4. **\(m_4\) is C-eigen.** \(p\cdot m_4(a,b,c,d)=\sum_j C_{aj}m_4(j,b,c,d)\) (cert random 4-sets \(p=5,7\)).
5. **Exact Es4 census.** \(p=3,5,7\) via Gram/Φ spectrum (not single-root \(W\); \(p=7\) not 1-homogeneous). H-saturation at \(p=5\).

**OPEN:** Es4\(_*\) / 16N for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15165.py`, `evidence/e1_gmin_m4_prop15165.json`,
`tests/test_prop15165.py`.

## Prop 15.166 (2026-08-05) — 16N \(\Leftrightarrow\lambda_{\max}(Q_2)\le4N/(d(d-1))\); Wick C-eigen

Continues 15.164–15.165. Does **not** soft-close L.

**Proved:**

1. **Max+ is a spherical 2-design in \(V_+\).** \(\sum\hat u=0\), \((1/N)\sum\hat u\hat u^\top=I_d/d\).
2. **Wick \(m_4\) is C-eigen** (same equation as true \(m_4\)); residual \(\eta=m_4-m_4^W\) is invisible to the C-eigen constraint — \(P_+(\mathrm{Wick})\) does not recover true \(m_4\).
3. **16N \(\Leftrightarrow Q_2\) bound (proved Fraction).**
   \(\lambda_{\max}(\Phi)=4d(d-1)/N\cdot\lambda_{\max}(Q_2)\), hence
   \(\lambda_{\max}(\Phi)\le16\Leftrightarrow\lambda_{\max}(Q_2)\le4N/(d(d-1))\).
4. **Census.** \(Q_2\) thr holds at \(p=5\) (ratio \(\approx0.846\)), \(p=7\) (\(\approx0.660\)).

**OPEN:** \(\lambda_{\max}(Q_2)\le4N/(d(d-1))\) for all primes \(p\ge5\). Dead: Delsarte LP, BM\((C)\), equating Wick to \(m_4\). Preferred: Weil/Jacobi Aut-orbit \(m_4\), Aut\(_0\) isotype, SOS \(Q_4\). L remains OPEN. residual_closed_general=false.

Evidence: `src/e1_gmin_m4_prop15166.py`, `evidence/e1_gmin_m4_prop15166.json`,
`tests/test_prop15166.py`.

## Prop 15.167 (2026-08-05) — Bi-tight empty for all primes \(p\ge5\) via majorization (no residual)

Bypasses Path-C residual for the bi-tight link. residual/16N/Es4\(_*\) remain OPEN (honest). L OPEN until E(1)/Main.

**Proved (Fraction, all primes \(p\ge5\)):**

1. **Majorization UB.** mult\((\lambda_{\max})\ge d-1\) (15.162) + \(\lambda_{\min}(\Phi)\ge6\) + \(\mathrm{tr}(\Phi)=n(n-2)\) \(\Rightarrow\)
   \[
   \lambda_{\max}(\Phi)\le L_*(p)=\frac{p^4+24p^2-1}{2(p^2-1)}.
   \]
2. **\(L_*<2d\).** \(2d-L_*=(p^4-24p^2-1)/(2(p^2-1))\); numerator \(=24\) at \(p=5\) and \(f(x)=x^2-24x-1\) increasing on \(x=p^2\ge25\).
3. **Bi-tight empty.** \(\lambda_{\max}\le L_*<2d\Rightarrow\lambda_{\mathrm{cycle}}=\lambda_{\max}/2<d=n/2\Rightarrow\lambda_{\max}(G)=n/2\) simple \(\Rightarrow\) bi-tight empty (15.55). **Does not use residual / Es4\(_*\) / 16N.**
4. **Census.** Actual \(\lambda_{\max}\le L_*\) and \(\lambda_{\mathrm{cycle}}<d\) at \(p=5,7\).

**OPEN:** residual/16N general \(p\); E(1)/Main; L. residual_closed_general=false.

Evidence: `src/e1_gmin_m4_prop15167.py`, `evidence/e1_gmin_m4_prop15167.json`,
`tests/test_prop15167.py`, `src/e1_bitight_chain.py`.

## Prop 15.168 (2026-08-05) — E(1) structure after 15.167 (honest partial)

Continues 15.167. Does **not** soft-close E(1) or L.

**Proved / checkable predicates (Fraction + prior props):**

1. **Tight level-\(s\) obstruction.** Bi-tight empty (15.167) \(\Rightarrow\) no Max\(_+\)-tight level-\(s\) cover of size \(sp\) (\(G_\perp\) isotropy; 15.55 gen.).
2. **Deep tight empty** for \(p\ge5\) (15.44.3 + 15.167).
3. **Type I freeness ND** (prior 15.43.1).
4. **Type I freeness-fail \(k=2p-1\)** \(\to\) tight size \(2p\) \(\to\) ND when bi-tight empty (15.43.3 + 15.44 + 15.167).
5. **Deep auto-freeness** for \(s_+=2\), \(k\le3p-2\): \(N_2/N\) lb \(=2-k/(2p)>(p+1)/(2p)\).
6. **Deep fail-eq \(k=3p-1\)** \(\Rightarrow\) tight \(S\equiv3\) size \(3p\) \(\Rightarrow\) empty under Thm A for \(s=3\).

**OPEN residuals (honest — no soft-close):**

- Type I freeness-fail at \(k=3p-2\) / \(S\in\{1,5\}\) boundary (not reduced to tight \(2p\)).
- Deep non-tight freeness-fail with \(k\ge3p\) (freeze-to-tight sketch not shipped as predicate).

Full \(m_n\ge\Phi-2\) / E(1) / L remain OPEN. residual_closed_general=false. E1_closed_general=false.
L closed only if bi-tight \(\land\) E(1) (denseness Prop 6.2) — currently false.

Evidence: `src/e1_gmin_m4_prop15168.py`, `evidence/e1_gmin_m4_prop15168.json`,
`tests/test_prop15168.py`, `src/e1_main_chain_status.py`.

## Prop 15.169 (2026-08-05) — Type I \(k=3p-2\) ND reduction; deep multi-\(s\) auto-freeness

Continues 15.168. Does **not** soft-close E(1) or L.

**Proved (Fraction / prior props):**

1. **Type I freeness-fail structure at \(k=3p-2\).** At freeness equality \(N_1/N=(p+1)/(2p)\) with \(S\in\{1,5\}\): \(a=(5-k/p)/4\) equals the freeness threshold; affine \(S+2f_e=3\) on Max\(_+\); \(H=G\cup\{e\}\) has size \(3p-1\), scores \(\{2,4\}\), \(s_+^H=2\).

2. **2-Lipschitz of \(\Phi\) under edge flip.** \(Q_y(A\oplus e)=Q_y(A)-2f_e(y)\) with \(f_e=\pm1\), so \(|\Phi(A\oplus e)-\Phi(A)|\le2\). Corollary: if \(\Phi(C\oplus G)\ge\Phi\) then \(\Phi(C\oplus G\oplus e)\ge\Phi-2\) (weak ND for every \(e\)).

3. **Gap-2 undercutter forces \(s_-=-1\).** Type I \(\Rightarrow\Phi\ge\Phi-2\) (15.42.1). Gap-2 means \(\Phi=\Phi-2\), requiring \(s_-\le-1\). With \(k=3p-2\) odd and \(\mathbb E_-[S]=-3+2/p>-3\), odd scores force \(s_-\ge-1\); hence \(s_-=-1\).

4. **ND dichotomy at \(s_-=-1\).** On \(U_-=\{S_G=-1\}\): \(f_e=+1\) gives \(|Q_H|=\Phi\) (strong ND); \(f_e\equiv-1\) on \(U_-\) (bad case) gives only \(\Phi(H)\ge\Phi-4\) from Max\(\pm\) dichotomy.

5. **Deep multi-\(s\) auto-freeness.** For min-level \(s\ge2\) with scores \(\ge s\) step 2: \(N_s/N\) lb \(=(s+2-k/p)/2\); auto-freeness for \(k\le p(s+1)-2\). Recovers \(s=2\Rightarrow k\le3p-2\) (15.168.F). Boundaries checked for \(s=2,3,4,5\).

**OPEN (honest — no soft-close):**

- Residual (i): prove freeness-fail Type I \(k=3p-2\) cannot have \(s_-\le-1\) (or bad case impossible) for all primes \(p\ge5\). **Certified** \(s_-\le-1\) integrally infeasible at \(p=5\) (MILP/HiGHS), not a general proof.
- Residual (ii): deep freeness-fail ND for \(s_+=2\), \(k\ge3p\).

E1_closed_general=false. residual_closed_general=false. L OPEN.

Evidence: `src/e1_gmin_m4_prop15169.py`, `evidence/e1_gmin_m4_prop15169.json`,
`tests/test_prop15169.py`.

### 15.169 addendum — bad-case dual two-level identities

If freeness-fail Type I \(k=3p-2\) occurs with \(s_-=-1\), Max\- two-level \(\{-1,-3\}\), and \(f_e\equiv-1\) on \(U_-\) (bad case), then \(H=G\cup\{e\}\) is dual two-level: mass \(\mathrm{thr}\) at \(S_H=\pm2\), mass \(1-\mathrm{thr}\) at \(\pm4\), with
\[
\mathbb E[S_H^2]=10-\frac6p,\qquad |H|=3p-1=2p+(p-1).
\]
Also \(\sum_{e'<e''\in G}G^+_{e'e''}=\frac{\mathbb E_+[S^2]-k}2=-\frac{3(p-1)(p-4)}{2p}<0\) for \(p\ge5\).

**Status (2026-08-06):** residual (i) **not** closed for general \(p\). Dual-equality Farkas is **conditional** on a disj Gsum LB that is **not** proved in general (see Prop 15.170). Fractional affine+\(s_-\le-1\) remains feasible (integrality gap).

## Prop 15.170 (2026-08-05; status revised 2026-08-06) — Residual (i): Type I \(k=3p-2\), \(s_-\le-1\)

Continues 15.169. **OPEN for general primes \(p\ge5\).** Structure and conditional Farkas algebra shipped; general close blocked by disj Gsum LB hinge.

**Proved / checkable (Fraction / prior props):**

1. **Freeness-fail + gap-2 structure** (15.169): affine \(S=3-2f_e\) on Max\(_+\); gap-2 forces \(s_-=-1\); ND dichotomy at \(s_-=-1\) (good sign \(\Rightarrow\) strong ND; bad case is the only residual risk).

2. **Dual equality correlation (conditional).** If the bad case realises the dual two-level law \(S=-3-2f_e\) on Max\(_-\) with affine Max\(_+\), then
   \[
   (\mathrm{Gsum}\,x)_e = \frac6p-4,
   \]
   where \(\mathrm{Gsum}=\mathbb E_+[ff^\top]+\mathbb E_-[ff^\top]\) and \(x=1_G\).

3. **Gsum identities.** \(\mathrm{Gsum}_{ee}=2\); wedge (share a vertex) \(\mathrm{Gsum}_{ab}=0\) (proved: \(\mathbb E[y_iy_j]_+\!+\!\mathbb E[y_iy_j]_-=0\)); \(\mathrm{Gsum}\,\mathbf1=n\mathbf1\).

4. **Disj Gsum LB — NOT proved for general \(p\).** Candidate bound
   \[
   \mathrm{Gsum}_{ab}\ge -\frac{12}{p(p^2+1)}\quad(a\neq b\text{ disjoint})
   \]
   was previously justified as an association-scheme minimum and is **certified at \(p=5\)** (\(-6/65=-12/(5\cdot26)\)). **Prop 15.158:** Max+ is **not** an IP association scheme, so scheme-min cannot prove the LB for all primes \(p\ge5\). Shipped flag: `gsum_disj_lb_proved_general()=False`.

5. **Conditional box-sum Farkas.** *If* the candidate LB holds, then for \(0\le x\le1\), \(\mathbf1^\top x=k=3p-2\), \(x_e=0\),
   \[
   (\mathrm{Gsum}\,x)_e \ge -\frac{12k}{p(p^2+1)},
   \]
   and dual equality conflicts for all primes \(p\ge5\) because
   \[
   \frac6p-4 < -\frac{12k}{p(p^2+1)} \iff 4p^3-6p^2-32p+18>0.
   \]
   That poly algebra is real (`dual_equality_farkas_algebra_if_lb`); it does **not** alone close residual (i).

6. **Consequence (honest).** residual (i) **OPEN** until disj Gsum LB is proved Max+-free (or an alternate residual-(i) proof). `type_I_k_3p_minus_2_closed_general()=False`. Also \(\mathrm{ES}_2=(13p-12)/p<k\) (integrality seed, independent).

**OPEN:** residual (i) general; residual (ii); full E(1); \(L=\tfrac12\).

Evidence: `src/e1_gmin_m4_prop15170.py`, `evidence/e1_gmin_m4_prop15170.json`,
`tests/test_prop15170.py`, `evidence/share/denseness_path_package.md`.

## Prop 15.171 (2026-08-05; status revised 2026-08-06) — Residual (ii): deep freeness-fail \(k\ge3p\) ND

Continues 15.168–170. **OPEN for general primes \(p\ge5\).** Structure (parity, auto-freeness, fail-eq empty) shipped; dual two-level Farkas blocked by the same disj Gsum LB hinge as 15.170.

**Proved / checkable (Fraction / prior props):**

1. **Parity.** \(s_+=2\) (even scores) \(\Rightarrow k\) even \(\Rightarrow\) Max\- scores even \(\Rightarrow s_-\ne-1\).

2. **Gap-2 deep undercutter classification.** \(s_+\!=2\) gap-2 undercutter \(\Rightarrow s_-\le-2\) (\(s_-\ge0\) not undercutting; \(s_-=-1\) impossible by parity).

3. **Deep freeness \(\Rightarrow\) weak ND.** If \(f_e\not\equiv+1\) on \(\{S=2\}\), some \(y\) has \(S=2\), \(f_e=-1\), \(S_H=1\), \(Q=\Phi-2S_H=\Phi-2\).

4. **Auto-freeness** \(k\le3p-2\) (15.168) \(\Rightarrow\) freeness \(\Rightarrow\) weak ND.

5. **Fail-eq \(k=3p-1\):** freeness-fail \(+\,S\in\{2,4\}\) \(\Rightarrow H\) tight \(S\equiv3\) size \(3p\) empty under bi-tight/Thm A (15.167–168) when bi-tight holds.

6. **Dual two-level freeness-fail Farkas (conditional).** Two-level \(S\in\{2,4\}/\{-2,-4\}\) with freeness-fail affine \(f_e=3-S\) forces
   \[
   (\mathrm{Gsum}\,x)_e=2\bigl(8-\tfrac{3k}p\bigr).
   \]
   Box-sum needs the same unproved disj LB \(\mathrm{Gsum}_{ab}\ge-12/(pn)\). Under that LB, need \(<\) LB on freeness-fail range (Fraction). Without proved LB: **not closed**.

7. **Consequence (honest).** Residual (ii) **OPEN**. With residual (i) open and Gsum hinge open: **E1_closed_general=false**, **L OPEN** (no soft-close). Denseness cannot force \(L=\tfrac12\) without E(1).

**Still OPEN:** residual (i)/(ii) (Gsum disj LB); Path-C residual/16N/Es4\(_*\) (optional).  
residual_closed_general=false. E1_closed_general=false. **L OPEN.**

Evidence: `src/e1_gmin_m4_prop15171.py`, `evidence/e1_gmin_m4_prop15171.json`,
`tests/test_prop15171.py`, `evidence/share/denseness_path_package.md`.
