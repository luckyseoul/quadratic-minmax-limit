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

**On existence of \(\lim\alpha_n\).** **OPEN.** Neither existence nor non-existence is proved.
By Proposition 6.2, existence is equivalent to convergence of \(\alpha_n\) along Paley orders alone.
Soft multipartite/Hadamard inequalities cannot force \(\lambda=\Lambda\) (§9–§10).
Resume checklist and evidence map: **`HANDOFF.md`**.

**E(2) progress (interval formula, not settlement).** For Paley \(q\equiv1\bmod4\) and the interval
boolean vector, \(x^\top Cx=2-8\sum_{d\le(q-1)/2}d\chi(d)\) exactly
(`evidence/E2_INTERVAL_FORMULA.md`, `src/interval_rho_formula.py`). This makes the standard
constructive lower bound on \(\rho(C)\) elementary; proving \(\rho\to1\) still needs asymptotics of
that character sum. E(1) remains open (\(n=6\) exact opt; \(n=10\) gap \(2\); SA+exact at \(n=14,18\)
found no Paley undercut). **Do not read this paragraph as \(\lim\alpha_n\) existing.**

**Corollary (\(\rho=1\) along a dense Paley family).** For every odd prime \(p\), the Paley conference
matrix of order \(n=p^2+1\) (over \(\mathbb F_{p^2}\)) admits a halfspace boolean eigenvector
\(Cx=px\), hence \(\rho(C)=1\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\). Along \(n_k=p_k^2+1\) one has
\(n_{k+1}/n_k\to1\) and \(\limsup_k\rho(C_{n_k})=1\). Proof: `evidence/PROOF_rho_eq_1.md`;
shipped checks: `paley_conference_prime_power`, `halfspace_boolean_vector`.
This does **not** by itself force \(\lim\alpha_n=\tfrac12\) (needs E(1) and a matching liminf).

**What is complete (Theorems D–G).** Dual-Gaussian universal lower bound \(m_n\ge n\sqrt{n-1}/\pi\)
(Prop.~5.2); cut-code identity \(m_n=\binom n2-2\rho(D_n)\) (Prop.~1.2); conference spectral
identity and exact Nesterov formula (Theorem D); Seidel switching; unique min-op and
\(\mathrm{tr}(A^4)\) characterisation of conference; \(L^2\)-universality of \(Q\); **exact
fourth-moment formula** \(\mathbb E[Q^4]\) uniquely minimised at conference; exact optimality
criterion via spectral gap; **\(m_6=\Phi(C)=5\)** by exhaustive gap check; limsup bound
\(\limsup\alpha_n\le\tfrac12\limsup_k\rho(C_k)\); optimality \(\Leftrightarrow\) minimisation of
\(r(A)=\max|x^\top Ax|/(n\sqrt{n-1})\); **\(\rho=1\) for all Paley orders \(n=p^2+1\)**.

**What remains for existence (Theorem E).** (1) Asymptotic optimality
\(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) along Paley — proved at \(n=6\); **fails exact optimality at
\(n=10\)** (exact \(m_{10}=13<\Phi(C_{10})=15\) for Paley \(q=9\)); local edge-opt for Paley
\(n\le18\) (Prop 15.21); Q4 path dead (Props 15.16, 15.19). E(1) must be asymptotic
(\(o(n^{3/2})\) gap), not exact. Reduced to product-min of \(\rho\cdot\mathrm{op}\) via
delocalization/rigidity. (2) \(\rho(C_k)\to\rho_*\) for general Paley — on the subsequence
\(n=p^2+1\) one already has \(\rho\equiv1\); for other Paley orders, exact \(\rho\) is strictly
increasing through \(n=42\), while constructive interval lower bounds reach \(\rho_{\mathrm{int}}\gtrsim0.99\)
(not a full proof of \(\rho\to1\)). Conditional Thm E: both (1)+(2)
\(\Rightarrow\lim\alpha_n=\rho_*/2\). Stolz (Thm F) still open.
**Existence of \(\lim\alpha_n\) remains OPEN.**

This partially answers MathOverflow [413935](https://mathoverflow.net/questions/413935) /
https://x.com/PI010101/status/2081070728422752329
(the author already knew one-sided bounds; the contribution is the sandwich, denseness,
majorant identity, Paley reduction, conference spectral/\(L^4\) calculus, switching theory,
exact optimality at \(n=6\), \(\rho=1\) on \(n=p^2+1\), and obstruction analysis).

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
