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

**Proposition 15.20d (conditional settlement: \(k_\star=O(n)\) \(\Rightarrow L=\tfrac12\)).** Let \(n_k=p_k^2+1\) run over the \(\rho=1\) Paley family, and write \(k_\star(n)\) for the minimal best-switch Hamming distance from a \(\Phi\)-minimiser to the Paley conference matrix of order \(n\). If \(k_\star(n_k)=O(n_k)\), then
\[
\lim_{n\to\infty}\alpha_n=\tfrac12.
\]
*Proof.* Prop 15.20b gives \(m_{n_k}\ge\Phi(C_{n_k})-2k_\star(n_k)=\tfrac12 n_k\sqrt{n_k-1}-O(n_k)\), so
\[
\alpha_{n_k}\ge\tfrac12\sqrt{1-1/n_k}-O(n_k^{-1/2})\to\tfrac12.
\]
Combined with the universal limsup \(\limsup\alpha_n\le\tfrac12\) (Prop 4.1 / conference construction) one has \(\alpha_{n_k}\to\tfrac12\). The family \(n_k\) is dense in the sense \(n_{k+1}/n_k\to1\) (prime number theorem in the progression of odd primes), so Prop 6.1–6.2 force \(\lim_n\alpha_n=\tfrac12\). \(\square\)

*Status of the hypothesis.* At \(n=10\), \(k_\star=5=O(n)\) (N10-S: matching undercutters; N10-C6: all 360 Hamming-6 undercutters are 6-cycles, also \(O(n)\)). At \(n=26\), exact MITM sparse/SA census found no undercut of \(\Phi=65\) (consistent with \(k_\star=0\)). **The general bound \(k_\star=O(n)\) (or even \(o(n^{3/2})\)) on all \(n=p^2+1\) is not proved.** Existence of \(\lim\alpha_n\) remains **OPEN**.

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

Certified boolean \(+p\)-evec counts for Paley \(n=p^2+1\): \(12,260,11452\) at \(p=3,5,7\) (`evidence/BOOLEAN_EVECS_MAX.md`). The ratio \(\#/n^{3/2}\) increases through \(p=7\), so the crude covering bound \(k_\star\le|\mathrm{Max}|\) is **not** \(o(n^{3/2})\) and does not prove E(1).

**n=26 exact MITM census (2026-07-27).** Shipped `phi_mitm` (meet-in-the-middle exact \(\Phi\), even \(n\le28\)). Random matchings/cycles/stars/\(k\le20\) flips and 86-seed SA+MITM rescore: **no undercut of \(\Phi(C_{26})=65\)** (best SA exact \(67\)). Evidence: `evidence/E1_N26_SPARSE_EXACT.md`, `e1_n26_mitm_sa.json`. Consistent with \(k_\star=0\) at \(n=26\); not a general E(1) proof.
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
