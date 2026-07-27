# Asymptotics of interval \(\rho_{\mathrm{int}}\) and \(\limsup\rho=1\)

**Status:** asymptotic identity checked to high precision; limsup claim along intervals.  
**Does not settle \(\lim\alpha_n\)** (still needs E(1)).  
**Code:** `src/interval_rho_formula.py`, `src/interval_rho_asymptotics.py`

---

## Setup

Paley conference \(C\) of order \(n=q+1\), \(q\) prime \(\equiv1\pmod4\), interval boolean \(x\),
\[
\Sigma_q=\sum_{d=1}^{(q-1)/2}d\,\chi(d),\qquad
x^\top Cx=2-8\Sigma_q,\qquad
\rho_{\mathrm{int}}=\frac{|2-8\Sigma_q|}{(q+1)\sqrt q}.
\]

## Fourier representation

Let \(f(k)=\min(k,q-k)\) on \(\mathbb Z/q\mathbb Z\), \(M=(q-1)/2\), and
\[
\hat f(m)=2\sum_{k=1}^{M}k\cos\Bigl(\frac{2\pi mk}{q}\Bigr).
\]
Then (Gauss sum \(\tau=\sqrt q\) for \(q\equiv1\pmod4\))
\[
\Sigma_q=\frac1{2\sqrt q}\sum_{m=1}^{q-1}\chi(m)\,\hat f(m).
\]
(Verified exactly against \(\Sigma_q\) for many \(q\).)

## Leading asymptotic

For odd \(m=o(q)\),
\[
\hat f(m)=-\frac{q^2}{\pi^2 m^2}+O\Bigl(\frac{q^0\cdot\mathrm{poly}(\log q)}{1}\Bigr)
\]
(continuum Fourier of the tent map; relative error \(O((m/q)^2)\)).  
Pairing \(m\leftrightarrow q-m\) (even \(\leftrightarrow\) odd) doubles the odd-\(m\) contribution. Hence
\begin{equation}
\label{eq:rho-asymp}
\rho_{\mathrm{int}}
=
\frac8{\pi^2}\sum_{\substack{1\le m<q\\ m\text{ odd}}}\frac{\chi(m)}{m^2}+o(1)
=
\frac8{\pi^2}\,L(2,\chi)\Bigl(1-\frac{\chi(2)}{4}\Bigr)+o(1).
\end{equation}
**Numeric check:** for primes \(q\le8009\), \(|\rho_{\mathrm{int}}-(8/\pi^2)L_{\mathrm{odd}}|<0.006\), and \(<5\cdot10^{-4}\) for \(q\ge1009\).

## \(\limsup\rho_{\mathrm{int}}=1\)

Along primes \(q\equiv1\pmod8\) (so \(\chi(2)=1\)) with \(\chi(p)=1\) for all primes \(p\le w\),
\[
L(2,\chi)\ge\prod_{p\le w}(1-p^{-2})^{-1}\cdot\bigl(1-O(w^{-1})\bigr)
\to\zeta(2)\quad(w\to\infty)
\]
by Dirichlet (infinitely many such \(q\) for each fixed \(w\)). Thus
\[
\limsup_q\frac8{\pi^2}L(2,\chi)\Bigl(1-\frac{\chi(2)}{4}\Bigr)
\ge\frac8{\pi^2}\zeta(2)\cdot\frac34
=
1.
\]
Combined with \(\rho_{\mathrm{int}}\le1\) and \eqref{eq:rho-asymp},
\[
\limsup_{\substack{q\equiv1(4)\\q\text{ prime}}}\rho_{\mathrm{int}}(q)=1.
\]
Hence \(\limsup\rho(C_q)=1\) along prime-field Paley orders (also follows from the \(\rho=1\) family \(n=p^2+1\)).

## What remains for \(\lim\alpha_n\)

- **E(2) full:** \(\rho(C_q)\to1\) for *every* Paley sequence (not only limsup). Interval \(\rho_{\mathrm{int}}\) itself oscillates (e.g. \(\sim0.77\) at some large \(q\)).
- **E(1):** \(m_n=\Phi(C_n)+o(n^{3/2})\) along \(\rho=1\) orders \(n=p^2+1\). **Still open.** Along that dense family, E(1) alone \(\Rightarrow\lim\alpha_n=\tfrac12\) (Prop 6.2).
- Maximizers of \(\rho=1\) Paley satisfy \(\mathbb E[yy^\top]=I\) (exact 2-design; checked \(p=3,5\)). Insufficient alone for E(1) (4th moments are not design-rigid; at \(n=10\), \(m_{10}=13<\Phi=15\)).

**Existence of \(\lim\alpha_n\) remains OPEN.**
