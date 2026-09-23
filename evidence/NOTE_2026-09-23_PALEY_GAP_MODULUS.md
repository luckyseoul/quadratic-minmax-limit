# One-step modulus for the Paley gaps

2026-09-23. Original limit OPEN. This does not prove that \(\gamma_p\)
converges.

## 1. The comparison

Let \(p<q\) be odd primes, \(n=p^2+1\), \(N=q^2+1\), and let \(C_n\),
\(C_N\) be the Paley conference matrices of those orders. Write
\[
\gamma_p=\frac{\Phi(C_n)-m_n}{n^{3/2}}
=\tfrac12\sqrt{1-\frac1n}-\alpha_n,
\]
and likewise \(\gamma_q\). CORE §2 gives \(m_n\le m_N\), and
\(m_n\le\Phi(C_n)=\tfrac12 n\sqrt{n-1}\), so
\(\alpha_n\le\tfrac12\sqrt{1-1/n}\).

**Proposition.** 
\[
\gamma_q\le\gamma_p+\tfrac12\left(
\sqrt{1-\frac1N}
-\sqrt{1-\frac1n}\left(\frac nN\right)^{3/2}
\right).
\]
The second summand is positive and depends only on the two orders.

*Proof.* 
\[
\gamma_q-\gamma_p
=\tfrac12\Bigl(\sqrt{1-\tfrac1N}-\sqrt{1-\tfrac1n}\Bigr)
+\alpha_n-\alpha_N.
\]
Restriction gives \(m_N\ge m_n\), hence
\(\alpha_N\ge\alpha_n(n/N)^{3/2}\) and
\[
\alpha_n-\alpha_N
\le\alpha_n\Bigl(1-\bigl(\tfrac nN\bigr)^{3/2}\Bigr)
\le\tfrac12\sqrt{1-\tfrac1n}
\Bigl(1-\bigl(\tfrac nN\bigr)^{3/2}\Bigr).
\]
Adding the square-root term produces the displayed majorant. \(\square\)

## 2. What a definite rise costs

Let \(M(n,N)\) denote the majorant in the proposition. It is increasing
in \(N\) for fixed \(n\), and
\[
M(n,N)=\tfrac34\cdot\frac{N-n}n+O\Bigl(\frac1n+\Bigl(\frac{N-n}n\Bigr)^2\Bigr)
\]
as \(n\to\infty\) with \((N-n)/n\to0\). Therefore a rise
\(\gamma_q\ge\gamma_p+\delta\) with \(\delta>0\) fixed forces
\[
\frac{N}n\ge 1+\tfrac43\delta+o(1).
\]
Along consecutive odd primes one has \(p_{k+1}/p_k\to1\) (the same fact
CORE §8 uses), so
\[
\limsup_k\bigl(\gamma_{p_{k+1}}-\gamma_{p_k}\bigr)\le0.
\]
An upward jump of fixed height cannot occur between one Paley order and
the next, once \(p\) is large. The allowed increments are not summable:
\(\sum(p_{k+1}-p_k)/p_k\) diverges like \(\log p\). Slow oscillation of
\(\gamma_p\) is still possible.

## 3. What restriction does not give

For fixed \(p\), as \(q\to\infty\),
\[
M(p^2+1,\,q^2+1)\to\tfrac12.
\]
The resulting bound \(\gamma_q\le\gamma_p+\tfrac12\) is the trivial one,
since \(0\le\gamma_q<\tfrac12\). Monotonicity of \(m_n\) does not prove
\(\gamma_q\le\gamma_p+\varepsilon(p)\) for every prime \(q\ge Q(p)\) with
an \(\varepsilon(p)\to0\) independent of \(q\). A tail comparison of that
strength needs an input other than restriction.

The crude padding bound is weak. CORE §2 yields
\(m_N\le m_n+dn+d(d-1)/2\) with \(d=N-n\). For a consecutive prime gap
this \(d\) is about \(2p(q-p)\), and that particular estimate has
normalized error of size about the gap itself.

An increment of size \(\Theta(\sqrt n)\) is not, by itself, a
non-summable error. With \(\delta_n=m_{n+1}-m_n\),
\[
\alpha_{n+1}-\alpha_n
=\frac{\delta_n-m_n\bigl((1+1/n)^{3/2}-1\bigr)}{(n+1)^{3/2}}.
\]
The main term \(m_n\bigl((1+1/n)^{3/2}-1\bigr)\) is
\(\tfrac32\alpha_n\sqrt n+O(n^{-1/2})\). It is exactly the increment that
holds \(\alpha\) fixed. If a minimizer admits a one-vertex extension with
\[
\delta_n\le m_n\bigl((1+1/n)^{3/2}-1\bigr)+r_n
\]
and \(\sum r_n n^{-3/2}<\infty\), then \(\alpha_n\) converges, and so
does \(\gamma_p\). A bound \(\delta_n\le\tfrac32 m_n/n+C\) is the case
\(r_n=O(1)\). The absolute CORE padding \(\delta_n\le n\) does not give
this. The extension formulation
\(\min_r\max_x\bigl(|Q(x)|+|x\cdot r|\bigr)\) is the object that would.

## 4. What this does not do

Existence of \(\lim\gamma_p\), the value of that limit, the factor-2
rate \(\Psi\le 2\sqrt2\, m_n+o(n^{3/2})\), and the factor-3 ray are
untouched. Gap-2 covers are not used.
