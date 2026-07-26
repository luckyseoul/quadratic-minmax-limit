# Theorem: ρ = 1 for Paley conference of order \(p^2+1\)

Let \(p\) be an odd prime, \(q=p^2\), \(n=q+1\), and \(C\) the Paley conference matrix over \(\mathbb F_q\).

## Construction

Write \(\mathbb F_q=\mathbb F_p(\omega)\). Let \(L\colon\mathbb F_q\to\mathbb F_p\) be the \(\mathbb F_p\)-linear form \(L(a+b\omega)=b\).
Let \(S=\{0,1,\ldots,\lfloor p/2\rfloor\}\subset\mathbb F_p\), so \(|S|=(p+1)/2\).
Define \(x\in\{\pm1\}^n\) by \(x_\infty=1\) and \(x_u=\sigma(L(u))\) where \(\sigma(c)=+1\) iff \(c\in S\), else \(-1\).

## Proof that \(Cx=px\)

### \((Cx)_\infty=p\)

\((Cx)_\infty=\sum_{u\in\mathbb F}x_u=p|S|-p|S^c|=p\bigl((p+1)/2-(p-1)/2\bigr)=p\).

### \((Cx)_v\) for \(v\in\mathbb F\)

\((Cx)_v=1+\sum_{d\neq0}\chi(d)\,x_{v-d}=1+\sum_{d\neq0}\chi(d)\,\sigma(L(v)-L(d))\).

Let \(c=L(v)\). Group by fibers \(\ell=L(d)\):

\[
\sum_{d\neq0}\chi(d)\,\sigma(c-L(d))
=\sum_\ell\sigma(c-\ell)\sum_{\substack{d:L(d)=\ell\\ d\neq0}}\chi(d).
\]

**Fiber character sums:**

- \(\ker L=\mathbb F_p\subset\mathbb F_q\). For \(t\in\mathbb F_p^\times\), \(\chi(t)=t^{(q-1)/2}=\eta(t)^{p+1}=1\).
  So \(\sum_{d\in\ker\setminus\{0\}}\chi(d)=p-1\).
- For \(\ell\neq0\), fiber \(=d_0+\mathbb F_p\), \(d_0\notin\mathbb F_p\). The sum \(\sum_{t\in\mathbb F_p}\chi(d_0+t)\) is a
  complete character sum of a non-square-discriminant quadratic, equal to \(-1\).

**Assembly:**

\[
\text{sum}=(p-1)\sigma(c)+\sum_{\ell\neq0}(-1)\sigma(c-\ell)
=p\sigma(c)-\sum_\ell\sigma(\ell).
\]

\(\sum_\ell\sigma(\ell)=|S|-|S^c|=1\). Thus sum \(=p\sigma(c)-1\), and
\((Cx)_v=1+p\sigma(c)-1=p\,x_v\).

## Conclusion

\(Cx=px\) with \(x\in\{\pm1\}^n\), so \(|x^\top Cx|=pn=n\sqrt{n-1}\),
hence \(\rho(C)=1\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\).

## Corollary

Along \(n_k=p_k^2+1\) (odd primes \(p_k\)), \(\rho(C_{n_k})=1\) for all \(k\).
Thus \(\limsup\rho(C_n)=1\) along this Paley subsequence, and \(n_{k+1}/n_k\to1\).

## Verification

Shipped: `paley_conference_prime_power`, `halfspace_boolean_vector` in `src/minmax_quadratic.py`.
Checks for \(p=3,5,7\) (\(n=10,26,50\)): `evidence/rho1_verify.json` and pytest `test_rho_eq_1_paley_prime_power`.
