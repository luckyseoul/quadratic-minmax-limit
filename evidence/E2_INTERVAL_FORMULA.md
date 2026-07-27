# Exact interval formula for Paley \(\rho\) (E(2) constructive lower bound)

**Status:** proved identity; does **not** settle \(\lim\alpha_n\).  
**Code:** `src/interval_rho_formula.py`  
**Census:** `evidence/e2_interval_rho.json` (shift-optimized numerical \(\rho_{\mathrm{int}}\))

Existence of \(L=\lim\alpha_n\) remains **OPEN**.

---

## Setup

Let \(q\) be prime, \(q\equiv1\pmod4\), \(\chi\) the Legendre symbol mod \(q\), and \(C\) the
Paley conference matrix of order \(n=q+1\) (vertices \(\{\infty\}\cup\mathbb F_q\)).
The **interval boolean vector** is
\[
x_\infty=+1,\qquad
x_t=+1\iff t\in\{0,1,\ldots,\lfloor q/2\rfloor\}\subset\mathbb F_q.
\]
Write \(\rho_{\mathrm{int}}(C,x)=|x^\top Cx|/(n\sqrt{n-1})\le\rho(C)\).

---

## Theorem (interval character-sum formula)

\[
x^\top C x
=
2-8\,\Sigma_q,
\qquad
\Sigma_q
:=
\sum_{d=1}^{(q-1)/2} d\,\chi(d).
\]
Hence
\[
\rho_{\mathrm{int}}(C,x)
=
\frac{|2-8\Sigma_q|}{(q+1)\sqrt q}.
\]

### Proof

**Step 1 (infinity contribution).** The first row/column of \(C\) is \(+1\) on
\(\mathbb F_q\). With \(x_\infty=+1\),
\[
\sum_{t\in\mathbb F_q}2\,x_t
=
2\sum_{t}x_t.
\]
The interval has \(m=\lfloor q/2\rfloor+1=(q+1)/2\) pluses and \((q-1)/2\) minuses, so
\(\sum_t x_t=1\). Contribution \(=2\).

**Step 2 (field–field contribution).** Off-diagonal field entries are
\(C_{a,b}=\chi(b-a)\). Hence
\[
\sum_{a\neq b}\chi(b-a)\,x_ax_b
=
\sum_{d\neq0}\chi(d)\,c(d),
\qquad
c(d)=\sum_{a\in\mathbb F_q}x_ax_{a+d}.
\]

**Step 3 (interval autocorrelation).** For the cyclic interval signing,
\[
c(d)=q-4\min(d,q-d)\qquad(d=1,\ldots,q-1).
\]
(Verified by direct count of interval overlaps; shipped tests recheck for many \(q\).)

**Step 4 (character sum collapse).** \(\sum_{d=1}^{q-1}\chi(d)=0\), so
\[
\sum_{d\neq0}\chi(d)\,c(d)
=
-4\sum_{d=1}^{q-1}\chi(d)\,\min(d,q-d).
\]
Pair \(d\) with \(q-d\). Since \(q\equiv1\pmod4\), \(\chi(-1)=1\), and
\[
\sum_{d=1}^{q-1}\chi(d)\,\min(d,q-d)
=
2\sum_{d=1}^{(q-1)/2}d\,\chi(d)
=
2\Sigma_q.
\]
Thus the field contribution is \(-8\Sigma_q\).

**Step 5.** Total \(x^\top Cx=2-8\Sigma_q\). Divide by \(n\sqrt{n-1}=(q+1)\sqrt q\). \(\square\)

---

## What this does *not* prove

- **Not** \(\rho(C)\to1\): one still needs
  \(\Sigma_q/q^{3/2}\to-1/8\) (or \(\limsup|\Sigma_q|/q^{3/2}=1/8\)) from analytic number
  theory / \(L\)-functions. Numerics show \(-\Sigma_q/q^{3/2}\in[0.08,0.122]\) on primes to
  \(10^3\), with record \(\rho_{\mathrm{int}}\approx0.978\) at \(n=1622\) (shift-optimized census).
- **Not** E(1) (asymptotic optimality of \(m_n\) vs \(\Phi(C)\)).
- **Not** existence of \(\lim\alpha_n\).

Conditional Theorem E still requires E(1)+E(2). This note only makes the standard
interval lower bound on \(\rho\) **exact and elementary**.

---

## Verification

```bash
python3 -m pytest tests/test_minmax.py -k interval_rho -v
python3 src/e2_interval_rho.py   # multi-q census with shifts
```
