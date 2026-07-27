# Theorem N10-C: algebraic classification of the 144 optimal matchings

**Status:** certified by exhaustive enumeration + group orbit  
**Date:** 2026-07-27  
**Script:** `src/n10_matching_classify.py`  
**JSON:** `evidence/n10_matching_classify.json`  
**Depends on:** Theorem N10-S (`evidence/N10_STRUCTURE.md`)

Existence of \(L=\lim\alpha_n\) remains **OPEN**. This note classifies the
Hamming-5 undercutters of Paley \(C_{10}\); it does not settle the limit.

---

## Setup

Paley conference matrix \(C\) of order \(n=10=3^2+1\) over \(\mathbb F_9\),
vertices \(\{\infty\}\cup\mathbb F_9\) (library encoding: index \(0=\infty\),
\(1..9\) = field elements \(0..8\) as \(a+b\cdot3\)). \(\Phi(C)=15\), \(\rho=1\),
and \(m_{10}=13\). Theorem N10-S: the only 5-edge undercutters of \(C\) are
**144 perfect matchings** of \(K_{10}\).

---

## Theorem N10-C (certified)

### (i) Maximizer-drop criterion

Let \(M\) be a perfect matching of \(K_{10}\) and
\[
S_M(x)=\sum_{\{i,j\}\in M}C_{ij}\,x_i x_j.
\]
Write \(Q_C(x)=\tfrac12 x^\top C x\). After flipping the edges of \(M\),
\[
Q_A(x)=Q_C(x)-2\,S_M(x).
\]
There are exactly **12** maximizers of \(\Phi(C)\) in the half-cube \(x_0=+1\)
(6 with \(Q_C=+15\), 6 with \(Q_C=-15\)).

**Criterion.** \(M\) is an optimal undercutter (\(\Phi(A)=13\)) if and only if
for every maximizer \(x\),
\[
\mathrm{sign}\bigl(Q_C(x)\bigr)\,S_M(x)\;\ge\;1.
\]
Equivalently: every old maximizer drops to \(|Q_A|\le13\).

**Census (exact \(\Phi\) on all 945 perfect matchings):**

| filter | count |
|--------|------:|
| \(\Phi\) after flip \(=13\) | 144 |
| maximizer-drop criterion (all 12) | 144 |
| criterion using only the 6 positive maximizers | 144 |
| set equality of the three collections | yes |

So the six maximizers with \(Q_C=+15\) already form a complete certificate.

### (ii) Sign product (necessary, not sufficient)

For every optimal matching,
\[
\prod_{\{i,j\}\in M}C_{ij}=-1.
\]
There are 504 perfect matchings with this sign product; only 144 are optimal.
Sign product is therefore necessary but not sufficient.

### (iii) Single \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)-orbit

Identify the vertex set with \(\mathrm{PG}(1,9)=\{\infty\}\cup\mathbb F_9\).
The group \(\mathrm{P}\Gamma\mathrm{L}(2,9)=\mathrm{PGL}(2,9)\rtimes\mathrm{Gal}(\mathbb F_9/\mathbb F_3)\)
acts by Möbius transformations \(z\mapsto(az+b)/(cz+d)\) and the Frobenius
\(u\mapsto u^3\) on field points.

- \(|\mathrm{PGL}(2,9)|=720\) (verified by enumeration of matrix classes).
- The orbit of any one optimal matching under \(\mathrm{P}\Gamma\mathrm{L}(2,9)\)
  has size **144** and equals the full set of optimal matchings.

Thus the 144 optimal matchings form a **single orbit** under the natural
Paley automorphism group of the point set.

---

## Consequence for E(1)

- The Hamming-5 optima are not sporadic: they are the \(\mathrm{P}\Gamma\mathrm{L}\)-orbit
  of any one maximizer-dropping perfect matching.
- The absolute gap remains \(2\); relative gap \(2/n^{3/2}\to0\) if \(O(1)\) gaps
  persist along \(\rho=1\) orders — still consistent with E(1).
- Random matchings at \(n=26\) raise \(\Phi\) (N10-S); the orbit description does
  not by itself produce a lifting construction. E(1) remains open.

---

## Verification

```bash
cd /home/nick/quadratic-minmax-limit
python3 src/n10_matching_classify.py
python3 -m pytest tests/test_minmax.py -k 'n10_matching_classify or n10_matching' -v
```

Expected: 144/144 agreement; single orbit of size 144; tests pass.
