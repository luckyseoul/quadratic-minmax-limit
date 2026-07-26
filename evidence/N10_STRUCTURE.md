# Structural gap at \(n=10\): perfect-matching undercutters of Paley

**Status:** certified by exhaustive enumeration (exact \(\Phi\), half-cube).  
**Date:** 2026-07-26  
**Scripts:** `src/n10_structure.py`, `src/n10_matching_optima.py`  
**JSON:** `evidence/n10_structure.json`, `evidence/n10_matching_optima.json`

---

## Setup

Let \(C\) be the Paley conference matrix of order \(n=10=3^2+1\)
(`paley_conference_prime_power(3)`). Then \(C\) is conference, \(\rho(C)=1\)
(halfspace boolean eigenvector), and
\[
\Phi(C)=\tfrac12\,n\sqrt{n-1}=15.
\]
Exact table: \(m_{10}=13\) (recorded multi-worker Gray; \(m_{10}<\Phi(C)\)).

Hamming distance = number of undirected edge-sign disagreements.

---

## Theorem N10-S (certified)

1. **Single-edge local optimality.** Every undirected edge of \(C\) is
   maximizer-balanced in the sense of Prop 15.21 (`solution.md`). Hence every
   single edge flip \(A\) satisfies \(\Phi(A)\ge 17=\Phi(C)+2\).
   (12 maximizers in the half-cube; 45/45 edges balanced.)

2. **Hamming threshold.** Among all Seidel matrices at undirected Hamming
   distance \(k\) from this fixed \(C\),
   \[
   \begin{array}{c|cccccc}
   k & 0 & 1 & 2 & 3 & 4 & 5 \\ \hline
   \min\Phi & 15 & 17 & 15 & 17 & 15 & 13
   \end{array}
   \]
   In particular \(\min_{d_H(A,C)\le 4}\Phi(A)=15\) and the first undercutting of
   Paley to the exact optimum \(m_{10}=13\) occurs at **\(k=5\)**.

3. **Perfect-matching classification.** There are \(\binom{45}{5}=1{,}221{,}759\)
   five-edge sets. Exactly **144** of them yield \(\Phi=13\), and **every one
   of those 144 is a perfect matching** of \(K_{10}\) (degree multiset
   \((1^{10})\)). No non-matching 5-edge set undercuts to 13.

4. **Matching census.** \(K_{10}\) has \(10!/(2^5\,5!)=945\) perfect matchings.
   Flipping each on \(C\) and scoring exact \(\Phi\) yields the histogram
   \[
   \Phi\in\{13,17,21,25\}
   \quad\text{with counts}\quad
   (144,\,405,\,360,\,36).
   \]
   So \(144/945=16/105\) of all perfect-matching flips are exact optima.
   Every optimal matching includes the infinity vertex of the Paley
   construction; partners of \(\infty\) are uniform over the nine field
   elements (16 optima per partner).

5. **Spectral fingerprint of Hamming-5 optima.** For a representative
   matching-optimum \(A\),
   \[
   \Phi(A)=13,\quad
   r(A)=\frac{2\Phi(A)}{n\sqrt{n-1}}=\frac{13}{15},\quad
   \delta(A)=\mathrm{tr}(A^4)-n(n-1)^2=320,
   \]
   and \(\|A\|_{\mathrm{op}}/\sqrt{n-1}\approx 1.311>1\) (not conference).

6. **More distant optima.** SA (86 workers × 12k iters, exact \(\Phi\)) finds
   \(\Phi=13\) matrices at Hamming distance **11–16** from the full Seidel
   switching class of \(C\) (class size \(2^9=512\)), with the same
   \(r=13/15\) and often the same \(\delta=320\). So the matching basin is not
   the only optimum; all observed optima share the product ratio \(r=m_{10}/(\tfrac12 n\sqrt{n-1})\).

---

## Consequence for E(1)

- Exact conference optimality **fails** at \(n=10\), but the absolute gap is
  only \(2\):
  \[
  \frac{\Phi(C)-m_{10}}{n^{3/2}}=\frac{2}{10^{3/2}}\approx 0.063.
  \]
  A uniformly \(O(1)\) (or even \(o(n^{3/2})\)) gap along the \(\rho=1\) family
  would still give E(1): \(\alpha_{n_k}\to\tfrac12\) along \(n_k=p_k^2+1\).

- Local edge-flip control (Props 15.20–15.21) **cannot** prove E(1) by itself:
  the optimum sits at Hamming distance 5 (and other optima much farther).
  Any proof must handle matrices outside the \(o(n)\)-edge ball of conference
  matrices, or exploit a global invariant (e.g. minimisation of \(r\)).

- The perfect-matching undercutters are a concrete algebraic-combinatorial
  object at order 10; they do **not** automatically lift. Random
  perfect-matching flips of Paley \(C_{26}\) (\(\Phi=65\)) produced
  \(\Phi\ge 73\) on 86 independent samples (`evidence/n26_matching_probe.json`)
  — matching flips **raise** \(\Phi\) at \(n=26\), opposite to the \(n=10\) gap.

---

## Verification

```bash
cd /home/nick/quadratic-minmax-limit
# full campaign (~4 min, 86 workers)
WORKERS=86 python3 src/n10_structure.py
WORKERS=86 python3 src/n10_matching_optima.py
# load-bearing unit tests (no full binom(45,5) re-sweep)
python3 -m pytest tests/test_minmax.py -k n10_structure -v
```

Expected: maximizer balance true; k=5 first reach; 144 matching optima; tests pass.
