# E(1) star-reduction / submatrix probe

**Date:** 2026-07-27  
**Status:** Structural facts certified; **does not close E(1)**. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Setup (Prop 15.25)

Every Seidel matrix is switching-equivalent to one with first row \((0,+1,\ldots,+1)\). On that slice
\[
\Phi(A)=f(B):=\max_{x\in\{\pm1\}^{n-1}}\bigl(|Q_B(x)|+|\textstyle\sum x_i|\bigr),
\]
and \(m_n=\min_B f(B)\). For Paley \(C\) already in first-row form (as constructed by `paley_conference_prime_power`), the principal submatrix \(C'=C[1\!:\!,1\!:]\) satisfies \(f(C')=\Phi(C)\).

## Certified at \(n=10\) (\(p=3\))

| Object | Value |
|--------|------:|
| \(\Phi(C)\) | 15 |
| \(\Phi(C')\) | 12 \(=m_9\) |
| \(f(C')\) | 15 \(=\Phi(C)\) |
| \(m_{10}\) | 13 |

So \(C'\) is \(\Phi\)-optimal of order 9, while the star objective \(f\) at \(C'\) still equals the conference value 15. The gap \(m_{10}=13<15\) is realised by some other \(B^*\) with \(f(B^*)=13\).

### Matching undercutters \(\to B^*\) distance

All **144** perfect-matching undercutters of \(C\) (N10-S), after first-row renormalisation, give
\[
d_H(B,C')=12
\]
constantly (not \(O(1)\)). Examples:
- match \(\{(0,1),(2,4),(3,5),(6,7),(8,9)\}\): \(d_H(B,C')=12\), \(f(B)=13\).

**Obstruction:** star reduction **scrambles** Hamming distance to the Paley submatrix. Sparse undercutters of \(C\) need not yield sparse \(B\) near \(C'\). Recursive control of \(k_\star(n)\) via \(k_\star(n-1)\) is therefore not straightforward.

### Random flips of \(C'\) (5000 trials per \(k\))

| \(k\) | best \(f\) | hits \(f\le13\) |
|------:|----------:|----------------:|
| 0 | 15 | 0 |
| 1 | 17 | 0 |
| 2–7 | ≥15 | 0 |
| 8 | 13 | 1/5000 |
| 9–20 | ≥15 | 0 |

Sparse random search on \(C'\) almost never finds the \(f=13\) layer; the star-optimal \(B^*\) is rare among near-Paley submatrices.

## Max-only covering is insufficient (reconfirmed)

Greedy edge cover of the 12 boolean \(+p\)-eigenvectors with \(S_F(y)\ge1\) for all \(y\in\mathrm{Max}_+\) terminates at **\(|F|=3\)**, but \(\Phi(C\oplus F)=21>15\): non-maximizers spike. Max-covering alone cannot prove E(1) (already blocked in `E1_RIGIDITY_ATTACK.md`).

On the other hand, the 144 undercutting matchings are **exactly** the perfect matchings with \(S_M(y)\ge1\) for every \(y\in\mathrm{Max}_+\) (count match). Spike control is free for those matchings at \(n=10\), not for arbitrary Max-covers.

## Far minimiser best-\(k\) (SA sample)

Random-start SA at \(n=10\) routinely finds \(\Phi=13\) matrices with best-switch distance to Paley in \(\{9,\ldots,16\}\) (mean \(\sim13\)), never the card-min \(k_\star=5\) from cold start. Confirms:
- far optima exist;
- \(k_\star=\min\) is realised by structured matchings, not by typical SA basins.

For \(W=A\circ C\) in best-switch form, \(\Phi(W)=e-2k_{\mathrm{best}}\) on the samples tested (all-ones achieves the max of \(Q_W\)).

## Prop 15.26: local maximality after matching flips

For any matching \(M\) on a \(\rho=1\) conference with \(p\ge3\), every boolean \(+p\)-evec \(y\) satisfies \(y_i(Ay)_i\ge p-2>0\) (coordinate-local max of \(x\mapsto x^\top Ax\)). Proved in `solution.md`.

### Global Max-determination

| Order | Object | Max-determines \(\Phi\)? |
|------:|--------|:-------------------------|
| 10 | all 945 perfect matchings | **yes** (exhaustive) |
| 26 | 30 random perfect matchings | **19/30 only** (MITM exact) |
| 10 | random non-matching \(k=5\) | \(\sim40\%\) |

At \(n=26\), matching flips that fail Max-determination still have \(\Phi\ge75>65=\Phi(C)\) (no undercut in the sample). Min observed matching \(\Phi=75\).

## What this means for E(1)

1. Prop 15.25 is correct and certified, but **does not by itself** force \(m_n=\Phi(C)-o(n^{3/2})\).
2. Star-reduction distance to \(C'\) is **not** a proxy for \(k_\star\) on the original order.
3. Prop 15.26 (local max) is general; global Max-determination is **not** (fails at \(n=26\)).
4. Remaining load-bearing routes (unchanged):
   - prove \(k_\star=O(n)\) / matching–path–cycle dichotomy on \(n=p^2+1\);
   - or direct gap \(\Phi(C)-m_n=o(n^{3/2})\);
   - or permanent relative gap / non-existence pair.

**Do not mark Main Theorem / HANDOFF settled.**
