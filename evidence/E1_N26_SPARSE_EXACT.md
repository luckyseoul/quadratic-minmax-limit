# E(1) at \(n=26\): exact MITM sparse-flip survey

**Date:** 2026-07-27  
**Status:** Certified exact-\(\Phi\) samples; **no undercut of \(\Phi=65\)**. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Method

Shipped `phi_mitm` (meet-in-the-middle exact \(\Phi\) for even \(n\le28\)):
fix \(x_0=+1\), enumerate each half of the remaining coordinates, combine.
At \(n=26\): \(2^{12}\times2^{13}\) states, \(\sim0.3\)s per matrix.

## Results (Paley \(C_{26}\), \(\Phi=65\))

| Family | samples | min \(\Phi\) | undercuts |
|:-------|--------:|------------:|----------:|
| Paley itself | 1 | 65 | — |
| Random perfect matchings | 30 | **75** | 0 |
| Random cycles \(C_4,\ldots,C_{26}\) | 40 each | ≥69 | 0 |
| Random \(k\)-edge sets, \(k\in\{1,2,3,5,8,10,13,15,20\}\) | 30 each | ≥67 | 0 |
| Stars degree \(1..12\) | 8 each | ≥67 | 0 |

Single-edge flips: \(\Phi=67=\Phi(C)+2\) on all samples (consistent with maximizer balance / Prop 15.21).

## Interpretation

- No sparse (path / cycle / star / random \(k\le20\)) undercutter of Paley at \(n=26\) in this exact census.
- Matching undercutters of \(n=10\) **do not lift**: random PM flips raise \(\Phi\) by \(\ge10\).
- Consistent with \(k_\star=0\) (Paley optimal) at \(n=26\), but **not a proof** of E(1).
- Combined with N10-C6: at the only known undercutting \(\rho=1\) order (\(n=10\)), closest undercutters have \(k_\star=5=O(n)\).

## Code

- `phi_mitm` in `src/minmax_quadratic.py`
- Prior intensive SA: `evidence/e1_n26_intensive.json` (no undercut; used exact rescore)

**Do not mark Main Theorem settled.**
