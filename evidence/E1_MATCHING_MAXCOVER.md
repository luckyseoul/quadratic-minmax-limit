# Matching Max-covers on the ρ=1 Paley family

**Date:** 2026-07-27  
**Status:** Certified numerics; structural. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Fact

Any strict undercutter of a \(\rho=1\) conference must satisfy \(S_F(y)\ge1\) for every \(y\in\mathrm{Max}_{+}\) (integer drop of \(Q\) by at least 2). So every undercutter is a **Max-cover**.

## Perfect matchings as Max-covers

| \(p\) | \(n\) | Random PM best \(\min_y S_M(y)\) | SA max \(\min_y S\) (12 seeds) | Undercutting PMs |
|------:|------:|--------------------------------:|--------------------------------:|:-----------------|
| 3 | 10 | \(+1\) (and 144 exact undercutters) | (exact: 144 of 945) | **144** with \(\Phi=13\) |
| 5 | 26 | \(-3\) (2000 samples) | **all \(\le -1\)** | **none found** |

At \(n=26\), SA (3000 iters, 2-swaps) maximising \(\min_y S_M(y)\) over perfect matchings never exceeds \(-1\). Consistent with: **no perfect matching is a Max-cover for \(p=5\)**, hence the \(n=10\) matching undercut does not lift.

Fractional matching-polytope LP still has value \(p=5\) (Prop 15.27), so the obstruction is integrality of the matching cover.

## Halfspace check

For the halfspace boolean evec \(h\): random PM gives \(\mathbb E[S_M(h)]\approx(p^2+1)/(2p)\); undercutting PMs at \(p=3\) have \(S_M(h)\in\{1,5\}\).

## Relevance to E(1)

- If for all \(p\ge5\) no matching Max-cover exists, matching dichotomy reduces to \(k_\star=0\) (conference exact-optimal) on those orders.
- Combined with \(k_\star=5=O(n^{3/2})\) at \(p=3\), Max-Lipschitz Prop 15.27 would give E(1) on the whole \(\rho=1\) family **if** the \(p\ge5\) claim is proved and non-matching Max-covers are shown not to undercut (they spike, or don't exist at card-min).
- **Not proved.** Do not mark settled.
