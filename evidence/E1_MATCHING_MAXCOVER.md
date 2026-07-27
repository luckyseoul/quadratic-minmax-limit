# Matching Max-covers on the ρ=1 Paley family

**Date:** 2026-07-27  
**Status:** Certified numerics; structural. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Fact

Any strict undercutter of a \(\rho=1\) conference must satisfy \(S_F(y)\ge1\) for every \(y\in\mathrm{Max}_{+}\) (integer drop of \(Q\) by at least 2). So every undercutter is a **Max-cover**.

## Perfect matchings as Max-covers

| \(p\) | \(n\) | Random PM best \(\min_y S_M(y)\) | SA max \(\min_y S\) (12 seeds) | Undercutting PMs |
|------:|------:|--------------------------------:|--------------------------------:|:-----------------|
| 3 | 10 | \(+1\) (and 144 exact undercutters) | (exact: 144 of 945) | **144** with \(\Phi=13\) |
| 5 | 26 | \(-3\) (5000 samples) | **all \(\le -1\)** (20×15k-iter SA) | **none found** |

At \(n=26\), an earlier weaker SA never exceeded \(\min S=-1\). A stronger 2-swap SA **does** find Max-covers (\(\min S=1\)); those covers are two-sided but exact MITM \(\Phi=65=\Phi(C)\) (non-max spike). See **`E1_MATCHING_COVER_SPIKE.md`** (supersedes the “no matching Max-cover” claim). Note \(n/2=13\) is odd, so \(S_M(y)\) is always odd.

Fractional matching-polytope LP still has value \(p=5\) (Prop 15.27), so the obstruction is integrality of the matching cover.

## Halfspace check

For the halfspace boolean evec \(h\): random PM gives \(\mathbb E[S_M(h)]\approx(p^2+1)/(2p)\); undercutting PMs at \(p=3\) have \(S_M(h)\in\{1,5\}\).

## Relevance to E(1)

- If for all \(p\ge5\) no matching Max-cover exists, matching dichotomy reduces to \(k_\star=0\) (conference exact-optimal) on those orders.
- Combined with \(k_\star=5=O(n^{3/2})\) at \(p=3\), Max-Lipschitz Prop 15.27 would give E(1) on the whole \(\rho=1\) family **if** the \(p\ge5\) claim is proved and non-matching Max-covers are shown not to undercut (they spike, or don't exist at card-min).
- **Not proved.** Do not mark settled.
