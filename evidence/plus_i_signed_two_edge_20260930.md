# Signed two-edge census of the plus-I lifts

**Date:** 2026-09-30
**Status:** finite exact check. Not an upper bound on \(\Phi(C)-m_n\), and not a status change. The limit stays open.

## What was already recorded

`evidence/ns_port_n26_undercut_20260912/` already proves that every one of the 52,650 two-edge flips of the order-26 plus-I lift raises the maximum from 61 to 65, and that every one of the 325 single flips raises it to 63. That unsigned statement is not repeated here as a new result.

Code search on `main` for the signed split (`+65`, `one-sided`, `both signs`) returned nothing, and the same statement is not in the tip commits of the other branches (`research/paired-disagreement-floor-20260924`, `residual/p13-u6-common-moments`, and the older archive, codex, maxplus, navier, and strategy branches).

## What is new

Flipping two edges changes every state by \(-4\), \(0\), or \(+4\). So from a maximizer of value \(\Phi\), the new positive maximum is \(\Phi+4\) exactly when some positive maximizer moves by \(+4\), and at most \(\Phi+2\) otherwise; likewise on the negative side.

Order 10, the plus-I lift of the Paley matrix of \(\mathbf F_5\), \(\Phi=13\), all 990 pairs, exhaustive on the \(2^{9}\) states with first coordinate \(+1\):

| new maximum | new minimum | pairs |
| --- | --- | --- |
| \(+17\) | \(-17\) | 840 |
| \(+17\) | \(-15\) | 75 |
| \(+15\) | \(-17\) | 75 |

None stay at 13, and none miss both sides.

Order 26, the plus-I lift of the Paley matrix of \(\mathbf F_{13}\), \(\Phi=61\). Every one of the 52,650 pairs moves some \(+61\) state to \(+65\) and some \(-61\) state to \(-65\). Since a two-edge flip cannot move any state by more than \(+4\), both sides are exactly \(\pm 65\). There is no one-sided pair.

## What this does not say

Both orders are strict local minima of the Boolean maximum under one- and two-edge flips. The gap below the conference envelope is still only a lower bound on \(\gamma_p\). It is not the \(o(p^3)\) upper bound on \(\Phi(C)-m_n\) that would force the limit to \(1/2\).
