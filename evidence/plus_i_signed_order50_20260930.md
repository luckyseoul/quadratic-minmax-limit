# Signed one- and two-edge census of the order-50 plus-I lift

**Date:** 2026-09-30
**Status:** finite exact check, conditional on the already recorded value \(\Phi(K)=169\). Not an upper bound on \(\Phi(C)-m_n\). The limit stays open.

## What was already recorded

`evidence/NOTE_2026-09-19_PLUS_I_COHERENT_LIFT.md` records that the plus-I lift of the Paley Seidel matrix of \(\mathbf F_{25}\) has exact \(\Phi=169\) (nuka RX 9070 XT hipBLAS, \(2^{49}\) projective states) against the Paley conference value 175. `evidence/plus_i_signed_two_edge_20260930.md` records the signed two-edge census at orders 10 and 26 only. The unsigned order-26 scan in `evidence/ns_port_n26_undercut_20260912` is not repeated here.

The order-26 block in that census is the exact \(m_{13}=20\) block of the plus-I note, not a matrix the note identifies as the Paley matrix of \(\mathbf F_{13}\).

## What is new

One edge changes every state by \(\pm 2\), and two edges change every state by \(-4\), \(0\), or \(+4\). From the recorded maximizers below, every edge and every pair is outward on both signs.

Witnesses used: 117 states at \(+169\) and 117 states at \(-169\), recomputed as \(x^T K x/2\) from `evidence/coherent_BI_lift_20260919/B25.txt`. These are witnesses, not a claim that the maximizer list is complete. The value \(\Phi(K)=169\) is the recorded sweep, not a new one.

| flip | pairs | positive side | negative side |
| --- | --- | --- | --- |
| one edge | 1225 | all reach \(+171\) | all reach \(-171\) |
| two edges | 749700 | all reach \(+173\) | all reach \(-173\) |

Outward degree on each witness is \(528=(1225-169)/2\). Pair multiplicities are positive on both signs (no holes). Lipschitz then forces the neighbor maxima to be exactly 171 and 173. Both sit strictly below the conference envelope 175.

## What this does not say

The order-50 lift is a strict local minimum under one- and two-edge flips, with both signs attained. The gap \(6=p-1\) remains a lower bound on \(\gamma_p\). It is not the \(o(p^3)\) upper bound on \(\Phi(C)-m_n\) that would force the limit to \(1/2\).
