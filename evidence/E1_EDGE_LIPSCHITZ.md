# E(1) via edge-counting Lipschitz (Prop 15.20b)

**Date:** 2026-07-27  
**Status:** Load-bearing **reduction** of E(1), not a full proof of \(\lim\alpha_n\).

## Statement

For Seidel \(A,C\) differing in \(k\) undirected edges, for every \(x\in\{\pm1\}^n\),
\[
|Q_A(x)-Q_C(x)|\le 2k,
\qquad
\Phi(A)\ge\Phi(C)-2k.
\]
(Shipped: `edge_hamming`, `phi_edge_lipschitz_lower`.)

Degree form (Prop 15.20c): disagreement max-degree \(D\) gives \(\Phi(A)\ge\Phi(C)-Dn\).

## E(1) criterion (improved)

On \(n=p^2+1\) with Paley \(C\) (\(\rho=1\), \(\Phi(C)=\tfrac12 n\sqrt{n-1}\)), let
\[
k_\star
=\min\bigl\{d_H(A',C):
A'\text{ a Seidel switching of a }\Phi\text{-minimiser},\,
\Phi(A')=m_n\bigr\}.
\]
Edge lip yields
\[
m_n\ge\Phi(C)-2k_\star.
\]
Hence **E(1)** \(m_n=\Phi(C)-o(n^{3/2})\) holds as soon as
\[
k_\star=o(n^{3/2}).
\]
(The older Frobenius lip needed the stronger \(k_\star=o(n)\).)

## Certified checks

| \(n\) | object | \(k\) | bound \(\Phi-2k\) | actual \(\Phi\) |
|------:|--------|------:|------------------:|----------------:|
| 10 | matching undercutter of Paley | 5 | 5 | 13 |
| 10 | all 144 \(\Phi=13\) matchings | 5 | 5 | 13 |
| 26 | Paley itself | 0 | 65 | 65 |

Relative gap at \(n=10\): \(2/n^{3/2}\approx0.063\). Any sequence with \(k_\star=O(n)\) has relative gap \(O(n^{-1/2})\to0\).

## What remains for settlement

One sufficient route is to prove \(k_\star=o(n^{3/2})\) for a closest exact
\(\Phi\)-minimiser on all large \(n=p^2+1\) (or the weaker
\(k_\star=o(n^2)\) using Max-Lipschitz). This is **not equivalent** to saying
that every near-minimiser is Hamming-close: far exact minimisers already occur
at \(n=10\), and general far edge-local minima are proved in
`NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md`.

This is **not** proved. SA+exact at \(n=14,18,26\) found no undercut of Paley \(\Phi\)
(so those runs are consistent with \(k_\star=0\)), but absence of a search hit is not a proof.

**Existence of \(\lim\alpha_n\) remains OPEN** until this rigidity step (or another full argument) lands.

## Tests

`tests/test_minmax.py`:
- `test_phi_edge_lipschitz_pointwise_and_global`
- `test_phi_edge_lipschitz_matching_undercut_n10`
- `test_e1_edge_lipschitz_criterion_string`
