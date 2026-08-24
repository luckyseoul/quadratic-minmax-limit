# R1 and residual-(ii) diagnostics after the profile-glue theorem

Date: 2026-08-24.  These are killed routes and finite diagnostics, not a
numbered proposition and not a leftover close.

## The minimum-shell graph is not an association-scheme shortcut

Join two Max+ vectors when their Hamming distance is the lattice minimum
(p+1), equivalently when their signed inner product is maximal below the
antipode.  At (p=5) the signed graph has degree 13 and the projectivized graph
has degree 26.  At (p=7), however, the projective degrees are

\[
0\ (2352\text{ vertices}),\qquad
22\ (8400),\qquad
50\ (700).
\]

Thus the first-shell graph is already nonregular at (p=7).  A one-adjacency
Poincare or strongly-regular-graph proof of R1 cannot be justified by the
minimum distance alone.

## Real functions constant on one Max-minus slice

For edge features (f_{ab}(y)=C_{ab}y_ay_b), fix
(U=\{y\in\operatorname{Max-}:f_{01}(y)=-1\}).  The numerical function-space
ranks are:

| p | edge features | full score rank | functions constant on U | explicit wedge rank |
|---:|---:|---:|---:|---:|
| 5 | 325 | 66 | 14 | 13 |
| 7 | 1225 | 276 | 25 | 25 |

At (p=7), every score function constant on (U) is generated, modulo the
global feature kernel, by one constant star, (f_{01}), and the wedges

\[
f_{0k}+\Delta_k f_{1k}\qquad(k\ge2).
\]

This is a floating-point rank diagnostic (`scripts/leftover_real_slice_span.py`),
not an all-prime rank theorem.

More importantly, function-space rigidity does not imply sparse graph
rigidity.  With (|G|=4p), (01\notin G), and (S_G=-2) on all of (U), exact
MILPs maximize the number of edges disjoint from ({0,1}) at

\[
14\quad(p=5),\qquad 18\quad(p=7).
\]

Both optima were proved in one branch-and-bound node.  Hence the hoped-for
literal double-star conclusion is false at both primes.

## Full official p=7 model

The official (k=4p) linear conditions were imposed on all deduplicated
Max+/Max- rows:

\[
S_G\ge2\ (\operatorname{Max+}),\quad
S_G\le-2\ (U),\quad
S_G\le-4\ (\operatorname{Max-}\setminus U),\quad
|G|=28,\quad01\notin G.
\]

HiGHS reached 603.614 seconds with no incumbent and no infeasibility
certificate.  The continuous relaxation is feasible: it has 538 fractional
variables, total mass 28, and maximum coordinate 0.3015749569.  Therefore
there is no LP/Farkas close behind the finite model; any valid obstruction
must use integrality or additional nonlinear structure.  The timeout is
UNKNOWN and is not evidence of emptiness.

Reproduction: `scripts/leftover_minus_slice_milp.py`.  Machine-readable
summary: `evidence/r1_residual_diagnostics_2026-08-24.json`.
