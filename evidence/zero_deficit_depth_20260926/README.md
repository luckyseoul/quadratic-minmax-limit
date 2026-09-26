# Zero-deficit transform-depth discriminator: first executions

2026-09-26. Evidence for the discriminator script committed as
`10fbea8` (`scripts/zero_deficit_transform_depth.py`), which had not
produced an output before this session. This directory records its first
completed runs, the accelerated re-implementation used to reach larger
orders, independent validation of the recursion semantics, and the
resulting data. Scope per `../SCOPE_FREEZE_2026-09-25_ONE_VERTEX_ONLY.md`
acceptance item 5: this is a discriminator/certificate computation for the
one-vertex extension route. It proves no all-orders statement.

## What is computed

For a signing A (order n, F=Phi(A)) and `r = F[(1+1/n)^(3/2)-1]` (the
neutral increment), the discriminator measures the minimax depth of the
pure ZERO-DEFICIT same-phase compatibility recursion of the Banaszczyk
one-vertex route: at each node (a principal block and phase), take all
projective global maximizers; choose a pivot; orient every maximizer to
pivot spin +1; children are all pairwise disagreement blocks, restricted
to genuine sigma-maximizers; recurse. The value is the worst branch length
under an optimal (depth-minimizing) pivot choice, with unit transform cost
per generation.

Verdict: `depth < r` means the pure zero-deficit same-phase branch cannot
exhaust an initial extension buffer r under an optimal elimination order
(the mechanism survives on the zero-deficit branch); `depth >= r` means it
can (the deficit-carrying branches are then required; unit-cost worst case).

## Results

All depths are exact integers from the recursion (not sampled).

| object | n | F | r_n | depth | margin r-depth | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| two-half exact minimizer | 6 | 5 | 1.301 | 1 | +0.301 | safe |
| two-half exact minimizer (m_7=9) | 7 | 9 | 1.996 | 1 | +0.996 | safe |
| two-half exact minimizer | 8 | 10 | 1.932 | 1 | +0.932 | safe |
| `A9_exact` (m_9=12) | 9 | 12 | 2.055 | 2 | +0.055 | safe |
| plus-I lift (`m_10=13` exact) | 10 | 13 | 1.998 | 2 | **-0.002** | **exhaustion** |
| `B13` block (m_13=20 exact) | 13 | 20 | 2.352 | 2 | +0.352 | safe |
| Paley conference | 14 | 21 | 2.290 | 2 | +0.290 | safe |
| `K15` exact minimizer | 15 | 27 | 2.745 | 2 | +0.745 | safe |
| `m16` witness | 16 | 30 | 2.856 | 2 | +0.856 | safe |
| Paley conference | 18 | 33 | 2.788 | 2 | +0.788 | safe |
| `m19` witness | 19 | 39 | 3.119 | 2 | +1.119 | safe |
| plus-I lift record (best known) | 25 | 60 | 3.636 | 3 | +0.636 | safe |
| Paley conference | 26 | 65 | 3.786 | 3 | +0.786 | safe |
| plus-I lift record (Phi=61) | 26 | 61 | 3.553 | 2 | +1.553 | safe |

Both phases are zero-deficit and give equal depths wherever both phases
attain F (the common case for the structured objects above); random
signings give trivial depth 0 (unique maximizer class) and are not
informative for the mechanism. The n=30 Paley conference has two
independent runs (reference script and parallel engine) that were still
in flight when this record was written and is therefore not claimed here.

## Findings

1. **First completed executions.** The committed discriminator now has
   outputs on 14 structured objects; before this session its only local
   artifact was a 0-byte file.

2. **The zero-deficit branch is safe at every tested structured object
   with n >= 13**, with margins +0.29 ... +1.12. Depth grows very slowly
   over the tested range (1 at n=6-8; 2 at n=9-19; 3 at n=25-26) while
   r_n grows like ~0.68 sqrt(n).

3. **The single exhaustion case is n=10** (the exact m_10=13 minimizer),
   by 0.002 in the unit-cost bound: depth 2 vs r=1.998. Since convergence
   needs only all sufficiently large n, this small-order failure is not an
   obstruction; it does show the property is not universal, and that
   deficit-carrying branches (or sub-unit pivot costs) are needed at small
   orders.

4. **The n=26 record object is the safest of all**: depth 2, margin
   +1.55 -- shallower than the Paley conference of the same order (depth
   3), so the better object is not the harder one for this branch.

4. **K15 twin.** The K15 discriminator was independently executed and
   recorded by the concurrent session as
   `../zero_deficit_transform_depth_k15_20260926.json` (commit `d99955b`).
   Its numbers agree exactly with `results/depth_K15.json` (depth 2, same
   witness chain 15 -> 8 -> 4): two independent executions.

5. **Triple-checked semantics.** The depth values from the committed
   reference script (which runs the recursion single-core) agree with
   (i) the packed-popcount parallel engine `scripts/zero_deficit_depth_fast.py`
   (subclasses the reference Engine; only the enumeration kernel differs)
   on all seven pinned witnesses, and (ii) an independent pure-Python
   brute-force implementation of the documented recursion on the critical
   small cases (n=9, 10, 13, 14, 15, 16).
## Limits

* Unit-cost, same-phase, zero-deficit restriction: the tool tests the
  extreme R=0 ancestry only. Branches carrying deficit (R>0), the other
  half of the Sep 25 target, are outside its scope.
* "Safe" is a positive statement about the zero-deficit branch only; it is
  not evidence for convergence beyond what the one-vertex program would
  need from this branch.
* The n=10 exhaustion margin is within the unit-cost worst case; actual
  pivot coefficients c_i, d_i <= 1 could make the true cost smaller.

## Reproduce

    python3 scripts/zero_deficit_depth_fast.py evidence/K15_exact_minimizer_20260911.json --workers 8
    python3 scripts/zero_deficit_depth_fast.py evidence/coherent_BI_lift_20260919/best_n25.json --workers 84
    python3 scripts/zero_deficit_depth_fast.py evidence/ns_port_n26_undercut_20260912/stage01.npz  # pending engine support

(the reference script is `scripts/zero_deficit_transform_depth.py` at the
pinned hash above; matrices built from the stored blocks as in `results/`).
The fast engine was re-validated against the CuPy-fallback revision
(`c5df76c`) on all seven pinned witnesses.

## Receipts

`results/*.json` are the raw discriminator outputs (14 objects), pinned in
`SHA256SUMS`. Script hashes:

    scripts/zero_deficit_transform_depth.py   8e06c760b936ef44c31c2560e51b994e0ac4c9ad8a8a4754c840769839705a21  (10fbea8 + c5df76c CuPy-fallback fix)
    scripts/zero_deficit_depth_fast.py        71e08a794d8a4d18653d95220edb472acabfcfed6007339288e95ccf74fc1e18

Environment: soulkiller, Python 3.14.4, numpy 2.4.4. The n=30 Paley
object was attempted with two engines and stopped unfinished (duplicated
compute; recorded here so it is not mistaken for a pending verification).
