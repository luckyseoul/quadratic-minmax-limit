# Matched-block gluing: exact minima, both-sided bounds, and what it means for the doubling route

**Status:** 2026-09-13 session note. Proved: Proposition 1 (both-sided bound),
Proposition 2 (lower-bound construction). Exact data: full enumerations of
three matched-block gluing families. No convergence claim; scoped to the
two-block (doubling) route.

**Convention (fixed once and used throughout).** For a zero-diagonal
antisymmetric `K` of order `N`,

    Phi(K) = max_{s in {+-1}^N} | sum_{i<j} K_ij s_i s_j |,

with `K_ij` the upper-triangle reading, evaluated over the **full cube**
(2^N states). A complementary-pair argument may reduce the state space by
one bit only for a single quadratic form; it must not be applied to split
sums whose two terms are not separately symmetric. All numbers below were
computed with a literal full-cube kernel (`scripts/glue_exact_raw.py`),
validated against the recorded table `m_2..m_6 = 1,3,4,4,5` and against an
independent literal double loop on a random merged matrix.

## 1. Identity and both-sided bounds

Let `A` (order `n`), `C` (order `m`) be signings, `eps in {+-1}^{n x m}` the
cross block, and `phi_A = max|Q_A|`, `phi_C = max|Q_C|` with
`Q_A(x) = sum_{i<j<=n} A_ij x_i x_j`. The merged signing
`K = [[A,eps],[-eps^T,C]]` satisfies the exact identity

    Phi(K) = max_{x,y} | Q_A(x) + Q_C(y) + x^T eps y |.        (1)

**Proposition 1 (bounds).** For every `eps`,

    phi_A + phi_C - nm  <=  Phi(K)  <=  phi_A + phi_C + nm.    (2)

*Proof.* Upper: triangle inequality over the three terms.
Lower: at a pair maximizing `|Q_A + Q_C|` the internal part equals
`+-(phi_A+phi_C)`; the cross part is at least `-nm` in absolute value
contribution. QED

**Proposition 2 (the infimum over `eps` is not the independent-budget
combination).** `min_eps Phi(K) <= phi_A + phi_C`, because for `eps`
minimizing the cross term at one fixed pair the cross contribution is
`<= 0`. In particular the recorded cross-term floor
`max_{x,y}|x^T B y| >= (sqrt(2/pi)-o(1)) N^{3/2}` constrains objects
maximized over `B`; it does **not** lower-bound `min_eps Phi(K)`. The two
functional operations (sup over `x,y` of a sum, versus inf over `eps`)
cannot be interchanged, and the independent-budget reading of the recorded
floor does not apply to the gluing minimax.

## 2. Exact minima (full enumeration, one validated kernel)

Blocks are exact optima of their own order (the `m_n` in the table below are
the recorded/independently reproduced exact minima; every declared
completion was scored).

| case | completions | min over `eps` | recorded `m_N` | ratio |
| --- | --- | --- | --- | --- |
| 4+4 (N=8) | 2^16 (all) | **10** | 10 | 1.00 |
| 4+5 (N=9) | 2^20 (all) | **12** | 12 | 1.00 |
| 5+5 (N=10) | 2^25 (all) | **13** | 13 | 1.00 |

Counts of optimal completions (Phi equal to the family minimum) and the
next values in the exact histogram:

| case | optimal | Phi+2 | Phi+4 | Phi+6 |
| --- | --- | --- | --- | --- |
| 4+4 | 184 | 12128 | 30768 | 17044 |
| 4+5 | 10800 | 174080 | 386400 | 320880 |
| 5+5 | 1840 | 54320 | 2935320 | 12768400 |

So at orders 8, 9, 10 the matched-block gluing family **attains** the
recorded optimum of the merged order: the doubling target is met with
equality, not merely `<= 2 sqrt(2) m_n`. This is a finite statement at three
orders.

## 3. What the doubling proof still needs

Equation (1) shows the task splits into two joint conditions:

1. the internal part `Q_A + Q_C` must be small on every pair where the cross
   part is large (and vice versa) -- the maxima of the two terms cannot be
   independent because they are evaluated at the same state pair;
2. the cross block must realize a negative correction on the pairs that
   maximize `|Q_A + Q_C|`.

The exact witnesses at N=8,9,10 show both are satisfiable at these orders.
The recorded order-16 completion runs (fixed `A8` sources; family floor 32
against target `2 sqrt(2) m_8 ~= 28.3`) show a recorded family that misses
the target there -- but that family fixes the source, while the minima above
allow the blocks to be chosen. The open content of the doubling estimate is
therefore a joint alignment/negativity statement at all orders with a
Dini-summable error; it is not a cross-term margin (Proposition 2).

## 4. Scope, corrections, reproducibility

- The exact minima here were recomputed after a convention error in an
  intermediate kernel (a half-cube state reduction applied where the split
  breaks the complementary symmetry). The full-cube kernel used for the
  table above reproduces: raw value 13 for the 5+5 witness, full-cube
  maximum 13, half-cube maximum 13 (equal here), and the literal loop 13.
- `Phi` values quoted for 4+4, 4+5, 5+5 match the independent literal
  double-loop evaluation on the assembled matrices.
- These are finite statements; no all-orders implication is claimed.

Reproduce:

    python3 scripts/glue_exact_raw.py --n 4 --m 4 --workers 64 --chunk-bits 12 \
        --out evidence/merge_obstruction_20260913/raw_glue_4p4.json
    python3 scripts/glue_exact_raw.py --n 4 --m 5 --workers 64 --chunk-bits 12 \
        --out evidence/merge_obstruction_20260913/raw_glue_4p5.json
    python3 scripts/glue_exact_raw.py --n 5 --m 5 --workers 64 --chunk-bits 12 \
        --out evidence/merge_obstruction_20260913/raw_glue_5p5.json

Artifacts: `raw_glue_4p4.json`, `raw_glue_4p5.json`, `raw_glue_5p5.json`
(full histograms and best completion inside each). Superseded: the earlier
`glue_*.json` / `glue_hist_*.json` files in this directory, retained only as
the audit trail of the convention error.

## 5. Addendum (same session, later): N=12 also attains the optimum

A follow-up search with the exact full-cube kernel (V100, `glue_exact_raw.py`
imported on soulkiller) found a 6+6 cross completion with

    Phi = 18 = m_12 (exact value re-verified by the same kernel),

stored as `witness_6p6.json` (36 cross entries, `verified_exact: true`).
An earlier, weaker probe in this session reported 22 for the same family;
that was a search artifact, not a lower bound -- the witness above
supersedes it. The family minima now read:

| case | family min | recorded `m_N` | source |
| --- | --- | --- | --- |
| 4+4 | 10 | 10 | full enumeration 2^16 |
| 4+5 | 12 | 12 | full enumeration 2^20 |
| 5+5 | 13 | 13 | full enumeration 2^25 |
| 6+6 | 18 | 18 | explicit witness, exact re-evaluation |

So the matched-block gluing family attains the optimum at N = 8, 9, 10
and 12. (The 6+6 family minimum could in principle be lower than 18 by the
both-sided bound; 18 is an attained value. The recorded `m_12 = 18`
bounds it from below, so 18 is exactly the family optimum.)

A full 2^36 enumeration of the 6+6 family was launched on the V100
(`scripts/glue_gpu_scan.py`, shardable) but not completed within this
session's budget; it can only confirm 18 and quantify the count of 18-level
completions. No convergence claim.
