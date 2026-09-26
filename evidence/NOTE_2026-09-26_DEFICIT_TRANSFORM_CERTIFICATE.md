# Deficit-aware transform sections with actual sign recovery

This is a bounded constructive discriminator for the one-vertex target, not
a convergence result. Implementation and regressions are prepared; the
regressions and research runs have **not been executed** at this checkpoint.
In particular there is no new extension row, finite record, or summable-error
theorem to report.

The September 26 handoff explicitly names deficit-carrying ancestry as the
next discriminator. The earlier zero-deficit same-phase depth statistic is
retained unchanged. The new program includes every stable deficit and both
phases, and returns an actual row when its sufficient construction succeeds.

## Reused results and the exact change

The stable-state reduction is
`NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md`. The transform and its
pair-compatibility formula are already in
`NOTE_2026-09-25_BANASZCZYK_STABLE_RECURSION.md`; they are not new claims.
The affine reserve/cost bookkeeping extends
`NOTE_2026-09-25_FIELD_WIDTH_POTENTIAL_RECURSION.md`.

The important implementation distinction is **section versus projection**.
The pair constraints alone describe a projection of the transformed body.
For recovery of an incident sign row, take its zero-coordinate section and
also retain the expanded original slabs. Those constraints cannot simply be
dropped. This supplies a constructive sufficient test for all the original
stable constraints, including positive deficits and mixed-phase ancestry.

Deduplication searches covered the active notes/code, all three linked
worktrees, preserved branch heads, and the text/code members of the four
September 6 campaign archives. The older weighted-anchor Banaszczyk rounding
concerns a different cross-block construction. The September 25 transform
formula and September 26 zero-deficit engines are the directly reused work.
No stronger novelty claim about convex-body transforms is made.

## Exact section formula

Let `K` be the intersection of finitely many open symmetric slabs
`|c.z| < b_c`, all with positive half-width. Fix coordinate `i`, and suppose
`b_c > |c_i|` for every slab. This is exactly the condition `e_i in K`.
Define

```
T_i K = ((K-e_i) intersect (K+e_i)) + (-2,2)e_i.
```

Identify `z_i=0` with the remaining coordinate space. The section
`(T_i K) intersect {z_i=0}` has the following exact description:

1. Every original slab is retained, with coordinate `i` deleted and width
   `b_c+|c_i|`. A zero-pivot slab therefore persists unchanged.
2. Orient each nonzero-pivot normal to have positive pivot. For each pair
   `c,d`, with `p=c_i>0` and `q=d_i>0`, add

   ```
   |(q c-p d).z| < q b_c+p b_d-2pq.
   ```

To prove this, fix the horizontal point `u`. Its fiber in `K` is an open
interval `I=(L,U)`, provided all zero-pivot constraints hold. The eroded
fiber is `(L+1,U-1)`. Its dilation by `(-2,2)` contains zero exactly when

```
U-L > 2,     L < 1,     U > -1.
```

The first condition gives the pair constraints (and the stipulated
self-pair conditions `b_c>|c_i|`); the last two give the retained slabs.
The same argument covers an unrestricted fiber. Thus the formula is finite
and exact, without a Gaussian approximation or a coefficient-size assumption.

These three interval inequalities also imply that `I` contains at least
one of `-1,+1`: if both were excluded while `L<1` and `U>-1`, its length
would be at most two. Consequently any point of the section can be lifted
to a point of `K` by choosing one of the two actual coordinate signs.

Repeat sections along distinct coordinates. If every chosen pivot remains
strictly feasible, the final zero-dimensional body contains zero. Recover
the signs in reverse order, checking each of the two possibilities against
the saved parent slabs. This proves the sufficient construction implemented
in `scripts/deficit_transform_certificate.py`.

## Deficit and cost accounting

For an integer desired extension norm `T>=F=Phi(A)`, put
`r=T-F+1/2`. Each stable state of energy `E` starts with

```
c=x,     b=r+(F-E).
```

Because extension scores are integers, the open inequality
`|a.x|<T-E+1/2` is exactly the requested `E+|a.x|<=T` at a Boolean row.
The half-unit is solely an exact way to encode an integer target.

Store each width as `b=M r+R-C`. Initially `(M,R,C)=(1,F-E,0)`.
The two section operations have different updates:

```
retained slab: M'=M, R'=R, C'=C-|c_i|;
pair child:   M'=q M_c+p M_d,
              R'=q R_c+p R_d,
              C'=q C_c+p C_d+2pq.
```

Divide all three scalars by the positive gcd when making a normal
primitive. Its sign may be reversed freely. These identities retain
positive deficits and the mixed-phase origin of every selected facet.
`C` is a **net** cost and can be negative, because retained slabs expand.
The September 25 upper field potential is not asserted for these retained
section slabs.

Slabs with the same primitive normal may be reduced to the tightest width
**at the current target only**. A slab may also be dropped if
`b>sum_j |c_j|`: it is then strict everywhere on the Boolean cube. This
second simplification need not preserve the continuous body, but preserves
the recovery proof. Every recovered parent point is in the cube, where a
dropped slab holds automatically. Equality is not sufficient for this
pruning rule.

## What an output establishes

- `ROW_FOUND`: reverse recovery and a check of every original stable
  constraint succeeded. If the program exhaustively enumerated the stable
  skeleton in this run, its reported extension norm is exact by the existing
  stable-state theorem. A supplied cache instead produces a certificate
  explicitly conditional on that cache's completeness; checking its listed
  states is not a completeness proof.
- `TRANSFORM_STOPPED`: the selected transform path has no allowed next
  pivot, or the prescribed next pivot fails. This does **not** exclude an
  extension row or another construction. A regression deliberately includes
  a stopped path whose original body nevertheless contains a Boolean row.
- `RESOURCE_LIMIT`: the configured pair or facet budget was reached.
  This has no mathematical infeasibility meaning. A state-enumeration limit
  is reported before enumeration starts.

The greedy pivot choice is deterministic but is not a minimax theorem.
The program establishes neither global minimality of its source nor success
at every order. The exact neutral target uses integer square roots and the
extension parity; it does not round a floating-point value near a boundary.

## External verification and current scope

Per the handoff's compute instruction, no mathematical tests or source runs
were launched on the controller or remote hosts. The first external command,
from a checkout containing this change, is:

```sh
python3 -m unittest discover -s tests -p test_deficit_transform_certificate.py -v
```

It checks exact section geometry against an independent interval
description, positive-deficit/mixed-phase propagation, strict boundaries,
sign recovery, incomplete-cache labeling, resource stops, Gray-code state
enumeration, and a full-cube check on one tiny synthetic fixture. It does
not replay recorded research objects or the old depth engines.

After that gate passes, a research invocation has the form:

```sh
python3 scripts/deficit_transform_certificate.py --matrix INPUT.json \
  --enumerate-stable --excess 1 --save-skeleton NEW-skeleton.json \
  --output NEW-result.json
```

The new files are created exclusively; earlier outputs are never
overwritten. Subsequent target/order experiments should use `--skeleton`
to reuse the stored enumeration, retaining its explicit completeness scope.
Before selecting a research input, fetch and check concurrent results,
stored matrices and live processes under the September 26 handoff's
duplication rule. No new research input or result is claimed in this note.

The convergence-producing implication is still to prove successful row
construction with a uniform summable excess for suitable exact minimizers.
This prepared discriminator has not established that implication. Review
at this checkpoint is an analytic/code author review plus syntax checks;
the external regression result is pending. It is not a new mathematical
milestone and does not require a major-milestone backup.
