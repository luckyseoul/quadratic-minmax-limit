# A direct physical Mobius parallel design at the p=31 top endpoint

**Status:** exact constructive certificate for the balanced top parallel
ledger and for one clean nonorigin cancellation.  Sixteen localized Mobius
halves give a ternary antisymmetric support of 478 nonfixed inversion orbits;
one fixed antipodal edge completes a 479-edge graph whose parallel counts are
exactly the required hard `14^14 15^2` and opposite `15^3 16^13` profiles.
The transverse compact-atom cells and even moments are not matched here, so
this is not a common-graph construction and does not close residual (ii).

## Changed premise and duplication check

The earlier `residual_branch_c_aux_sdr_cpsat.py` selects complementary target
pairs and checks only distinct auxiliary directions and their Paley types.
The later sparse script imports the same paired proxy and imposes lower bounds
and an `l1` relaxation, but it does not impose the exact balanced 32-direction
parallel vector.  The present certificate instead replays the ordinary
parallel profile of every individual half from all thirty of its physical
edges.  Its final equality is therefore independent of the complementary
profile identity used to discover the initial assignment.

## Exact raw and corrected ledgers

Projective directions use the canonical order returned by
`projective_functionals(31)`.  The frozen sixteen actual scaled pairs `(L,M)`
are in `HALVES` in the source module.  Their target directions are exactly the
sixteen hard directions.  Their auxiliary projective directions are distinct,
with fourteen hard and two opposite.

Summing the thirty-edge half profiles gives

```text
raw =
(15,14,14,15,16,16,16,14,16,14,14,16,16,15,15,14,
 14,16,16,16,16,14,14,16,14,16,16,16,14,14,14,14).
```

Thus the hard multiset is already `14^14 15^2`, while the opposite multiset
is `15^2 16^14`.  Direction 5 is opposite and has raw count 16.  Removing the
two edges in one cancelled orbit and adjoining one fixed antipodal edge in
that same spatial direction makes its count 15.  The resulting exact vector
is

```text
final =
(15,14,14,15,16,15,16,14,16,14,14,16,16,15,15,14,
 14,16,16,16,16,14,14,16,14,16,16,16,14,14,14,14),
```

with the required opposite multiset `15^3 16^13`.

For each half the replay independently constructs its edges

```text
u_t = j(1,t/(t+1)),  v_t = j(t,t),  t != -1
```

in `(L,M)` coordinates, counts the unique spatial direction of every edge,
and checks the result against the closed `0/1/2` one-half formula.  Hence the
displayed equality is an edge-level identity, not a type-count proxy.

## Clean physical cancellation from center scaling

Scaling a half's nonzero center multiplies every endpoint by that scalar and
does not change its parallel profile.  Use, in the order of `HALVES`,

```text
centers = (27,9,28,6,19,21,28,9,19,16,9,3,4,1,30,13).
```

The exact pairwise replay then has one and only one shared inversion orbit:

```text
halves       = (2,13)
orbit        = {(2,25),(29,1)}
coefficients = (-1,+1)
direction    = 5
origin       = no.
```

All other half pairs are disjoint.  In particular there is no same-sign
overlap and no triple overlap.  Cancelling the displayed orbit leaves 478
nonfixed orbits.  Add the fixed edge

```text
{(5,30),(26,1)},
```

whose spatial direction is also 5.  Selecting the physical side prescribed
by each surviving normalized orbit coefficient gives 479 distinct graph
edges.  A fresh normalized edge--Radon calculation returns the final vector
above.  The graph-edge serialization has SHA-256

```text
c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d
```

and it contains all sixteen forced origin edges; the cancellation is not an
origin shortcut.

## Exact scope

This construction closes three previously separate feasibility questions for
one explicit nonzero center profile: auxiliary distinctness/type, ternarity
with the sharp support loss two, and every parallel quota.  It also refutes a
parallel-profile obstruction at the top endpoint.  It does **not** show that
the transverse row left after subtracting each literal star decomposes into
the prescribed compact atoms, and it does not impose the simultaneous global
degree-six/eight forms.  Those are the remaining common-graph conditions.

The limitation is already visible without an atom solver.  A hard row has at
most `P-3` compact three-edge atoms after its literal star is removed.  An
opposite row has six all-positive triangles and `Q-9` compact atoms, again
`Q-3` three-edge atoms in total.  Hence every row necessarily has residual
transverse `l1` norm at most `3(P-3)` or `3(Q-3)`.  Direct evaluation of this
frozen graph violates that necessary bound in all 32 directions: the smallest
excess is 122 and the largest is 194.  Thus the displayed graph is decisively
a support-plus-parallel skeleton, not a hidden full common graph.  This check
does not obstruct another choice among the many direct parallel designs.

For downstream exact work, `centered_physical_graph()` returns all 479
canonical unordered edges, the fixed edge, the hash, and every hard center in
the canonical projective-row scaling.  `transverse_compact_l1_diagnostic()`
replays the limitation above from that public graph record.

## Replay

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_direct_mobius_parallel_design.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_p31_direct_mobius_parallel_design.py
```
