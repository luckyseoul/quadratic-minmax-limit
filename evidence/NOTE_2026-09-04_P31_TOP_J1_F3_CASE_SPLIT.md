# Exact case split for the p31 top `j=1, f=3, d=0` localized-Mobius ledger

## Scope

This note classifies every fixed-direction parity support and every forced
auxiliary multiplicity pattern in the `p=31` top localized-Mobius branch with

\[
 j=1,\qquad f=3,\qquad d=0.
\]

It is conditional on the established top parallel ledger.  It does **not**
show that a survivor admits compatible relative scales, centers, physical
collision loci, compact atoms, or one common graph.  It therefore does not
close residual (ii).

## Ledger and notation

There are `m=16` localized halves and `30` raw physical edges per half.  The
final graph has `479` edges.  Three of these are fixed antipodal edges, so the
nonfixed support has size `476`; consequently the `480` raw half occurrences
lose exactly two cancellation units.

Let `P_D` be the final parallel profile, `R_D` the raw sixteen-half profile,
`f_D` the number of fixed edges in direction `D`, and `kappa_D` the number of
cancellation units in direction `D`.  Then

\[
 R_D=P_D+2\kappa_D-f_D,
 \qquad \sum_D f_D=3,
 \qquad \sum_D\kappa_D=2.                 \tag{1}
\]

Write `b_D=f_D mod 2`.  Since three fixed edges are distributed either as
`3`, `2+1`, or `1+1+1`, `wt(b)` is one or three.  Let

\[
 v_D=P_D+1_{D\text{ hard}}\pmod2,
 \qquad a=v+b.                              \tag{2}
\]

The auxiliary-direction occurrence parity is `a`.  At this endpoint `v` is
supported on fourteen hard directions and three opposite directions, hence
`wt(v)=17`.

For any fixed support `b`, put `k=|b intersect supp(v)|`.  Equation (2) gives

\[
 \operatorname{wt}(a)=17+\operatorname{wt}(b)-2k.       \tag{3}
\]

Because there are only sixteen auxiliary occurrences, (3) is already a
strong capacity restriction.

## Adaptive-selector parity

The adaptive kernel-selector product from the preceding boundary-parity note
gives

\[
 h+f+z=0\pmod2,                             \tag{4}
\]

where `h` is the number of hard auxiliary occurrences and `z` is the parity
of origin-containing cancellation units.  Here `f=3`, while (2) and the even
fourteen-element hard support of `v` give

\[
 z=1+|b\mathbin{\cap}\text{hard}|\pmod2.    \tag{5}
\]

If `wt(a)=16`, every direction in `supp(a)` occurs exactly once.  This is an
auxiliary SDR, so no origin collision is possible and (5) must be zero.

If `wt(a)=14`, write the multiplicities as `n_D=a_D+2r_D`.  Since their sum is
sixteen, `sum r_D=1`.  There is a unique repeated auxiliary direction `A`:

- if `A` is in `supp(a)`, it occurs three times and the other thirteen support
  directions occur once;
- otherwise, all fourteen support directions occur once and `A` occurs
  twice.

This one repeated group can account for at most one origin cancellation unit.
Thus (5) determines the exact origin count, not only its parity.

## Complete split

### `wt(b)=1`

Capacity forces the single direction `B` into `supp(v)`, giving `wt(a)=16`
and an auxiliary SDR.  Equation (5) then forces `B` to be hard.  Therefore:

- `B` in the fourteen hard directions of `supp(v)` survives the
  parity/support gate and requires two nonorigin cancellations;
- `B` in the three opposite directions of `supp(v)` is selector-impossible;
- `B` outside `supp(v)` is capacity-impossible.

For a surviving `B`, the actual fixed allocation is either `f_B=3`, or
`f_B=1,f_C=2` for one `C!=B`.  Formula (1), which depends on the actual `f`
and not just `b`, must still be solved.

### `wt(b)=3` and `k=2`

Here `wt(a)=16`, so again the auxiliaries form an SDR and both cancellations
are nonorigin.  Equation (5) says that `b` must contain an odd number of hard
directions.  Writing `H_v,O_v` for the hard and opposite parts of `supp(v)`,
and `H_out,O_out` for their complements, the surviving type patterns are

\[
 2O_v+H_{out},\qquad H_v+O_v+O_{out},\qquad 2H_v+H_{out}. \tag{6}
\]

There are respectively `6`, `546`, and `182` labelled supports, for `734` in
total.  The other `1,306` supports with `k=2` violate (5).

### `wt(b)=3` and `k=3`

Now `b` is a three-subset of `supp(v)`, `wt(a)=14`, and the unique repeated
auxiliary pattern above is forced.  If `r=|b intersect hard|`, then:

- `r` even (`r=0` or `2`) forces one origin and one nonorigin cancellation;
- `r` odd (`r=1` or `3`) forces two nonorigin cancellations.

The labelled support counts are

| `r` | support count | origin / nonorigin cancellations |
|---:|---:|---:|
| 0 | 1 | 1 / 1 |
| 1 | 42 | 0 / 2 |
| 2 | 273 | 1 / 1 |
| 3 | 364 | 0 / 2 |

In particular, if all three fixed directions are opposite, then `b` is
exactly the three opposite directions in `supp(v)`.  The auxiliary parity
support is exactly the fourteen hard directions, and the two cancellations
split as one origin plus one nonorigin.  This is a genuine surviving
mechanism, not a contradiction.

## Boundary-signature form of the remaining gate

Let `G` be the XOR of the sixteen center-independent half signatures.  A
fixed edge in direction `D` contributes `e_D`; hence the three fixed edges
contribute `b`.  An origin correction at the repeated auxiliary `A`
contributes `e_A`.  A nonorigin orbit pair with endpoint annihilator
directions `U,V` contributes

\[
 e_U+e_V,\qquad U\ne V.                     \tag{7}
\]

If its spatial parallel direction is `D`, necessarily `D` differs from both
`U` and `V`.  Therefore the exact remaining boundary equation is

\[
 G+b=\bigoplus_i C_i.                       \tag{8}
\]

Two nonorigin corrections make the right side have even weight in
`{0,2,4}`.  One origin plus one nonorigin correction makes it have odd weight
in `{1,3}`.  These are necessary tests to impose together with (1) and the
exact Mobius collision loci.  The aggregate product alone cannot sharpen
them further.

## Exhaustive algebra replay and negative computational facts

The deterministic module exhausts all `4,992` possible parity supports of
weight one or three and all `5,984` allocations of three indistinguishable
fixed edges.  It finds `1,428` supports and `1,862` actual fixed allocations
surviving this exact parity/support gate.

Two exploratory calculations are recorded only to prevent overreading:

- replacing the two opposite auxiliaries in the frozen `j=0` witness by one
  repeated direction gave zero exact-profile hits;
- a global 14,912-choice CP-SAT model for (1) and the forced multiplicities
  ran for 306 seconds and returned `UNKNOWN` without a feasible point.

The first is only a tiny local repair family.  `UNKNOWN` is neither SAT nor
UNSAT.  Neither observation is evidence that the global `f=3,d=0` family is
empty.

Replay:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_top_j1_f3_case_split.py
```

Source:

- `src/e1_gmin_m4_p31_top_j1_f3_case_split.py`
- `tests/test_p31_top_j1_f3_case_split.py`

The next unresolved implication is simultaneous satisfaction of the raw
profile equation (1), the boundary equation (8), the exact origin/nonorigin
collision loci, and all compact-atom rows.  Residual (ii), E1, `L=1/2`, and
the original MathOverflow limit remain **OPEN**.
