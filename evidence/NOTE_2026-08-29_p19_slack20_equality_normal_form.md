# p=19 slack-20 equality normal form

**Date:** 2026-08-29
**Proposition:** 15.694
**Status:** proved; strict structural reduction of the four slack-20 rows,
not endpoint closure

Proposition 15.693 forces every slack-20 witness to use exactly five repair
deletions. Write

```text
S = A disjoint_union D,   |A|=11,   |D|=5,
```

where `A` is an affine arc. If a point of `D` lay on no `A`-secant, it could
be restored and would give a forbidden four-deletion repair. Thus every
deleted point has positive secant multiplicity. The standard line charge
and total slack give

```text
20 = slack(S) >= 4 sum_{x in D} mu_A(x) >= 4*5 = 20.
```

Every inequality is equality. In particular, all five deleted points have
`mu_A(x)=1`.

For an affine line `l`, put `a=|A cap l|` and `d=|D cap l|`. Since `A` is an
arc, `a<=2`. Comparing the exact line slack with its charged contribution
leaves only

```text
(a,d) = (0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(2,2).
```

Consequently `D` is itself an affine 5-arc, every boundary line has at most
four points, and a line through two deleted points contains either zero or
two core points. The five charged deleted/secant incidences have exactly
three possible global patterns:

```text
five 3-lines;
one 4-line and three 3-lines;
two 4-lines and one 3-line.
```

Let the profile have `t` undetermined directions and adjoin any two of their
infinity points to `A`. The result is a 13-arc `K`. Every deleted point has
secant index one relative to `K`, and each of the other `t-2` undetermined
infinity points has the line at infinity as its unique `K`-secant. Hence

```text
c1(K) >= 5+(t-2) = 7 for t=4, and 8 for t=5.
```

Al-Zangana's exhaustive classification has 2,733 projective 13-arc classes
(501 incomplete and 2,232 complete) and `c1<=9`. This sharply filters the
classes that can occur but does not by itself contradict either value 7 or
8. The exact p=19 remainder therefore remains seven profiles.

Two lossless exact diagnostics now encode this normal form:

- `scripts/p19_second_boundary_profile_cryptominisat.py`
- `scripts/p19_slack20_repair_cpsat.py`

Bounded trials returned `UNKNOWN`; those timeouts are not evidence and are
not used in Proposition 15.694.

Primary classification input: E. B. Al-Zangana, *The Geometry of the Plane
of Order Nineteen and its Application to Error-Correcting Codes*, PhD thesis,
University of Sussex, 2011, Chapter 4, Section 4.21, pp. 103--104.

Reproduction:

- `src/e1_gmin_m4_prop15694.py`
- `evidence/e1_gmin_m4_prop15694.json`
- `tests/test_prop15694.py`
