# Four-point residual boundary exclusion for `p>=11`

Proposition 15.652 excludes every odd-degree boundary `D` of size four at
residual size `|H|=4p+1` for every odd prime `p>=11`.  It uses only the
type-split slack budget of Proposition 15.632 and finite affine-plane pair
counting.  There is no asymptotic step and no solver certificate in the
all-prime argument.

Together with the empty-boundary exclusion (15.632) and the complete
two-point exclusion (15.650--15.651), the first still-open boundary size is
at least six for `p>=11`.  The size-four cases at `p=5,7`, larger boundary
shapes, residual (ii), R1, and the limit remain open.

## 1. Exact floors for at most four odd fibres

Write `M(p,b,eta)` for the exact quadratic parity-majorant expectation from
15.632 and

```text
F(p,b,eta) = 2 ceil(p M(p,b,eta)).
```

For every odd `p>=7`, exact positive degree-two quadrature gives

| `b` | phase `0` | phase `1` |
|---:|---:|---:|
| 0 | `0` | `2p` |
| 1 or 2 | `p+1` | `p-1` |
| 3 or 4 | `2p-6` | `2p` |

Here are explicit optimal majorants.  For `b=1`, interpolate the two parity
values.  For `b=2`, the phase-zero and phase-one polynomials are
`-t^2+2t` and `(t-1)^2`.  For `b=3,4`, phase zero is `(t-2)^2` and phase one
is `(t-1)^2` at `b=3` and the constant `1` at `b=4`.

Optimality is certified by expressing expectation of every quadratic as a
nonnegative combination of its values at parity contacts.  The nontrivial
quadrature weights are

```text
b=3, phase 0, nodes 1,2,3:
  3(p-3)/(4p), 3/p, (p-3)/(4p)

b=3, phase 1, nodes 0,1,2:
  (p-3)/(4p), 0, 3(p+1)/(4p)

b=4, phase 0, nodes 1,2,3:
  (p-5)/(2p), 3/p, (p-1)/(2p)

b=4, phase 1, nodes 0,2,4:
  (p-7)/(8p), 3(p+1)/(4p), (p+1)/(8p).
```

They are nonnegative for `p>=7`, sum to one, and reproduce moments of
degrees zero, one, and two exactly.  Thus they reproduce the expectation of
every feasible quadratic.  At each positive-weight node the displayed
majorant meets the required parity value, proving the lower bound; direct
expectation proves equality.  The source companion checks all identities in
exact rational arithmetic.

## 2. Four finite boundary points

Suppose infinity is not in `D`.  In a direction `d`, the four finite points
have one of the following fibre partitions:

| partition | odd fibres `b_d` | colliding pairs `c_d` |
|---|---:|---:|
| `1+1+1+1` | 4 | 0 |
| `2+1+1` | 2 | 1 |
| `2+2` | 0 | 2 |
| `3+1` | 2 | 3 |
| `4` | 0 | 6 |

Every unordered pair determines exactly one projective direction, so

```text
sum_d c_d = C(4,2) = 6.                              (1)
```

All `b_d` are even.  The parity phase is one in the quadratic type
`epsilon_d=c_H` and zero in the opposite type.  Call the former the bad
type.  Put `m=(p+1)/2`; each type has `m` directions and exact budget
`m(p+1)`.

In the bad type, `b=0,4` costs `2p`, whereas `b=2` costs `p-1`.  If `x`
bad-type directions have `b=2`, the budget forces

```text
x(p-1) + (m-x)2p <= m(p+1),
hence x >= (p-1)/2.                                  (2)
```

Every `b=2` direction consumes at least one of the six pairs in (1).  This
immediately contradicts (2) for `p>=17`.

At `p=13`, all six pairs must be six distinct bad-type `b=2` directions.
The seven good-type directions consequently all have `b=4`, costing
`7(2p-6)=140` against budget `98`.

At `p=11`, five bad-type `b=2` directions leave at most one collision for
the good type.  At least five of its six directions therefore have `b=4`;
the smallest possible sixth cost is the `b=2` cost.  The resulting lower
bound is

```text
5(2p-6) + (p+1) = 92 > 72.
```

Thus four finite boundary points are impossible for every odd prime
`p>=11`.

## 3. Infinity plus three finite boundary points

Now suppose infinity is in `D`.  The finite partitions are

| partition | odd fibres `b_d` | colliding pairs `c_d` |
|---|---:|---:|
| `1+1+1` | 3 | 0 |
| `2+1` | 1 | 1 |
| `3` | 1 | 3 |

There are only three finite pairs, so at most three directions have `b=1`.
The parity phase is now independent of direction: it is zero for `c_H=+1`
and one for `c_H=-1`.

For `c_H=+1` and `p>=11`, `b=1` costs `p+1`, while `b=3` costs
`2p-6>p+1`.  A type budget equals the cost of making all its directions
`b=1`, so every one of the `p+1` directions would have to satisfy `b=1`.
Only three can.

For `c_H=-1`, `b=1` costs `p-1` and `b=3` costs `2p`.  Relative to making
all directions `b=1`, each `b=3` direction costs `p+1`; each type has room
for at most one.  Hence at least `p-1` total directions must have `b=1`,
again contradicting the upper bound three for every `p>=7`.

Both product signs are therefore impossible with infinity in `D` for
`p>=11`.  Combining Sections 2 and 3 proves the proposition.

## 4. Exact scope

- all size-four boundaries: closed for every odd prime `p>=11`;
- infinity plus three finite points with `c_H=-1`: also closed at `p=7`;
- `p=7`, infinity present, `c_H=+1`: still open in this proposition;
- other `p=5,7` size-four shapes: still open;
- size six and larger, residual (ii), R1, global QVAR, and `L`: open.

Reproduce the machine companion with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15652.py
pytest -q tests/test_prop15652.py
```
