# Six-point residual boundary exclusion for `p>=11`

Date: 2026-08-26. This is Proposition 15.657. It combines the exact
type-split affine slack budget of Proposition 15.632, the small parity-floor
quadrature of Proposition 15.652, and a finite-plane pair-incidence bound.
It excludes every residual boundary of size six for every odd prime
`p>=11`. It does not close the exceptional primes `p=5,7`, larger boundary
sizes, residual (ii), R1, or the limit.

## 1. Pair-deficit budget

Let `s` be the number of finite points in the odd-degree boundary. In a
projective direction `d`, let the nonzero fibre multiplicities be
`n_1,n_2,...`, and let `b_d` be the number of odd multiplicities. Then

\[
 s-b_d=2\sum_i\lfloor n_i/2\rfloor
       \le 2\sum_i\binom{n_i}{2}.
\]

Every unordered pair of distinct affine points lies in exactly one
projective direction. Summing over all `p+1` directions gives

\[
 \boxed{\sum_d(s-b_d)\le 2\binom{s}{2}=s(s-1).}
\]

For six finite boundary points this budget is 30. If infinity is one of the
six boundary vertices, there are five finite points and the budget is 20.
The source checks the local inequality for every integer partition of five
and six.

## 2. Exact floors through six odd fibres

For `b=0,...,4`, Proposition 15.652 already gives the exact phase-zero and
phase-one directional floors

| `b` | phase 0 | phase 1 |
|---:|---:|---:|
| 0 | 0 | `2p` |
| 1,2 | `p+1` | `p-1` |
| 3,4 | `2p-6` | `2p` |

For `b=5,6` and odd `p>=11`, exact positive degree-two quadrature gives

| `b` | phase 0 | phase 1, `p=11` | phase 1, `p=13` | phase 1, `p>=17` |
|---:|---:|---:|---:|---:|
| 5,6 | `2p` | `2p-4` | `2p-2` | `2p` |

For phase zero, the constant candidate `q=1` has a moment-matching positive
quadrature on nodes `1,3,5`. For phase one and `p<=15`, the candidate
`q(t)=(t-3)^2` has expectation `3(p-5)/(2p)` and a positive quadrature on
`2,3,4`. For `p>=15`, `q=1` has a positive quadrature on `0,2,4`. The JSON
artifact records all weights and verifies moment matching, contact,
majorization, and agreement with the exact LP implementation over
`Fraction`.

## 3. Infinity plus five finite points

Here `b_d` is odd, and the residual parity phase is independent of direction
type: it is one exactly when `c_H=-1`. We may use the weaker phase-independent
floors

\[
 f(1)\ge p-1,\qquad f(3)\ge2p-6,\qquad f(5)\ge2p-4.
\]

Start from `p+1` directions of cost `2p-4`. Replacing `b=5` by `b=3`
saves at most 2 at deficit 2; replacing it by `b=1` saves at most `p-3`
at deficit 4. For `p>=11`, the latter has the larger saving per deficit.
The total deficit budget 20 therefore saves at most `5(p-3)`. The required
cost is at least

\[
 (p+1)(2p-4)-5(p-3)=2p^2-7p+11.
\]

Subtracting the total slack budget `(p+1)^2` leaves

\[
 p^2-9p+10>0\qquad(p\ge11).
\]

## 4. Six finite points for `p>=13`

Now `b_d` is even. Uniformly over both phases and every odd `p>=13`,

\[
 f(0)\ge0,\quad f(2)\ge p-1,\quad
 f(4)\ge2p-6,\quad f(6)\ge2p-2.
\]

Start from cost `2p-2` in every direction. The savings for deficits
`2,4,6` are at most `4,p-1,2p-2`, respectively. For `p>=13`, the last has
the largest saving per deficit, namely `(p-1)/3`. The pair budget 30 can
therefore save at most `10(p-1)`. Hence the total cost is at least

\[
 (p+1)(2p-2)-10(p-1)=2p^2-10p+8.
\]

Its excess over `(p+1)^2` is

\[
 p^2-12p+7>0\qquad(p\ge13),
\]

starting with gap 20 at `p=13` and increasing thereafter.

## 5. Six finite points for `p=11`

There are six directions of each quadratic type and each type has budget
72. The two types have opposite phases. At phase one, the floors at
`b=0,2,4,6` are `22,10,22,18`. Starting from six `b=6` directions costs
108; only a `b=2` replacement lowers the cost, by 8 at deficit 4. Reaching
72 therefore requires at least five such replacements and deficit at least
20.

At phase zero the floors are `0,12,16,22`. Starting from cost 132 requires
saving at least 60. The best saving per deficit is the `b=0` replacement,
`22/6`; because all deficits are even, this forces deficit at least 18.
The two types therefore require total deficit at least `20+18=38`, beyond
the geometric budget 30. The source independently enumerates the tiny type
count-profile LP and obtains the same minima.

## 6. Scope and literature check

Thus every six-vertex boundary is excluded for every odd prime `p>=11`.
The result is structural and produces no proposed integer sequence, so no
OEIS submission search is relevant. Targeted searches of the Johnson-slice,
finite-incidence, and Paley degree-parity literature found nearby tools but
no statement combining the affine slack floors with the pair-deficit and
quadratic-type budgets. Li--Zhou (arXiv:2512.19312) studies parity-uniform
*induced* Paley subgraphs and associated MDS codes; it does not study these
edge-set boundaries, directional Johnson slacks, or this exclusion. This is
a duplicate/context search record, not an unqualified priority claim.
