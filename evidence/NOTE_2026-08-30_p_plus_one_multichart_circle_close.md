# First `p+1` shell: signed multi-chart reduction and circle close

Date: 2026-08-30

Status: Propositions 15.722--15.724 are proved for every odd prime
`p>=17`. They exclude outside slack zero and every positive slack
`R<=max(3,floor(sqrt(p)-5/2))` at the first general boundary `|D|=p+1`,
and repair the middle floor-plus-two rule. Only slack beyond that cutoff
remains open.

## Exact phase transport

For a signed PSL map `g`, let `delta_g(v)` be its Paley endpoint
multiplier. Transporting a relative flip set leaves one multiplier at every
odd-degree vertex:

```text
c_(gH) = c_H product_(v in boundary(H)) delta_g(v).
```

For `g(z)=(az+b)/(cz+d)`, away from a pole this is `chi(cx+d)`.  If
`c!=0`, the pole and infinity both have multiplier `chi(c)`; if `c=0`,
there is no pole and every finite point and infinity have multiplier
`chi(d)=chi(a)`.  The affine case must not be obtained by substituting
`c=0` into `chi(c)`.

For an all-finite boundary with root polynomial `f`, inversion about `r`
therefore gives `c_r=c_H chi(f'(r))` on the boundary and
`c_r=c_H chi(f(r))` outside it. The product of all derivative characters is
one by the Vandermonde identity.

## Outside pair slack

For `P=p+1` affine points, the exact normalized slack is

```text
R = (sum_d b_d-P)/4
  = sum_(n=2r) r(r-1) + sum_(n=2r+1) r^2.
```

`R=1` would be one trisecant and no other line of occupancy at least three.
Deleting a point of that trisecant gives a `p`-arc on a conic. The deleted
off-conic point retains more than one conic secant after the missing conic
point is removed, contradicting uniqueness.

The input here is Segre's odd-order `q`-arc theorem for a `q`-arc in
`PG(2,q)`, not merely the classification of `(q+1)`-arcs.

At `R=2`, the rich lines are one 4-secant or two 3-secants. Deleting one
or two suitable points leaves a `p`- or `(p-1)`-arc. At `R=3`, deleting
three private rich-line points leaves a `(p-2)`-arc. Ball--Lavrauw's
classifications of complete `(q-1)`- and `(q-2)`-arcs extend these to a
`p`-arc and hence a conic for prime `p>=17`. Every restored point is
off-conic and retains respectively at least `(p-3)/2`, `(p-5)/2`, or
`(p-7)/2` conic secants, more rich lines than its total slack allows. Thus
`R=2,3` are impossible.

More generally, delete `n_l-2` points from each rich line and make the
deletion set `T` inclusion-minimal. Then `A=D\T` is an arc and

```text
1 <= |T| <= sum_l(n_l-2) <= R.
```

Minimality makes every deleted point lie on a line retaining two points of
`A`. If `R<=floor(sqrt(p)-5/2)`, the prime-field conic threshold of
Ball--Lavrauw puts `A` on a conic. A deleted point is off that conic and
retains at least `(p-1)/2-|T| > R` conic secants, a contradiction. Combining
the two arguments excludes every positive
`R<=max(3,floor(sqrt(p)-5/2))` without a finite search.

At `R=0`, Segre gives a conic with external infinity line. Its profile is
`m*b=0+m*b=2`. The phase-one budget permits at most one `b=0` direction, so
the conic direction character and Paley norm character disagree at at most
two projective points. A nonproportional pair would give a genus-one
character sum of magnitude at most `2 sqrt(p)`, while the agreement gives at
least `p-3`. Thus the conic is a Miquelian circle with exact type alignment.

## Full-circle phase

Normalize a circle point to infinity:

```text
D = {infinity} union (a+b F_p),  eps=chi(b),  m=(p+1)/2.
```

Sending an outside point to infinity multiplies the boundary sign by
`(-1)^m eps`. The resulting affine circle has its `b=2` tangent directions
in type `eps`. Alignment therefore forces `c_H=(-1)^m` in every circle-point
chart. This statement is independent of the two PSL circle orbits.

## Isolated-vertex contradiction

The complement of the circle has `p^2-p` vertices. A graph with `4p+1`
edges has at most `8p+2` nonisolated vertices, and all `p+1` boundary
vertices are nonisolated. Hence at least `p^2-8p-1>0` outside vertices are
isolated.

Send one isolated outside vertex to infinity. The chart has

```text
I=0,
m phase-zero b=0 directions,
m phase-one b=2 directions.
```

The phase-one floor and common-residue equation force `m-1` xnor baselines
and one elevation, with parallel counts `x,...,x,x+1`. Write the phase-zero
counts as `y+k_d`. Exact finite-edge counting gives

```text
m(x+y-7)=u-4,
```

so `u=4` and `x+y=7`. The xnor baseline congruence
`(p-1)/2 | I+x-4` then gives `x=4,y=3`. At least four phase-zero directions
have quotient zero and scaled mean eight. Since `b=0`, their slack is
`A=2B` with `B` a nonzero nonnegative integral quadratic, so
`4p E[B]=8`. Proposition 15.688 requires `4p E[B]>=p-3>=14`. This excludes
the full circle and completes the `R=0` close.

The reused two-coordinate calculation in 15.673 has both XOR and XNOR
signs.  Its sign parameter drops out of `(p-1)c=I+P_d-4`, so the divisibility
above is genuinely sign-independent; the two baselines are not identified.

## Floor-plus-two repair

On a Boolean cube, a nonnegative integral quadratic with parity
`1+sum_(i in R) z_i` has mean at least one for `|R|>=3` and at least `3/2`
for `|R|>=5`. The strict step follows from the degree-two Fourier cutoff and
half-integrality. Combined with the exact Johnson paired-cube operator,

```text
T A(X) = (A(X)+p E[A])/(p+1),
```

this excludes every middle odd-profile cell with scaled mean `2p+2`, except

```text
(p,b,phase)=(17,5,1), (17,11,0).
```

Both exceptions are real: `A=(t-3)^2` has `E[A]=18/17` and scaled mean 36.
They must remain in all later profile ledgers.

In the phase-zero far-contact step, the minimum active-coordinate counts are
`k-1` for odd `k` and `k-3` for even `k`, not `k` and `k-2`.  They are still
at least five throughout the range where the `3/2` cube bound is used.

## Duplicate and literature check

The repository search found the xnor coefficient congruence already in
15.672/15.673 and the sharp integral lift floor already in 15.688. The new
step is their isolated-circle synthesis after exact signed phase transport;
no earlier proposition makes that combination. Searches of the cited
MathOverflow problem, GitHub, and the nearby Paley/Miquelian/Johnson-slice
literature found general background but no duplicate of this residual
boundary argument. No new numerical sequence is asserted or submitted to
OEIS.

The low-slack extension imports S. Ball and M. Lavrauw, *Planar arcs*,
published Theorem 5 and Corollaries 10--11 (Theorem 3 and Corollaries 8--9
in arXiv v4).

## Artifacts

- `src/e1_gmin_m4_prop15722.py`, `tests/test_prop15722.py`,
  `evidence/e1_gmin_m4_prop15722.json`;
- `src/e1_gmin_m4_prop15723.py`, `tests/test_prop15723.py`,
  `evidence/e1_gmin_m4_prop15723.json`;
- `src/e1_gmin_m4_prop15724.py`, `tests/test_prop15724.py`,
  `evidence/e1_gmin_m4_prop15724.json`.
