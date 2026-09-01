# Equal-endpoint RG2 and the balanced Paley-skew shield

**Status:** Propositions 6.5--6.6 are proved reductions. They do not prove
the multiplier-two ray, the multiplier-three ray, or convergence.

The general arbitrary-block relative-gauge/max-plus composition framework
already appears in the public
[`Robby955/mo-413935-research`](https://github.com/Robby955/mo-413935-research/blob/main/paper/composition_framework.tex)
project. The new claim here is only the equal-endpoint Hadamard/skew
specialization and its balanced Paley shielding theorem, not the general
framework.

## 1. The original problem gate

For an optimal order-`n` signing `A`, write `M=Phi(A)=m_n`. Proposition
6.3 reduces convergence to Dini-summable amplification at multipliers two
and three. The equal-endpoint Hadamard doubling frame is obtained from a
skew signing `R` by

```text
B_ij = [[A_ij, -R_ij], [R_ij, A_ij]].
```

Both endpoint signings are exactly `A`, and the four-state minimax is

```text
K(A,R) = (1/2) max_(x,y)
         (|Q_A(x)+Q_A(y)| + |x^T R y|).            (1)
```

For every cut, the two hereditary endpoint inequalities from Proposition
6.4 are automatic. If `D` is the sum of the two independently optimized
within-part energies and `X` is the `A`-cross energy, the two full states
have energies `D+X,D-X`; hence `D+|X|<=M`. The same proof applies to the
two minima. Endpoint selection is therefore retired as a live obstruction.

The exact remaining diamond is

```text
|Q_A(x)+Q_A(y)| + |x^T R y|
  <= 2 sqrt(2) M + n^(3/2) Omega(n)                (2)
```

uniformly in the Boolean pair, with vanishing dyadic Dini tail for the
supremum envelope of `Omega`.

## 2. Why separate budgets and the disk are not the target

Every skew signing obeys

```text
max_(x,y) |x^T R y| = max_x ||R x||_1
                    >= n E|S_(n-1)|
                    = (sqrt(2/pi)+o(1)) n^(3/2).   (3)
```

Indeed, each coordinate of `R X` has the law of a sum of `n-1`
independent signs, and expectation commutes with the coordinate sum. An
uncoupled proof using only `|Q_A(x)+Q_A(y)|<=2M` would require

```text
max |x^T R y| <= 2(sqrt(2)-1)M + o(n^(3/2)).
```

The random-method upper bound on `M` makes the right-hand leading constant
`2(sqrt(2)-1)sqrt(log 2)`, strictly below `sqrt(2/pi)`. Thus the proof must
use statewise energy/cross anticorrelation.

The stronger disk surrogate

```text
I_A(x,y)^2 + C_R(x,y)^2 <= M^2+o(n^3),
I_A=(Q_A(x)+Q_A(y))/2, C_R=x^T R y/2,
```

would itself imply `liminf m_n/n^(3/2)>=1/sqrt(2pi)`, well beyond the
proved `1/pi` floor. Its zero-error form is false even for an optimizer at
`n=5`: every order-five signing has `E Q_A(X)^2=10`, and all energies are
even, so `Phi(A)>=4`. Put `-1` on a five-cycle and `+1` on the five
diagonals; it attains `Phi(A)=m_5=4`. The five positive maximizers with
negative coordinates at `{r,r+2}` are the rows of an invertible circulant
matrix `V`. Zero disk
error would give `V R V^T=0`, hence `R=0`, impossible for a skew signing.
This is a route kill for the zero-error disk, not for a Dini-error diamond.

## 3. Balanced near-conference skew signing

The fixed-progression effective prime-number theorem gives constants
`c,N_0>0` such that every `n>=N_0` has a prime

```text
n <= q <= n+n exp(-c sqrt(log n)),  q=3 mod 4.     (4)
```

Let `T_q(a,b)=chi(b-a)` be the Paley tournament matrix. Then

```text
T_q 1=0,                 T_q T_q^T=qI-J,
||T_q||op=sqrt(q).                                  (5)
```

Compress to `n` coordinates, call the matrix `T`, and put `k=q-n`. The
deleted columns give

```text
||T1||_2 <= sqrt(qk),    ||T1||_1 <= sqrt(nqk).    (6)
```

A tournament can be made regular for odd `n`, or near-regular for even
`n`, by repeatedly transferring one outdegree unit from a maximum-degree
vertex to a minimum-degree vertex. Reverse their edge when possible; in the
opposite orientation a directed two-path exists and reversing its two edges
performs the same transfer. The total degree-deviation potential drops at
each step. Thus at most

```text
s <= ||T1||_1 <= sqrt(nqk)                         (7)
```

edge reversals yield a skew signing `R_n=T+F` with

```text
||R_n 1||infinity <= 1,
|u^T R_n v| <= n sqrt(q)+4s                        (8)
```

for every Boolean pair. If `delta=(q-n)/n`, its global relative error is

```text
eps(n)=sqrt(1+delta)-1+4sqrt((1+delta)delta).
```

It is `O(exp(-c' sqrt(log n)))`. The same decreasing majorant bounds its
supremum envelope, and

```text
sum_(j>=0) exp(-c' sqrt(log(2^j n))) -> 0.
```

Hence `eps` has precisely the Dini property required by Proposition 6.3.

## 4. Conjugation by an optimizer and three geometric shields

Choose `z` with `|Q_A(z)|=M` and conjugate

```text
R=diag(z) R_n diag(z).
```

Put

```text
a=(sqrt(2)-1)/pi,
rho=a^2=0.01738...,
tau=(1-sqrt(1-4rho))/2=0.01769....                 (9)
```

The universal Proposition 5.2 floor gives the pairwise margin

```text
Delta_A(x,y)=2sqrt(2)M-|Q_A(x)+Q_A(y)|
            >=2a sqrt(1-1/n)n^(3/2).              (10)
```

### Anchor shield

If `u` differs from `1` in `h` coordinates, row balance and the
`T+F` decomposition give

```text
|u^T R_n v| <= n+2sqrt(qhn)+8s.                   (11)
```

For `h<=rho n`, (11) is at most (10) plus
`O(n^(-1/2)+sqrt(delta))n^(3/2)`. Thus every pair incident to either
Hamming ball of radius `rho n` around `{z,-z}` satisfies the diamond with a
Dini error.

### Diagonal and antidiagonal shield

For `S={i:u_i!=v_i}` and `r=|S|`, skewness gives the exact cut identity

```text
u^T T v = -2 v_S^T T_(S,S^c) v_(S^c).
```

Therefore

```text
|u^T R_n v| <= 2sqrt(q r(n-r))+4s.                (12)
```

This is below (10), with Dini error, whenever
`r(n-r)<=rho n^2`, equivalently `r<=tau n` or
`r>=(1-tau)n`.

### Two-sided anchor product shield

Let `h_x=d_H(x,{z,-z})` and similarly for `h_y`. After independent global
signs, expand both vectors around `1`. Row balance gives

```text
|u^T R_n v|
 <= 2(h_x+h_y)+4sqrt(q h_x h_y)+16s.              (13)
```

Because each `h` is at most `n/2`, (13) is below (10), with Dini error,
when `h_x h_y<=(rho/4)n^2`. In particular it shields the Cartesian product
of the two antipodal anchor balls of radius `(a/2)n`.

The global estimate (8) also closes every pair with
`Delta_A(x,y)>=n^(3/2)`.

## 5. Exact residue and next target

For the explicit orientation above, the only pairs not certified by these
inequalities satisfy all of

```text
h_x,h_y > rho n,
d_H(x,y)(n-d_H(x,y)) > rho n^2,
h_x h_y > (rho/4)n^2,
|Q_A(x)+Q_A(y)| > 2sqrt(2)M-n^(3/2).              (14)
```

This is an infinite-family theorem and an exact residual reduction, not a
finite state census. It does **not** show that (14) is empty.

One natural coupled replacement has also been checked and does not shrink
(14). Orient the positive and negative edge subgraphs of a switched `A`
separately by near-Eulerian orientations, combine them as `T`, and put
`R=A∘T`. This gives both `||T1||infinity<=2` and
`||R1||infinity<=2`. If `x=1-2 1_S`, `y=1-2 1_U`, and
`c=|S intersection U|`, however, the resulting pointwise estimate is only

```text
|x^T R y| <= 4(|S||U|-c^2+|S|+|U|-2c)
          <= 4(|S||U|+|S|+|U|).                  (15)
```

Using the universal floor on `M`, (15) can certify at most
`|S||U|+|S|+|U| <= (a/2)n sqrt(n-1)`. For
`n>(2/a)^2`, every pair in (14) already has product above that threshold.
Thus simultaneous Euler balance is strictly weaker on the live residue than
the Paley spectral product shield; do not reopen it as a residue-closing
mechanism.

If `M/n^(3/2)<=1/(2sqrt(2))-c_0` for fixed `c_0>0`, the last threshold in
(14) is negative for large `n` and hence vacuous. In that regime the global
Paley norm does not remove an additional nonanchor pair; the three geometric
Hamming shields are the substantive gain.

The next multiplier-two target is to choose the Paley principal embedding
and degree-balancing reversals so that the diamond (2) holds on (14) as
well. Independent cross control is impossible by (3). The asymptotic disk
is not disproved, but it is a strictly stronger lower-bound problem; only
its zero-error form is false. Even after multiplier two is closed,
multiplier three (or the `1:2` split in Proposition 6.3) remains necessary
for the original convergence question.

## Replay

```bash
PYTHONPATH=src pytest -q tests/test_direct_rg2_equal_endpoint.py
```

The replay tests the exact four-state formula, automatic hereditary cuts,
skew cut identity, skew-norm floor, and the `n=5` zero-disk counterexample.
