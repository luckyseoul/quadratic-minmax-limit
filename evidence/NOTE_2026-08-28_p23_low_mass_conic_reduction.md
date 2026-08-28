# The `p=23` second endpoint is reduced from seven residues to 203 exact profiles

Date: 2026-08-28. This is Proposition 15.684. It is a strict reduction of the
`p=23,s=20` second all-finite endpoint, not a closure. Every positive
phase-zero residue is impossible. Of the 1,247 exact residue-zero profiles,
1,044 are impossible and exactly 203 remain.

The two main ingredients are a sharp low-value analysis of nonnegative
integral quadratics on `J(23,12)` and the exhaustive complete-arc
classification of `PG(2,23)`.

## 1. Exact residue ledger

Here

```text
p=23, m=12, P=p+1=24, s=20, pair-deficit budget=s(s-1)=380.
```

The exact type minima leave only phase-one residue `u_1=11`. Pairing it with
phase zero leaves

```text
u_0 = 0,2,3,4,5,6,8.
```

A positive residue forces a quotient-zero direction. Its scaled mean is
`c=2u_0`, which is smaller than the least positive-`b` phase-zero floor 24,
so this direction has `b=0` and factors as `A=2B`, where `B` is a nonzero,
nonnegative, integer-valued quadratic on `J(23,12)`.

Proposition 15.681 gives `c>=12`, immediately excluding `u_0=2,3,4,5`.
The cases left by that bound are

```text
u_0=6: c=12,
u_0=8: c=16.
```

## 2. Paired-cube exclusion through height three

Write `H=max B`. The `p=23` stabilizer identity gives

```text
c >= 4H.                                             (1)
```

Choose a maximum point `X`. Fix one element of `X` and pair the other eleven
elements with the eleven elements outside `X`. Proposition 15.681's paired
cube has dimension 11 and mean

```text
T B(X) = (H+c/4)/24.                                 (2)
```

The restricted function `f` is integer-valued, nonnegative, degree at most
two, takes values in `{0,...,H}`, and has `f(0)=H`.

Three elementary cube facts suffice:

1. A nonzero degree-two function has support density at least `1/4`.
2. If a nonnegative integral degree-two function has mean exactly `1/4`,
   every nonzero value must be one. It therefore cannot take value two.
3. If `f` takes values in `{0,1,2,3}` and `f(0)=3`, then

   ```text
   f = (f mod 2) + 2*(binom(f,2) mod 2).
   ```

   Integer-valued cube functions have integer Möbius coefficients. The two
   bits above are nonzero binary functions of degrees at most two and four.
   The elementary Reed--Muller distance bound therefore gives

   ```text
   E[f] >= 1/4 + 2*(1/16) = 3/8.                    (3)
   ```

For `c=12`, equation (1) gives `H<=3`. The means in (2) for `H=1,2,3` are
`1/6,5/24,1/4`, contradicting respectively the `1/4,1/4,3/8` floors.

For `c=16`, equation (1) gives `H<=4`. The means for `H=1,2,3` are
`5/24,1/4,7/24`; the first and third are below their floors, while the middle
is the forbidden value-two equality case. Only `H=4` needs another idea.

## 3. The height-four shell kernel

At `c=16,H=4`, the exact stabilizer identity at a maximum is

```text
E[B] = (22/23) q(6) + B(X)/23.
```

Both sides equal `4/23`, hence `q(6)=0`. Nonnegativity makes this pointwise:

```text
B(Y)=0 whenever |X intersect Y|=6.                  (4)
```

The degree-at-most-two function space on `J(23,12)` has dimension

```text
dim V_2(J(23,12)) = C(23,2) = 253.
```

The shell in (4) is `J(12,6) x J(11,6)`. Its filtered degree-two image has
dimension

```text
1 + 11 + 10 + 54 + 44 + 11*10 = 230.
```

Restriction is onto this filtered product space because coordinate monomials
of total degree at most two span it. Thus the kernel has dimension 23.

Multiplication by `t-6`, where `t=|X intersect Y|`, maps the 23-dimensional
space `V_1(J(23,12))` into this kernel. The map is injective: if an affine
`L` vanishes off `t=6`, then its vanishing on the `t=5` shell makes its
coefficients constant inside each of the two blocks, and vanishing also at
`t=7` gives two roots for the resulting affine function of `t`. Hence `L=0`.
Dimension equality and (4) give

```text
B(Y) = (|X intersect Y|-6) L(Y)                     (5)
```

for an affine function `L`.

At `X`, equation (5) gives `L(X)=2/3`. For one replacement, write

```text
b_ij = B(X-i+j),        L(X-i+j)=b_ij/5.
```

For distinct `i_1,i_2` and `j_1,j_2`, the affine parallelogram identity gives

```text
B(X-i_1-i_2+j_1+j_2)
  = 4*(3*(b_i1j1+b_i2j2)-10)/15.                   (6)
```

The left side is an integer. The right side cannot be: divisibility by 15
would require `3n-10` to vanish modulo three. This excludes `c=16,H=4` and
therefore every positive residue.

## 4. Exact residue-zero census

Completion-bounded exact enumeration gives 426 phase-zero rows and 11
phase-one rows within the pair budget. Pair compatibility and slack modulo
four give 1,247 phase-labelled profiles and 485 global shapes:

```text
slack   profiles   shapes
    0        363      124
    4        264       95
    8        189       72
   12        136       54
   16         94       39
   20         68       30
   24         49       23
   28         35       17
   32         21       11
   36         13        8
   40          7        5
   44          4        3
   48          1        1
   52          1        1
   56          1        1
   60          1        1
```

The canonical phase-labelled profile ledger has SHA-256

```text
19ea72e792303d42863d327114eea6edde0abb3039a578b991387ead83fa5cc0
```

For a line containing `n` boundary points, its contribution to pair slack is

```text
delta(n) = 2*(C(n,2)-floor(n/2)).                   (7)
```

Thus slack zero is an arc, slack four means one 3-secant, and so on.

## 5. All 363 arc profiles

A 20-arc in `PG(2,23)` has `tau=23+2-20=5` tangents through every point.
Ball--Lavrauw Theorem 11 supplies a degree-10 dual tangent envelope. Call a
direction high if it has at most four secants. It then has at least 12
tangents, so its dual direction line divides the envelope twice.

Every exact arc profile has between three and five high directions. Let `d`
be their number and `E` the number of their secant edges. After removing the
`d` direction squares, the residual degree is

```text
r = 10-2d.
```

If an arc point is incident with `e>=1` high edges, then `d-e` high directions
are tangent there. Its other `5-d+e` tangents remain as distinct double roots
on its point-pencil, with total multiplicity

```text
2*(5-d+e) = r+2e > r.
```

That point-pencil is therefore a residual line component. After removing it,
every other arc point has at least `r` surviving tangent multiplicity on a
curve of degree `r-1`; its point-pencil is forced too. This produces more line
components than the residual degree. It excludes all 320 profiles with
`E>0`.

The other 43 profiles have `d=3,E=0`: the three high directions are all
undetermined. Adjoin two of their infinity points. The resulting 22-arc is
conic-contained by the classification below. Doing this with two pairs gives
two conics sharing the original 20-arc, hence the same conic. It would contain
all three collinear infinity points, impossible. All 363 arc profiles are
excluded.

## 6. Complete arcs and the conic-core repair lemma

Coolsaet and Sticker's exhaustive classification lists complete arcs in
`PG(2,23)` only at sizes

```text
10,12,13,14,15,16,17,24.
```

There are no complete arcs of sizes 18 through 23, and the unique 24-arc is a
conic. Since the plane is finite, every arc of size at least 18 extends to a
complete arc and is therefore conic-contained.

Equation (7) also gives a useful repair bound. A line of occupancy `n>=3`
can be repaired by deleting `n-2` points, and

```text
n-2 <= delta(n)/4.
```

Process every bad line and take the union of the chosen deletions. Deletion
cannot create a new bad line. Hence pair slack `4r` permits deleting at most
`r` points to obtain an arc.

Once the repaired arc is conic-contained, slack below 24 is impossible. To
see this, let `C` be its conic and suppose `h`, with `1<=h<=4`, of the original
20 points lie off `C`. Every off-conic point lies on at least 11 full secants
of `C`. The `20-h` retained conic points omit `4+h` points of `C`, destroying
at most that many secant pairs. Thus each off-conic point sees at least
`7-h` secants whose two conic points remain.

If one such line contains `r` off-conic points, its occupancy is at least
`2+r`, and (7) contributes at least `4r`. Summing off-conic/secant incidences
gives

```text
pair slack >= 4h(7-h) >= 24.                        (8)
```

If `h=0`, the original set lies on the conic and is itself an arc. Therefore,
after reaching a conic core, every positive slack below 24 is contradictory.

Apply this as follows:

```text
slack  4: delete <=1, obtain an arc of size >=19; all 264 profiles excluded.
slack  8: delete <=2, obtain an arc of size >=18; all 189 profiles excluded.
slack 12: delete <=3, obtain size >=17; one undetermined infinity point
          reaches size 18; 135 of 136 profiles excluded.
slack 16: delete <=4, obtain size >=16; two undetermined infinity points
          reach size 18; 93 of 94 profiles excluded.
```

Together with the arc argument, this excludes

```text
363 + 264 + 189 + 135 + 93 = 1044
```

of the 1,247 residue-zero profiles.

## 7. Exact remainder

The two low-slack exceptions are

```text
slack 12:
  phase 0: 8*b0 + 4*b18, deficit 168
  phase 1: 11*b2 + b18,   deficit 200
  floor-secants: {t1:5, t9:11, t10:8}; no undetermined direction.

slack 16:
  phase 0: 7*b0 + b2 + 3*b18 + b20, deficit 164
  phase 1: 11*b2 + b18,                 deficit 200
  floor-secants: {t0:1, t1:4, t9:12, t10:7}; one undetermined direction.
```

The full remaining histogram is

```text
{12:1, 16:1, 20:68, 24:49, 28:35, 32:21, 36:13,
 40:7, 44:4, 48:1, 52:1, 56:1, 60:1}.
```

It sums to exactly 203. The `p=23` endpoint remains open on precisely these
profiles.

## 8. Literature, OEIS, and reproduction

The external inputs are:

- K. Coolsaet and H. Sticker, *A full classification of the complete k-arcs
  of PG(2,23) and PG(2,25)*, J. Combin. Des. **17** (2009), 459--477,
  doi:10.1002/jcd.20211.
- H. Sticker, *Classification of Arcs in Small Desarguesian Projective
  Planes*, Ghent PhD thesis, 2012, Section 5.1, which independently tabulates
  the same complete-arc spectrum and counts.
- S. Ball and M. Lavrauw, *Planar arcs*, J. Combin. Theory Ser. A **160**
  (2018), 261--287, Theorem 11, doi:10.1016/j.jcta.2018.06.015.

Exact OEIS searches for the profile-count block
`363,264,189,136,94,68`, the reduction block `1247,1044,203`, and the large
classification-count block `112449,4341514,1828196` returned no sequence.
Individual values occur in unrelated entries, so no OEIS claim is made.

Reproduce with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15684.py
PYTHONPATH=src pytest -q tests/test_prop15684.py
```

The generated machine-readable record is
`evidence/e1_gmin_m4_prop15684.json`.
