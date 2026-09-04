# Proposition 15.772: the third generic p=1 mod4 layer

Date: 2026-09-04.

Classification: **infinite-family theorem**, for every prime `p=1 mod4`,
`p>=29`: residual (ii) is empty at `t=(p-3)/2`, equivalently `k=5p-3`.
The global residual-(ii), E1, and original-MO gates remain open.

## 1. The new implication and the repaired premise

Write `q=(p-1)/2`, `m=q+1`. The new layer is `t=q-1`, so `|H|=5p-2`.
There are at least `p^2-10p+5>0` isolated vertices. Signed PSL transport
therefore gives the usual all-finite chart with infinity degree `I=0`;
every directional parity boundary has even size. Choose the hard sign `h`
to have phase one. The exact budgets are

```text
a_L=2u+(p+1)k_L,          0<=u<=q,
sum_hard k_L=2q-u,
hT=(p+1)P_L-3p-a_L.                                      (1)
```

The genuinely new local result is
[the complement-triple punctured gap theorem](NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md).
For all odd `p>=29`, a nonnegative integral quadratic of parity
`r=|X intersect C|`, `|C|=3`, and mean `2p E[A]=2p-6+delta`,
`0<=delta<=4`, has only these possibilities:

```text
delta=0: A=(r-2)^2,                              offset 2;
delta=4: A=2-x_i-x_j-x_k+2x_i x_j,               offset 4.
```

The second line has exactly three choices of the pair in `C`. In
particular, gap two is impossible. The proof does not assume
`(A-(r-2)^2)/2>=0` on the omitted `r=0` layer: that assumption is false
in general. Positive quadrature, a prime-free neighboring-slice lower
bound, and the degree-two fixed-weight kernel replace it. This also
supplies the missing gap-two justification in Proposition 15.770;
its already stated range does not change.

## 2. Exhaustive hard residue classification

The exact phase-one floors are `p-1` at `b=2`, `p+1` at `b=p-1`,
`2p-6` at `b=p-3`, and `2p` elsewhere. If `u<q`, a quotient zero
would give `a_L=2u<p-1`, impossible. Thus every quotient is positive,
and at least `u+2` equal one. At `u=q`, the quotient sum is `q<m`,
so at least one quotient is zero. These low rows suffice:

| Residue | Forced low mean | Surviving low cells |
|---|---|---|
| `0` | `p+1` | exact complement literal, offset `5` |
| `1,...,q-4` | `p+1+2u` | none: genuine subsharp lifts |
| `q-3` | `2p-6` | exact complement triple, offset `2` |
| `q-2` | `2p-4` | XNOR plus sharp lift, offsets `3,5` |
| `q-1` | `2p-2` | new pair-plus-literal or literal plus sharp lift, offsets `4,6` |
| `q` | `p-1` | exact XNOR, offset `4` |

The `b=2,p-1` baselines are pointwise parity minima, so their half
differences really are nonnegative integral quadratics. Below mass `p-3`
Proposition 15.688 applies. At mass `p-3`, the sharp Boolean classification
from 15.769--15.770 gives only omitted pairs and all-equal triples.
The XNOR signed target is `4+z_i z_j`, with offset four. The complement
literal target is `4+z_a`, with offset five. The sharp lift increments
have offsets `-1,+1`, giving the displayed pairs of offsets. Support
overlaps do not change this calculation; no disjointness is assumed.

At `u=q-2`, the other possible `b=p-3` cell has gap two and is excluded
by the new punctured theorem. At `u=q-1`, it instead has gap four and
must be retained. The `b=2` cell there would require mass `p-1`.
For `p=1 mod4`, the paired-cube and stabilizer inequalities at that mass
force `H=(p+3)/4>=8`; all maximizing paired cubes then have mean `1/2`,
contradicting their proved maximum bound three. The height-one case
has influence bound less than six and its density is absent from the
fixed four-bit catalog. Thus this local lift does not exist.

No mean-`2p` classification is being assumed: such a quotient-one mean
appears only at `u=q`, where the forced quotient-zero baseline already
fixes the common row identity for every higher row.

## 3. Normalize every row before using offsets

Let a forced low row have quotient `ell` (one, except zero at `u=q`),
parallel count `P`, and coefficient offset `c`. Equation (1) implies

```text
P_L=P+k_L-ell          for every hard direction,
hard edges=m(P-ell)+(2q-u).                                (2)
```

Nonnegativity of opposite edges bounds
`0<=P<=floor((5p-2-(2q-u)+m ell)/m)<=9`. Isolated-chart coefficient
comparison gives `P=c mod q`, while `q>=14` and `c` is one of
`2,3,4,5,6`. Therefore **P=c**. This common-row argument precedes any
choice of normalization and applies to elevated rows as well as low rows.
Different offsets cannot mix at the same low mean; two distinct families
of offset four are allowed to mix.

The opposite row of parallel count `Q` has

```text
a_opp(Q)=(p+1)Q+hT-3p,
opposite edges=|H|-hard edges.                             (3)
```

Equations (1)--(3) give the complete contradiction ledger:

| Hard low family | P | hT | Forbidden Q and mass | Next Q and mass | Number of next rows at least |
|---|---:|---|---|---|---:|
| complement literal | `5` | `p+4` | `2,6` | `3,p+7` | `5` |
| complement triple | `2` | `8-3p` | `6,14` | `7,p+15` | `9` |
| XNOR sharp lift | `3` or `5` | `(p+1)P-5p+4` | `8-P,12` | `9-P,p+13` | `8` |
| new all-low equality | `4` or `6` | `(p+1)P-5p+2` | `8-P,10` | `9-P,p+11` | `7` |
| quotient-zero XNOR | `4` | `5` | `3,8` | `4,p+9` | `6` |

In each line lower `Q` gives negative mass. The forbidden small positive
mass is below both the nonzero phase-zero boundary floor and the sharp
`b=0` lift floor. Thus all `m` opposite counts reach the next `Q`.
Subtracting this baseline from the opposite total leaves `m-d` units,
where `d` is the final column, forcing at least `d` exact next rows.

At any next mass, phase-zero floors leave only nonzero boundaries
`b=2,p-1`. Their pointwise baselines `(x_i-x_j)^2` and `1-x_j` leave positive
lift mass below `p-3`, so they are impossible. Thus a next row must have
`b=0`, namely `A=2C`, with scaled lift mass shown in the table.

## 4. The new local mass p+11 and closure

The old masses `p+7,p+9,p+13,p+15` are excluded by
15.751, 15.752, 15.770, and 15.768 respectively. At the new mass `p+11`,
height at least two implies `H>=(p-9)/4>=5`. Stabilizer averaging bounds
the average maximizing paired-cube mean by `(p+11)/(2(p-1))<3/4`.
Thus some such cube has mean exactly `1/2` and maximum at most three,
contradiction. At height one, the corrected Johnson influence bound is
less than eight. The resulting at-most-seven-variable slice junta
extends to a cube, where degree-two Boolean influence reduces to at most
four active variables. The already pinned four-bit density catalog omits
`(p+11)/(4p)`. The exact inequalities and positive shifted polynomials
are recorded by `p1_p_plus_eleven_local_exclusion` in the local module.

All residues in Section 2 and all boundary sizes are now exhausted.
Therefore every prime `p=1 mod4`, `p>=29` has no residual-(ii) witness
at `k=5p-3`. The p1 frontier moves from `t>=q-1` to `t>=q`.
The p23 and p3 frontiers are unchanged, and no global closure follows.

Replay the proof identities and source/evidence agreement with
`PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 tests/test_complement_triple_gap.py tests/test_prop15770.py tests/test_prop15772.py`.
The mesh checks are bounded complementary verifications of the new
quadrature, contact kernel, equality tables, and row arithmetic, not a
prime, graph, or full-slice census.
