# `p=31`: exact audit of the eleven public complete 22-arc classes

**Result status: exhaustive finite certificate for the eleven publicly
sourced classes, not a theorem and not endpoint closure.**

Proposition 15.727 shows that an endpoint configuration at `p=31,R=10`
would require a 22-arc `A` with at least ten outside points of secant index
one:

```text
c_1(A) >= 10.
```

Coolsaet's corrected exhaustive classification reports twelve projective
classes of complete 22-arcs in `PG(2,31)`. Keri's older public coordinate
file contains only the eleven classes known before that correction. The
present certificate audits every representative in that public file. It
does not have coordinates for, and makes no assertion about, the corrected
twelfth class.

## Public input and coordinate convention

The checked-in input is a verbatim transcription of Keri's public A-list:

- local: `evidence/data/31q3x19.txt`;
- public coordinates:
  <https://old.sztaki.hu/~keri/n-arcs/31q3x19.txt>;
- A-list instructions:
  <https://old.sztaki.hu/~keri/n-arcs/instructions.pdf>;
- original table reporting eleven classes:
  <https://old.sztaki.hu/~keri/n-arcs/n-arcs_in_PG(2,31).pdf>.

For a prime field the A-list entries are ordinary residues. Each displayed
`2 x 18` matrix with rows `(a_i)` and `(b_i)` expands to

```text
e_1, e_2, e_3, (1,1,1), (1,a_i,b_i) for 1 <= i <= 18.
```

The corrected count of twelve is from K. Coolsaet, *The Complete Arcs of
PG(2,31)*, J. Combin. Des. **23** (2015), 522--533,
[doi:10.1002/jcd.21410](https://doi.org/10.1002/jcd.21410). The
[UGent record](https://biblio.ugent.be/publication/7076091) exposes no public
coordinate supplement; its sole deposited full-text file is restricted.

## Exact audit

For each of the eleven representatives, the script verifies:

1. there are 22 distinct projective points;
2. its 231 point pairs determine 231 distinct lines, so it is an arc;
3. all 971 outside projective points have positive secant index, so the arc
   is complete;
4. the outside-index histogram has total 971 and first moment
   `C(22,2)*(31-1)=6930`;
5. every index-one point and its unique secant are enumerated exactly.

In Keri file order the resulting `c_1` sequence is

```text
0, 2, 0, 0, 0, 0, 0, 2, 2, 0, 0.
```

Thus all eleven public classes have `c_1<=2<10` and fail the necessary
condition from 15.727.

There is additional rigidity in the three classes with `c_1=2`. In each
case the two index-one outside points lie on the same unique secant:

| public class | secant line | arc pair | index-one outside points |
|---:|---|---|---|
| 2 | `(1,28,8)` | `(1,13,28)`, `(1,14,9)` | `(1,24,5)`, `(1,28,22)` |
| 8 | `(0,1,16)` | `(1,0,0)`, `(1,21,20)` | `(1,19,24)`, `(1,20,22)` |
| 9 | `(1,2,20)` | `(1,6,4)`, `(1,9,13)` | `(1,22,21)`, `(1,25,30)` |

Consequently the largest matching of pairwise-disjoint unique secants has
size one in classes 2, 8, and 9, and zero in the other eight classes. This
is stronger than needed for these eleven representatives, but still says
nothing about the missing class.

## Scope boundary

This certificate is exhaustive over the eleven records in the public Keri
file. It is **not** exhaustive over Coolsaet's corrected twelve-class
classification. In particular:

```text
twelfth representative available = false
twelfth class audited             = false
all twelve classes excluded       = false
p=31 endpoint closed              = false
```

Closing this finite side route requires an authoritative representative for
the twelfth class, or a theorem bounding `c_1` for all complete 22-arcs
without enumerating representatives. The all-prime endpoint reduction
remains the primary route.

## Reproduction

```bash
python scripts/p31_complete_22arc_public_audit.py
python -m pytest -q tests/test_p31_complete_22arc_public_audit.py
```

The generated certificate is
`evidence/p31_complete_22arc_public_11_audit.json`.
