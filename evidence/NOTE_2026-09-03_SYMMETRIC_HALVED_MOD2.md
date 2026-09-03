# The full and punctured halved symmetric map modulo two

Date: 2026-09-03

Status: the full halved map is surjective for every odd prime. Uniform
surjectivity after at most $|Delta|$ arbitrary column deletions is false.
The image after the actual structured Mobius deletion, the exact dual
distance, the direction-weight Boolean fibre, residual (ii), E1, and the
original limit remain open.

## 1. The reduced binary map

Put $p=2h+1$, $d=p+1$, and

    Delta = (F_p^2 minus {0})/{+1,-1},    |Delta|=dh.

After the fixed antipodal choices have been forced, the fixed-cell
remainder is divided by two as in the fixed-edge-elimination note. The
remaining columns are indexed by pairs $([a],[delta])$ in
$Delta$ times $Delta$. Modulo two, all Paley signs disappear.

Write $C$ for the component on paired nonfixed target cells. Apply the
explicit fixed-edge inverse to the divided fixed-cell component and call
the resulting word $Phi$. Row operations do not change image membership,
so the halved map may be studied as

    D = (C,Phi):
        F_2^(Delta times Delta) -> F_2^(d h^2) plus F_2^Delta.       (1)

This is the simultaneous fixed-cell-modulo-four and
paired-cell-modulo-two gate which remains after the first fixed parity
has been solved.

## 2. Reused block-design input

Section 6 of
NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md proves the following
facts. They are imported here, not reproved by a second implementation.

For a nonfixed source orbit with midpoint $a$ and half-difference $delta$,

    Phi(a,[delta]) = 0                                      if a || delta,
    Phi(a,[delta]) = 1_{{[delta+c a]: c in F_p}}            otherwise. (2)

The nonzero words are the $|Delta|$ antipodal pairs of affine lines not
through the origin. If $M$ is their point-block incidence matrix, then

    M M^T = M^T M = I over F_2.                              (3)

Consequently these block words are a basis of $F_2^Delta$.

That note also gives a $C$-kernel lift of every block. Let $B$ be an
antipodal affine block of direction $A$, fix any one of the $h$ midpoint
magnitude classes $[a]$ in $A$, and set

    k_(a,B) = sum_([delta] in B) e_([a],[delta]).            (4)

There are $p$ summands. Every one has $Phi$-word $B$, so their odd sum
still has $Phi(k_(a,B))=1_B$. In the row annihilating $a$ there is no
paired-cell contribution. In every other row, $L(delta)$ runs once
through $F_p$; every nonzero square occurs twice and cancels. Thus

    C k_(a,B)=0,             Phi k_(a,B)=1_B.               (5)

The $h$ choices of $[a]$ give disjoint lifts, although only one lift per
block is needed below.

## 3. Full surjectivity

The inversion-symmetric lattice theorem already proves that

    C: F_2^(Delta times Delta) -> F_2^(d h^2)

is onto. Given a desired pair $(y,z)$ in the codomain of (1), first choose
$x$ with $Cx=y$. By (3)-(5), the restriction of $Phi$ to $ker C$ is onto
$F_2^Delta$. Choose $k$ in $ker C$ with
$Phi(k)=z+Phi(x)$. Then

    D(x+k)=(y,z).

Therefore

    rank D = d h^2 + d h = d h(h+1),                        (6)

and the full unpunctured halved map is surjective for every odd prime.
This closes the unrestricted second parity gate only. It does not put a
binary point in the required punctured, constant-weight fibre.

## 4. Exact puncture criterion

Let $U$ be any set of used nonfixed orbits and let $D_U$ retain precisely
the columns outside $U$. The following are equivalent:

1. $D_U$ is onto.
2. $C_U$ is onto and $Phi(ker C_U)=F_2^Delta$.
3. No nonzero word $w$ in the full row code $Row(D)$ has
   $supp(w)$ contained in $U$.

The third formulation is the exact dual-support criterion: a left
functional becomes a new relation after puncturing exactly when all
columns on which it evaluates nontrivially were deleted.

For a particular reduced target $r_U$, full surjectivity is stronger than
needed. Exact image membership is

    r_U in im(D_U)
    iff lambda(r_U)=0 for every lambda with lambda D_U=0.    (7)

Thus a rank drop is an obstruction only when one of its new dual
functionals pairs nontrivially with the actual target.

## 5. A puncture of size p h destroys surjectivity

Fix a projective functional $L$ and a nonzero square $beta$. Consider the
halved fixed-transverse row indexed by $(L,0,beta)$. Its evaluation word is
supported exactly on

    X_(L,beta)
      = {([a],[delta]): L(a)=0 and L(delta)^2=beta}.         (8)

There are $h$ nonzero antipodal midpoint classes in $ker L$. The two
affine fibres $L(delta)=+sqrt(beta)$ and
$L(delta)=-sqrt(beta)$ combine to one antipodal affine block with $p$
classes. Hence

    |X_(L,beta)| = p h = |Delta|-h <= |Delta|.              (9)

Delete exactly $X_(L,beta)$. The row in (8), which was nonzero for the full
map, is now identically zero. The punctured rank drops, so $D_U$ is not
onto. Define the puncturing distance here to be the minimum support size of
a nonzero word in $Row(D)$; this is the row-code distance which controls
column deletions. The construction proves only

    puncturing distance d_row(D) <= p h.                    (10)

It disproves robustness under every arbitrary deletion of size at most
$|Delta|$. It does not prove equality in (10), and it does not classify
all dual words of weight at most $|Delta|$.

## 6. What the Mobius midpoint bound does and does not say

The fixed-edge-elimination theorem independently proves that one localized
Mobius half has at most two orbit midpoints in any fixed projective
direction. Therefore a union of $m$ such halves meets the rectangle
$X_(L,beta)$ in at most $2m$ columns. For the conditional value
$m=(p+1)/2$, this is at most $p+1$, much smaller than $ph$ when $p>=7$.
So that structured union cannot contain the explicit counter-puncture
(8).

This intersection bound does not establish punctured surjectivity. A
different low-weight dual word could still have support contained in the
actual used set. The next exact linear task is therefore either:

* classify the nonzero words of $Row(D)$ of weight at most $|Delta|$; or
* directly prove that no such support is contained in the actual Mobius
  set $U$.

No equality assertion for this distance is made here.

## 7. Direction weights remain an integer obstruction

Before the fixed-row change of basis, the $P_L$ coordinate of a halved
column is one exactly when its difference is parallel to $ker L$.
Consequently the mod-two equation includes

    sum_(O parallel L) b_O = n_L  (mod 2).                 (11)

The physical graph requires the stronger integer equality

    sum_(O parallel L) b_O = n_L                           (12)

for every one of the $p+1$ directions. Full or punctured parity
surjectivity cannot replace (12). The whole-slab exchanges from
NOTE_2026-09-03_SYMMETRIC_UNUSED_SLICE_EXCHANGE.md preserve these weights;
they cannot repair a solution in the wrong direction slice.

After the punctured parity test, the live obstruction is therefore the
intersection of the affine binary solution set with all exact
direction-weight slices, followed by the full integral target equation.
Nothing here proves that intersection nonempty.

## 8. Fail-when-wrong replay

The implementation builds only the small raw halved matrices for
$p=3,5,7$. It obtains respectively

| $p$ | $|Delta|$ | raw rows | rank $D$ | $ph$ | rank loss |
|---:|---:|---:|---:|---:|---:|
| 3 | 4 | 12 | 8 | 3 | 1 |
| 5 | 12 | 42 | 36 | 10 | 1 |
| 7 | 24 | 104 | 96 | 21 | 1 |

In every case the paired component retains rank $dh^2$ after the explicit
counter-puncture. These finite computations check the formulas only.

Reproduction:

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_symmetric_halved_mod2.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_halved_mod2.py
