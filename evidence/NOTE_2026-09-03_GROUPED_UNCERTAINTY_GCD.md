# A homogeneous-form gcd for the grouped uncertainty obstruction

Date: 2026-09-03

Status: exact reduction and odd-support theorem.  The even-support inequality
is **OPEN**.  In particular this note does not close residual (ii).

## 1. Setup

Let

    Delta=(F_p^2\{0})/{+1,-1}

and let `S` be a nonempty subset of `Delta`, of size `s`.  Choose one
representative `v_i` of each antipodal class.  A projective direction `A` is
the kernel of a nonzero linear functional `u`, defined up to scalar.  Put

    ell_i(u)=u(v_i).

The paired nonorigin affine blocks in direction `A` are exactly the fibres of
the nonzero squared projection values

    ell_i(u)^2 in F_p^*.

Thus `A` is silent for the block incidence transform precisely when every
nonzero squared-projection fibre has even cardinality.  Classes with zero
projection lie on the origin line and impose no parity condition on the
original antipodally doubled support.

## 2. Detecting forms

For each support class define the split homogeneous binary form

    F_i(u)=ell_i(u) product_(j!=i)(ell_i(u)^2-ell_j(u)^2).       (1)

It has degree `2s-1`.  For a split form `F`, write `odd(F)` for the product of
the projective linear factors which occur in `F` with odd multiplicity.

Fix a direction `A`.  Suppose first that `i` is one of `r` radial classes,
so `ell_i(A)=0`.  The factor `ell_i` has order one and the difference of
squares with each other radial class has order two.  Hence

    ord_A(F_i)=1+2(r-1)=2r-1.                                 (2)

If instead `i` belongs to a nonzero squared-projection fibre of size `m`,
exactly `m-1` difference-of-square factors vanish simply, and

    ord_A(F_i)=m-1.                                            (3)

The order in (2) is always odd.  The order in (3) is odd exactly when `m` is
even.  It follows that

    A is silent  iff  A divides odd(F_i) for every i.           (4)

All factors in (1) split over `F_p`: they are the radial factor `ell_i` and
the signed-chord factors `ell_i-ell_j` and `ell_i+ell_j`.  Consequently, if

    D=gcd_i odd(F_i),                                          (5)

then there is no extension-field or multiplicity ambiguity:

    number of silent directions = deg D.                       (6)

The point at infinity is included by keeping (1) homogeneous.

## 3. The odd-support branch

Multiplying (1) over `i` gives, up to a nonzero scalar,

    product_i F_i
      =(product_i ell_i)
       product_(i<j)(ell_i^2-ell_j^2)^2.                       (7)

If `s` is odd, every factor of `D` occurs oddly in the left side of (7).
The squared factor on the right contributes nothing to the odd part, so

    D divides odd(product_i ell_i).

The latter has degree at most `s`.  Equations (5)-(6) therefore prove

    z <= s                                                     (8)

for odd `s`.

This is the homogeneous-form version of the shorter radial argument.  In a
silent direction the number of radial support classes is congruent to `s`
modulo two.  If `s` is odd, every silent direction contains at least one
support class on its origin line, and distinct radial directions account for
at most `s` directions.

## 4. Why the same product does not settle even support

When `s` is even, a common factor `D` occurs an even number of times in
`product_i F_i`; it disappears into the square in (7).  Thus (7) gives no
degree bound on `D`.  This cancellation is exact, not a missing sign.

There is nevertheless a sharp necessary condition for a counterexample.  If
`z>s`, then each degree-`2s-1` form is divisible by the common degree-`z`
factor `D`.  The `s` quotients have degree at most

    2s-1-z <= s-2,

so they lie in a space of homogeneous binary forms of dimension at most
`s-1`; hence they must be linearly dependent.  For the first currently live
case `s=8,z=9`, a counterexample forces a common degree-nine squarefree factor
and eight dependent homogeneous sextic quotients.  This is a reduction, not
an exclusion.

## 5. Why Lev's rational theorem cannot simply be imported

In the binary group algebra `F_2[F_p^2]`, silence in direction `H` means that
the support indicator lies in the relative augmentation ideal for `H`.
Because `2` does not divide `p`, the algebra is semisimple, and simultaneous
silence is membership in the intersection (equivalently product) of those
directional ideals.  The desired inequality is therefore a minimum-support
statement for a binary abelian code whose defining zero set is a union of
dual frequency lines.

Lev's theorem, [Point distribution and perfect directions in
`F_p^2`](https://arxiv.org/abs/1903.01518), instead starts with a nonzero
rational or real weight, supported on the same point set, whose line sums
vanish exactly.  Binary even line sums provide only a kernel vector modulo
two.  A restricted incidence matrix can have 2-primary torsion, so this
kernel need not lift to a characteristic-zero kernel on the same support.
Orthogonally projecting an integral lift onto the rational zero-line-sum
space also does not fix the problem: it can introduce nonzero even
coefficients outside the original support, destroying the support bound
needed by Lev's uncertainty inequality.

Accordingly, a valid Lev bridge would itself have to prove that the relevant
restricted incidence kernel lifts without enlarging support.  No such lift is
proved here, and characteristic-zero uncertainty is not used as evidence for
the even-support claim.

## 6. Replay

The implementation records (2)-(8), checks several local partitions, and
freezes the `s=8,z=9` quotient-dimension consequence:

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      python -m pytest -q tests/test_grouped_uncertainty_gcd.py

Files:

    src/e1_gmin_m4_grouped_uncertainty_gcd.py
    tests/test_grouped_uncertainty_gcd.py
