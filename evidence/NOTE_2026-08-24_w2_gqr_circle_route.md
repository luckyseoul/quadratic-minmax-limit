# W2 via the GQR circle code and a Frobenius-pair completion

Date: 2026-08-24

Status: one new all-prime linear theorem; nonlinear sign-completion lemma
open.  W2 remains open.

## 1. Why the inversive-plane code is the right ambient object

Let `C` be the Paley conference matrix on

\[
P^1(\mathbb F_{p^2}),\qquad n=p^2+1,
\]

and let `H0` be the binary direction space of the Max-minus ensemble.
Props. 15.603--15.609 identify `H0` with the binary code spanned by the
nonsquare orbit of Miquelian circles.  Van Lint's generalized
quadratic-residue (GQR) code gives the same identification: the two PSL
orbits of circles are the minimum-weight supports of the two extended GQR
codes, and each orbit spans its code.

Fix the edge

\[
e=\{\infty,0\},\qquad \ell_e(x)=x_\infty+x_0.
\]

Call a nonsquare circle *eligible* when it meets `e` in zero or two points.
Its incidence word lies in `H0 cap ker(ell_e)`, the target direction space
of the xor slice `U`.

## 2. New theorem: eligible circles span the entire edge slice

**Theorem.**  For every odd prime `p`, the binary incidence words of the
eligible nonsquare circles span

\[
H_0\cap\ker(\ell_e).
\]

The dimension is `(p^2-1)/2`.

### Proof

Let `K` be the span of the eligible nonsquare circles.  A tangent pencil of
`p` same-type circles sums over `F2` to the all-ones word: its circles share
their carrier and partition every other point.  Take a nonsquare circle
through both endpoints of `e` and a third point `x` on it.  Every circle in
the tangent pencil at `x` meets `e` evenly.  Hence `1 in K`.

Now form a bipartite graph.  Its left vertices are nonsquare circles through
`infinity` but not `0`; its right vertices are nonsquare circles through `0`
but not `infinity`; join two vertices when the circles are tangent outside
`e`.  If `A` and `B` are joined, the tangent-pencil relation at their common
point has `A` and `B` as its only edge-incompatible circles.  The other
`p-2` circles and the all-ones word are in `K`, so

\[
A+B\in K.
\]

It remains to see that this bipartite tangency graph is connected.  Put

\[
M=(\mathbb F_{p^2}^{\times})^2,
\qquad
A=\{\infty\}\cup(1+\mathbb F_p\sigma),
\]

where `sigma` is nonsquare.  Square dilations label the two vertex classes
regularly by

\[
A_m=mA,\qquad B_t=I(tA),\qquad m,t\in M.
\]

After dilating by `m^{-1}`, tangency of `A_m` and `B_t` is tangency of `A`
and `I(mt A)`.  The involution `x -> 1/(mt x)` must fix their unique
intersection point.  Therefore the connection set is

\[
T=\{(1+r\sigma)^{-2}:r\in\mathbb F_p\}\subset M.
\]

It contains `1`.  If `T` lay in a proper subgroup of `M`, a nontrivial
multiplicative character `theta` of `M` would be trivial on `T`.  Therefore

\[
\eta(z)=\theta(z^{-2})
\]

is a nontrivial even character of `F_(p^2)^times` which is identically one
on the affine line `1+F_p sigma`.

Here is the precise reduction to Katz's estimate.  After scaling the line
and translating its parameter, its character sum is a constant of modulus
one times

\[
\sum_{r\in\mathbb F_p}\eta(r+\delta),
\qquad \delta^p=-\delta\ne0.
\]

Because `eta(-1)=1`, there are characters `epsilon` of `F_p^times` and
`omega` of the norm-one torus `N` such that

\[
\eta(z)=\epsilon(Nz)\,\omega(z/z^p).
\]

The substitution `x=(r+delta)/(r-delta)` bijects `F_p` with
`N minus {1}`.  Since

\[
N(r+\delta)=\frac{4\delta^2}{x+x^{-1}-2},
\]

the line sum is, up to a constant, Katz's Soto--Andrade sum

\[
\sum_{x\in N}\epsilon^{-1}(x+x^{-1}-2)\omega(x).
\]

The omitted point `x=1` contributes zero under the standard extension of a
multiplicative character to zero.  Theorem 1 of Katz bounds this by
`2 sqrt(p)`, except when `epsilon^{-1}` and `omega` are both trivial or both
quadratic.  Those are exactly the two pairs in the kernel of
`(epsilon,omega) -> eta`, so either exception would make `eta` trivial.
They are unavailable here.  A sum identically equal to `p` is therefore
impossible for `p>=5`.  For `p=3`, the three distinct elements of `T` cannot
fit in a proper subgroup of the four-element group `M`.  Thus `T` generates
`M` and the tangency graph is connected.

Consequently all edge-incompatible nonsquare circles have the same class
modulo `K`.  Any word of `H0 cap ker(ell_e)` is a sum containing an even
number of edge-incompatible generators, whose common quotient class
cancels.  Hence it lies in `K`.  The reverse inclusion was immediate.

Exact rank checks give `4,12,24,60` at `p=3,5,7,11`.

## 3. Sparse real eigenvectors and circle flips

Every nonsquare circle `S` has a unique signing `v_S` up to sign such that

\[
\operatorname{supp}(v_S)=S,\qquad Cv_S=-p v_S.
\]

If `y` is a Max-minus sign vector, flipping `y` on `S` remains Max-minus
exactly when

\[
y|_S=\pm v_S|_S.
\]

Indeed, with `y' = y-2v_S` after choosing the matching orientation,
`Cy'=-py'`.  Conversely the difference of two Max-minus points at the
minimum Hamming distance `p+1` is such a signed circle vector.  Van Lint's
minimum-word classification identifies its binary support with a nonsquare
circle.

If `S` meets `e` in two points, every completion automatically has the `U`
sign on `e`, because the switched principal block on `S` is `-(J-I)`.
If `S` is disjoint from `e`, one needs a completion with

\[
C_{\infty0}y_\infty y_0=-1.
\]

Thus the theorem in Section 2 would close the W2 slice if every eligible
circle were realised by a flip inside `U` (or if a realised subset generated
the eligible-circle code).

## 4. The nonlinear problem collapses to one exceptional pair orbit

Fix the standard nonsquare circle

\[
S=\{\infty\}\cup\mathbb F_p\sigma.
\]

Its stabilizer is conjugate to `PGL(2,p)`.  On unordered pairs outside `S`
there is one Frobenius-conjugate orbit of size `p(p-1)/2` and `p-2` generic
orbits of size `p(p^2-1)/2`.

The named halfspace Max-minus vector of Prop. 15.613 is a completion of this
standard circle for every tested prime through `29`.  Exact orbit censuses at
`p=3,5,7,11,13` show that its stabilizer orbit meets every generic
outside-pair orbit.  On the Frobenius-conjugate orbit it has the correct `U`
sign exactly in the tested `p == 1 (mod 4)` cases; for `p == 3 (mod 4)` a
second completion is required.  This explains why the named-circle family by
itself misses the odd W2 factors at `p=7,11`.  The generic-orbit coverage has
an explicit cross-ratio/character formulation but is not yet recorded as an
all-prime proof.

The remaining geometric target is therefore:

> **Frobenius-pair completion lemma.**  For every odd prime `p>=5`, the
> standard nonsquare circle has a Max-minus sign completion with the `U`
> sign on one (hence every, by its stabilizer) Frobenius-conjugate outside
> pair.

Together with Section 2 and an all-prime proof of the generic named-orbit
coverage, this would be a scalable W2 route: it has no factor degree,
prime-order, or boundary-prefix parameter.

## 5. Exact computations and holdout

The scripts are:

- `scripts/w2_circle_completion_probe.py`: exact sparse-eigenvector check,
  chunked ensemble scan, and pair-orbit completion counts;
- `scripts/w2_circle_completion_milp.py`: exact binary linear feasibility
  formulation for the Frobenius-pair completion, with HiGHS and CP-SAT
  backends and exact integer validation.

For `p=5` and `p=7`, every eligible circle occurs as a difference of two
points in `U`; their ranks are respectively `12` and `24`.

The independent `p=11` holdout used all `37,457,112` normalized rows in
`/home/nick/e1work/maxplus_p11/maxplus_p11_eps1.npy` on NUKA.  For one
standard circle it found `369,302` normalized completions.  All 66 internal
pairs always have the `U` sign.  Every one of the 5,995 outside pairs has a
`U` completion; the minimum count is `179,760`.  Independently, the 551
eligible nonsquare circles have binary rank `60`, exactly the target.

At `p=19`, the named halfspace vector has no finite flippable nonsquare
circle and does not supply the exceptional orbit.  The exact feasibility
problem has 181 independent equations and 340 free binary variables after
fixing the standard circle and one Frobenius pair.  HiGHS and 16-worker
CP-SAT each timed out at 300 seconds on both pair orientations.  This is
`UNKNOWN`, not infeasible.  The shortened binary GQR code has dimension 161
and can toggle the pair, so there is no mod-two obstruction; the unresolved
step is the lift to an exact real `-p` eigenvector.

## 6. Literature and sequence checks

The relevant literature is:

- J. H. van Lint, *Generalized Quadratic-Residue Codes* (1979), especially
  Theorem 4.4 classifying the Miquelian circles as the minimum-word supports;
- N. Katz, *Estimates for Soto-Andrade Sums*, J. reine angew. Math. 438
  (1993), Theorem 1 for the `2 sqrt(p)` estimate;
- A. C. Kable, *Legendre sums, Soto-Andrade sums, and Kloosterman sums*
  (2002), for the finite `PGL(2,p)` harmonic-analysis dictionary.

OEIS/web checks on the new completion counts (`52`, `628`, `369302`) and
the full-ensemble counts found only unrelated numerical coincidences.  No
catalogued code/design sequence was identified.

## 7. Status

- Eligible-circle span theorem: **proved for every odd prime**.
- Every eligible circle realised inside `U`: certified `p=5,7,11`, not a
  p-law.
- Frobenius-pair completion lemma: certified `p=5,7,11`; `p=19` unknown.
- W2: **OPEN**.
