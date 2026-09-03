# Mobius exclusion of the explicit weight-|Delta| dual words

Date: 2026-09-03

Status: vertical fibres and scalar graphs are excluded from the union, and
hence from the actual post-cancellation support, of `(p+1)/2` nontrivial
localized Mobius halves. The coordinate proof is symbolic for every odd prime
whenever the displayed Mobius parameterization is used; its Paley application
has `p=3 mod 4`. It does not classify other low-weight row-code words or prove
punctured surjectivity.

## 1. One-half coordinates

Fix independent functionals `L,M` and a nonzero center `j`. Put `z=t+1`, so
`z` runs through `F_p^*`. From the parameterized Mobius edge, oriented before
passing to antipodal classes, the midpoint and half-difference have `(L,M)`
coordinates

    a_z     = (j/2) (z,   z-z^(-1)),
    delta_z = (j/2) (z-2, z-2+z^(-1)).                       (1)

All conclusions below are therefore independent of any coherence among the
chosen directions, auxiliaries, or centers.

## 2. Vertical fibres

The map `z -> [delta_z]` is injective. Indeed, suppose
`delta_z' = epsilon delta_z` for `epsilon` in `{+1,-1}`. The `L` coordinate
gives `z'=z` when `epsilon=1`. When `epsilon=-1`, it gives `z'=4-z`, while the
`M-L` coordinate gives `z'=-z`. These would force `4=0`, impossible for odd
`p`.

Thus one half meets any vertical fibre

    V_[d] = Delta times {[d]}                                 (2)

in at most one column. There is an exact intrinsic hit test. For a chosen
representative `d`, a hit exists precisely when, for one sign `epsilon`,

    4 (M(d)-L(d)) (L(d)+epsilon*j) = j^2.                    (3)

The sign is unique, and the parameter is

    z = 2 + 2 epsilon L(d)/j.                                (4)

For a union of `m` halves, equality in the bound `|V_[d] intersect union|<=m`
requires every half to pass (3) and the resulting midpoint classes to be
pairwise distinct.

## 3. Scalar graphs

For a nonzero scalar class `[c]` modulo sign, put

    G_[c] = {([a],[delta]) : [a]=[c delta]}.                  (5)

Taking the determinant of the two vectors in (1) gives

    det_(L,M)(a_z,delta_z) = j^2 (z-1)/(2z).                 (6)

Consequently `a_z` and `delta_z` are proportional only at `z=1`. At that
parameter `a_1=-delta_1`. Each half therefore meets the identity scalar graph
`G_[1]` in exactly one column and every nonidentity scalar graph in zero
columns. A union of `m` halves meets `G_[1]` in at most `m`; equality means
that their `z=1` diagonal columns are pairwise distinct.

## 4. Containment consequence and limits

Let

    h=(p-1)/2,  m=(p+1)/2=h+1,
    |Delta|=(p+1)h=m(p-1).                                  (7)

Every vertical fibre and scalar graph has `|Delta|` columns. Since
`m<|Delta|` for every odd prime, Sections 2-3 prove that none is contained in
the union of the `m` halves. The actual support after orientations and
cancellations is a subset of that union, so none is contained there either.

This removes the two explicit non-rectangle weight-`|Delta|` families from
the puncture obstruction. It does not prove the conjectural row-code minimum
distance `ph`, classify all words through weight `|Delta|`, exclude every
other row-code support from the structured Mobius set, prove punctured
surjectivity, or close residual (ii).

## 5. Replay-only checks

The implementation checks all independent `(L,M)` and all nonzero `j` only at
`p=3,7`. It verifies injectivity of the difference-class map, the exact hit
criterion (3), and the unique identity-graph intersection. These are
fail-when-wrong formula replays, not theorem evidence and not a prime census.

Reproduction:

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_symmetric_halved_mobius_explicit_words.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_halved_mobius_explicit_words.py
