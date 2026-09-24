# Geometry of the one-flip-stable extension skeleton

2026-09-23. Exact finite structure for the stable states isolated by
NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md. It gives a Hamming
packing rule and a local entropy bound in terms of light coordinates.
The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing, fix an orientation
sigma in {+1,-1}, and write

    w_i^sigma(x)=sigma x_i(Ax)_i.

A Boolean state x is sigma-stable when w_i^sigma(x)>=0 for every i.

## 1. Pair geometry

Let x and y be two sigma-stable states, and let

    T={i:x_i!=y_i},    t=|T|.

For every i in T,

    w_i^sigma(x) <= 2(t-1),                                 (1)

and symmetrically

    w_i^sigma(y) <= 2(t-1).                                 (2)

More precisely,

    sum_(i in T) w_i^sigma(x)
       <= 4 sigma Q_(A[T])(x_T),                            (3)

so the induced oriented energy on the disagreement set is nonnegative.

Proof. Since y is obtained from x by flipping exactly T, for i in T,

    sigma y_i(Ay)_i
      = -w_i^sigma(x)
        +2 sum_(j in T, j!=i) sigma A_ij x_i x_j.           (4)

The left side is nonnegative by stability of y. Each summand on the
right is +/-1, which proves (1). Summing (4) over i in T gives (3),
because every internal pair is counted twice. Equations (2) and its
summed analogue follow by symmetry.

## 2. Stability margin gives Hamming packing

Define the strict stability margin

    h_sigma(x)=min_i w_i^sigma(x).

If h_sigma(x)>0, every distinct sigma-stable y satisfies

    d_H(x,y) >= 1 + h_sigma(x)/2.                           (5)

In integral form,

    d_H(x,y) >= 1 + ceil(h_sigma(x)/2).                     (6)

Thus stable states with a macroscopic local-field margin are automatically
Hamming-separated. At even n, every local field is odd, so every stable
state is strict and no two distinct stable states are Hamming neighbours.

## 3. Light-coordinate localization

For h>=0 define the light set

    L_x(h)={i:w_i^sigma(x)<=h}.

Equation (1) gives the exact containment

    T subseteq L_x(2(t-1)).                                 (7)

Consequently the number N_x(r) of sigma-stable states y with
1<=d_H(x,y)<=r obeys

    N_x(r)
      <= sum_(j=1)^r binom(|L_x(2(r-1))|,j).                (8)

This is a deterministic local entropy bound for the stable skeleton.
It uses no optimizer assumption, spectral cap, random model, or
asymptotics.

## 4. Why this strengthens the one-vertex reduction

The one-vertex extension norm is exactly the maximum over the two oriented
stable skeletons. Equation (8) says that the nearby portion of each skeleton
is supported entirely on the light coordinates of any chosen stable state.
Therefore a random-row or discrepancy argument need not treat nearby stable
states as arbitrary independent constraints: their Hamming supports are
confined by the local-field profile of the center.

This does not yet bound the global number of stable states. It supplies an
exact geometric restriction on the state family that remains after the
previous all-unstable-state pruning.
