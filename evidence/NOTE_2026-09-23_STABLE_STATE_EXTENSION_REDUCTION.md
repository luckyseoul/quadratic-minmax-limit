# One-vertex extension needs only one-flip-stable states

2026-09-23. Exact finite reduction for the original one-vertex extension
problem. This removes every oriented state with an improving coordinate
from the extension maximum, before any probabilistic or discrepancy
estimate is applied. The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order n, let

    F=Phi(A),

and let a in {+1,-1}^n be a proposed incident row for one new vertex.
The extended signing has norm

    Phi(A+a)
      = max_x ( |Q_A(x)| + |a.x| )
      = max_(sigma=+/-1,x) [ sigma Q_A(x)+|a.x| ].          (1)

For sigma in {+1,-1}, call x sigma-stable if

    sigma x_i (A x)_i >= 0    for every i.                  (2)

Then, for EVERY row a,

    Phi(A+a)
      = max_(sigma=+/-1, x sigma-stable)
            [ sigma Q_A(x)+|a.x| ].                         (3)

Thus one-vertex extension is an exact discrepancy problem on the
one-flip-stable states, not on all 2^n Boolean states.

## 1. Monotone repair of the extension score

Fix sigma and x. If x is not sigma-stable, choose i with

    w_i = sigma x_i(Ax)_i < 0.

Flipping coordinate i gives

    sigma Q_A(x^(i)) - sigma Q_A(x) = -2 w_i.              (4)

The local field (Ax)_i is an integer sum of n-1 signs. Hence w_i is an
integer, so w_i<=-1 and the oriented energy gain in (4) is at least 2.

At the same flip,

    a.x^(i) = a.x - 2 a_i x_i,

therefore

    |a.x^(i)| - |a.x| >= -2.                               (5)

Combining (4)-(5),

    [sigma Q_A(x^(i))+|a.x^(i)|]
      -[sigma Q_A(x)+|a.x|] >= -2w_i-2 >= 0.               (6)

Repeat while a negative oriented local field remains. Each strict energy
improvement is at least two and sigma Q_A is bounded above by F, so the
process terminates at a sigma-stable state. The extension score never
decreases. This proves (3).

Zero local fields need not be flipped. They already satisfy stability.

## 2. Stable-state form of the deficit objective

For a sigma-stable state define

    D_sigma(x)=F-sigma Q_A(x) >= 0.

Then (3) is exactly

    Phi(A+a)-F
      = max_(sigma,x sigma-stable)
            [ |a.x|-D_sigma(x) ].                           (7)

Moreover stability gives the exact local-field identity

    ||Ax||_1 = 2 sigma Q_A(x),                              (8)

because all terms sigma x_i(Ax)_i are nonnegative and their sum is
2 sigma Q_A(x).

So the dangerous extension states simultaneously have:
1. one-flip stability,
2. small oriented deficit,
3. a fixed total absolute local-field mass determined by their energy.

These facts hold before choosing a.

## 3. Strictly sharper random-row criterion

Choose the new row a uniformly at random. For each fixed Boolean x,
a.x has the law S_n of a sum of n independent uniform signs.

Let H_sigma contain one representative from every antipodal pair of
sigma-stable states. From (7) and a union bound,

    P_a{ Phi(A+a)-F >= r }
      <= sum_(sigma=+/-1) sum_(x in H_sigma)
           P{ |S_n| >= D_sigma(x)+r }.                      (9)

Consequently, if

    sum_(sigma=+/-1) sum_(x in H_sigma)
       P{ |S_n| >= D_sigma(x)+r } < 1,                      (10)

then there exists an actual sign row a with

    Phi(A+a)-F < r.                                         (11)

This strictly sharpens the earlier exact-tail criterion that summed over
all Boolean states: every unstable state has been removed exactly, not
estimated or discarded heuristically.

For an exact order-n minimizer, taking

    r = F*((1+1/n)^(3/2)-1)

in (10) yields alpha_(n+1)<alpha_n.

## 4. Relation to the pointwise deficit-disagreement theorem

The separate pointwise theorem
NOTE_2026-09-23_POINTWISE_DEFICIT_DISAGREEMENT.md shows that an arbitrary
near-edge state can be repaired toward one-flip stability while controlling
weighted disagreement. The present theorem is stronger for the extension
objective itself: because an improving flip raises oriented energy by at
least the maximum possible loss 2 in |a.x|, the extension score can be
repaired monotonically all the way to a stable state.

This uses the integral complete-signing local fields essentially. For
general real edge weights an arbitrarily small improving energy step would
not dominate the possible change two in the new-row correlation.

## 5. Scope

No bound on the number of stable states is asserted. Equation (3) is,
however, an exact reduction of the cross-order one-vertex problem from the
entire Boolean cube to its oriented one-flip-stable skeleton, and (10) is a
strictly stronger sufficient extension test than the previous all-state
tail sum.
