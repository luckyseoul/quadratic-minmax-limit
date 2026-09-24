# Pointwise deficit controls best-response disagreement

2026-09-23. This sharpens the same-order mean-update information into a
pointwise, distribution-free inequality. It is valid for every complete
signing and every Boolean state. The original convergence problem remains
OPEN.

Let A be a complete symmetric zero-diagonal signing of order n and put

    M=A/sqrt(n),
    alpha=Phi(A)/n^(3/2),
    F=alpha*n.

Fix an orientation sigma in {+1,-1} and a Boolean state x. Write

    q=Q_(sigma M)(x),
    D=F-q >= 0,
    Y=sign(sigma M x),

with sign(0)=x_i at a zero field, and define the total weighted
best-response disagreement

    R=sum_(i:Y_i!=x_i) |(M x)_i|.

Then

    R <= D + sqrt(2 F D).                                  (1)

This is sharp with respect to the entire one-parameter mean-update family:
(1) is exactly the result of optimizing that family over its interpolation
parameter.

## 1. Exact derivation

For 0<=p<=1 put

    z=(1-p)x+pY.

Since z belongs to [-1,1]^n and the quadratic form is multilinear,

    Q_(sigma M)(z) <= F.                                   (2)

The exact expansion is

    Q_(sigma M)(z)
      =(1-p)^2 q
        +p(1-p)||M x||_1
        +p^2 Q_(sigma M)(Y).                               (3)

Also Q_(sigma M)(Y)>=-F. On a changed coordinate,
x_i (sigma M x)_i=-|(M x)_i|, while unchanged coordinates contribute
positively. Hence

    ||M x||_1-2q = 2R.                                     (4)

Combining (2)-(4) gives, for every 0<p<1,

    R <= F*p/(1-p) + D*(1+p)/(2p).                         (5)

The right side is minimized at

    p_* = sqrt(D/2)/(sqrt(F)+sqrt(D/2)),                    (6)

with the D=0 case interpreted as the p->0 limit. Substitution gives (1).
In particular D=0 forces R=0: every exact positive maximizer is already
a synchronous best-response fixed point, up to zero-field ties.

Conversely, the right side of (1) is the minimum of the right side of
(5), so (1) is equivalent to the full family of inequalities (5).

## 2. Near-edge scale

If D=o(F), then

    R/F <= D/F + sqrt(2D/F) -> 0.                           (7)

For the original asymptotic regime F=Theta(n), a state only O(sqrt(n))
below the oriented Boolean edge satisfies

    R=O(n^(3/4)).                                           (8)

Thus a near-edge state cannot carry Theta(n) weighted best-response
disagreement at that scale. This conclusion is pointwise; it needs no
law, operator cap, Gaussian approximation, or spectral regularization.

## 3. Two-phase distribution-free closure

Let X_sigma be arbitrary Boolean laws in the two orientations and set

    e=(1/(2n)) sum_sigma E Q_(sigma M)(X_sigma),
    f=(1/(2n)) sum_sigma E ||M X_sigma||_1.

Averaging (1) over the equal two-phase mixture and using concavity of
the square root gives

    (f-2e)/2
      <= (alpha-e) + sqrt(2 alpha (alpha-e)).               (9)

Equivalently,

    f-2e
      <= 2(alpha-e)+2 sqrt(2 alpha (alpha-e)).              (10)

Solving (9) for alpha recovers the exact optimized mean-update envelope

    alpha
      >= e-f/2 + (1/2)*sqrt(f^2+(f-2e)^2).                 (11)

So (9) loses nothing relative to optimizing the old common-p inequality
when only the phase-averaged e and f are retained, while (1) adds the new
statewise information.

The earlier adaptive-best-response closure

    alpha >= e+(f-2e)^2/(16 alpha)

is a weaker corollary at the averaged level. It remains valid, but (9)-(11)
are the canonical closure.

## 4. Scope

This is not a cross-order comparison and does not select a one-vertex
extension row. Its positive contribution is a sharp pointwise edge-stability
law at the O(sqrt(n)) deficit scale that is directly relevant to the
near-extremal states appearing in one-vertex extension and repeated-update
arguments.
