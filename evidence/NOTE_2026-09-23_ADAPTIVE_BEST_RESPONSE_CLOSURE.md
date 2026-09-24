# A cap-free adaptive best-response closure inequality

2026-09-23. This is a same-order positive theorem for the original complete
Boolean quadratic form. It converts weighted best-response disagreement
directly into an energy gain, without an operator cap and without counting
the number of changed coordinates. The original convergence problem remains
OPEN.

Let M=A/sqrt(n) for a complete symmetric zero-diagonal signing A, and put

    alpha = Phi(A)/n^(3/2),

so every z in [-1,1]^n obeys |Q_M(z)|<=alpha n.

Let X be ANY random Boolean vector. Write

    F = M X,
    Y = sign(F),

with sign(0)=X_i at zero field, and let

    S={i:Y_i!=X_i},
    R=sum_(i in S) |F_i|.

Then there exists a randomized partial best-response update X' on the SAME
Boolean cube such that, conditionally on X,

    E[Q_M(X') | X] >= Q_M(X) + R^2/(4 alpha n).             (1)

The randomization is only a proof device; hence some deterministic Boolean
outcome in its support attains at least the same conditional expectation.

## 1. Exact partial-flip expansion

Flip each coordinate in S independently with probability p, leaving every
other coordinate unchanged. If Delta=Y-X, then the updated vector can be
written X'=X+xi*Delta coordinatewise. Since M has zero diagonal,

    E[Q_M(X') | X]
      = Q_M(X)
        +p [ ||MX||_1-X^T M X ]
        +(p^2/2) Delta^T M Delta.                           (2)

On S one has X_i F_i=-|F_i|, so

    ||MX||_1-X^T M X = 2R.                                 (3)

Also Delta_i=-2X_i on S and zero elsewhere, hence

    (1/2) Delta^T M Delta = 4 Q_(M[S])(X_S).                (4)

Principal restriction and multilinearity give

    |Q_(M[S])(X_S)| <= alpha n.                             (5)

Therefore

    E[Q_M(X') | X] >= Q_M(X)+2pR-4p^2 alpha n.             (6)

The bilinear polarization bound
|u^T M v|<=4 alpha n for Boolean u,v, together with
X^T M X>=-2 alpha n, gives R<=3 alpha n. Thus

    p_X = R/(4 alpha n) <= 3/4                              (7)

is admissible. Substitution in (6) proves (1).

The alpha=0 case cannot occur for a complete signing of order n>=2.

## 2. Distribution-free closure for two phases

For either phase s=+/- let X_s have ANY law on the Boolean cube and define

    F_s=s M X_s,
    e=(1/(2n)) sum_s E Q_(sM)(X_s),
    f=(1/(2n)) sum_s E ||M X_s||_1.

Let R_s be the corresponding changed-field mass. From (3),

    (1/2) sum_s E R_s / n = (f-2e)/2.                      (8)

Apply (1) separately in both phases. Since every updated Boolean state still
has oriented energy at most alpha n,

    alpha n
      >= E Q_(sM)(X_s)+E R_s^2/(4 alpha n).

Average over s, divide by n, and use Jensen over the phase/law mixture:

    alpha >= e + (f-2e)^2/(16 alpha).                       (9)

Equivalently,

    |f-2e| <= 4 sqrt(alpha(alpha-e)),                       (10)

and in fact f-2e>=0 by (8).

This is a cap-free posterior closure inequality. It applies to Gaussian,
Gibbs, or arbitrary post-update laws and does not require a joint CLT,
spectral regularization, or an estimate for the Hamming count of changed
coordinates.

## 3. Relation to the September 23 interior-gap term

The finite interior-gap theorem gives a correction proportional to the
cubed changed-coordinate fraction. That is strongest when disagreement is
spread over a macroscopic set. Equation (9) is complementary: it uses the
TOTAL changed local-field mass and remains nontrivial even if that mass is
concentrated on relatively few coordinates.

For the tilted paired Gaussian phases, inserting only the existing lower
inputs e(t),f(t) into (9) does not beat the stronger mean-update constant
already recorded in NOTE_2026-09-23_TILTED_PAIRED_FIELD_LOWER.md. Its value
here is structural: repeated-update arguments now have a distribution-free
relation between e and f without first proving a lower bound on eta.

In particular, whenever e approaches alpha along any sequence of actual
laws, (10) forces f-2e->0. Thus near-extremal source energy automatically
suppresses weighted best-response disagreement.

No cross-order comparison is claimed here.
