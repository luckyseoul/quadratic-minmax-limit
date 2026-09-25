# Pair-factorization theorem for stable near-edge states

2026-09-24. Exact deterministic structure of the stable near-edge skeleton.
The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N and put

    F=Phi(A).

Fix one orientation sigma in {+1,-1}. Let x,y be sigma-stable Boolean states
with deficits

    D_x=F-sigma Q_A(x),   D_y=F-sigma Q_A(y).

Let

    T={i:x_i!=y_i},   S=T^c,

and write

    q_T=sigma Q_(A[T])(x_T),
    q_S=sigma Q_(A[S])(x_S).

Because y=-x on T and y=x on S, the internal energies are the same for x
and y on both blocks.

For a principal block U define its one-sided sigma maximum

    P_sigma(U)=max_z sigma Q_(A[U])(z).

## Theorem 1: exact near-optimal factorization

Put d=(D_x+D_y)/2. Then

    boxed:
    0 <= P_sigma(T)-q_T <= d,                               (1)

    boxed:
    0 <= P_sigma(S)-q_S <= d.                               (2)

Moreover x_T is sigma-stable for A[T], and x_S is sigma-stable for A[S].
The same is true for the corresponding restrictions of y.

Thus every pair of same-orientation stable d-near-edge states splits the
ambient signing into TWO induced subproblems on which the inherited states
are again stable and d-near their own one-sided optima.

Proof. Let C be the signed cross energy in orientation sigma for x:

    C=sigma sum_(i in T,j in S) A_ij x_i x_j.

Then

    sigma Q_A(x)=q_T+q_S+C,
    sigma Q_A(y)=q_T+q_S-C.

Averaging gives

    q_T+q_S=F-d.                                            (3)

Now freeze x_S and complete the T coordinates by a maximizer z of
P_sigma(T). Pair z with -z. The T-internal energy is unchanged while the
cross term changes sign, so one of the two completions has oriented full
energy at least

    q_S+P_sigma(T).

Since every full oriented energy is at most F,

    q_S+P_sigma(T)<=F.

Subtract (3):

    P_sigma(T)-q_T<=d.

The lower inequality is admissibility of x_T. The S statement is symmetric.

For stability, fix i in T. Write h_T for the field at i from T and h_S for
the field from S. Stability of x and y gives

    sigma x_i(h_T+h_S)>=0,
    sigma x_i(h_T-h_S)>=0,

because both the spin and the T-internal field reverse when passing from x
to y, while the S-field does not. Adding yields

    sigma x_i h_T>=0.

So x_T is sigma-stable in A[T]. The proof on S is identical. QED.

## Theorem 2: the cross cut is exactly the deficit imbalance

With the same notation,

    boxed:
    C=(D_y-D_x)/2.                                         (4)

Hence if D_x,D_y<=d_0,

    |C|<=d_0/2.                                             (5)

So a pair of active near-edge stable states does not merely have good
internal restrictions: their disagreement cut has raw signed discrepancy
only O(d_0).

For the Komlos active window d_0=O(N), this improves the generic
O(N^(3/2)) cut scale by a factor sqrt(N).

## Corollary 3: exact maximizers factorize exactly

If x and y are same-orientation exact maximizers, D_x=D_y=0. Then

    boxed:
    C=0,                                                    (6)

    boxed:
    q_T=P_sigma(T),   q_S=P_sigma(S),                       (7)

and the inherited states on T and S are stable exact one-sided maximizers
of the two principal blocks.

Thus two exact maximizers can coexist only across a ZERO signed cut, and
their common/opposite coordinate blocks inherit exact extremality.

There is an immediate parity consequence. If N is even, a zero cut has
t(N-t) sign terms. This number must be even. Since N is even, t and N-t
have the same parity, so t cannot be odd. Therefore

    boxed:
    at even order, any two same-orientation exact maximizers
    have even Hamming distance.                             (8)

The antipodal quotient preserves this parity because N is even.

## Corollary 4: self-similarity of the active stable skeleton

For any fixed d_0, the family

    S_sigma(d_0)
      ={x: x sigma-stable, F-sigma Q_A(x)<=d_0}

has the following hereditary pair property:

For every x,y in S_sigma(d_0), both the agreement block and disagreement
block inherit sigma-stable states lying within d_0 of the corresponding
one-sided principal optimum, while their cross cut has absolute signed
energy at most d_0/2.

This is substantially stronger than a Hamming separation statement. It says
that a large near-edge stable family would force the signing to contain a
large family of principal bipartitions that are simultaneously

1. internally stable on both sides;
2. near-optimal on both sides; and
3. almost perfectly balanced across the cut.

That is the new structural object to count in the Komlos obstruction.
