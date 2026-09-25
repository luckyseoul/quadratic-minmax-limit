# Half-per-deleted-vertex principal restriction theorem

2026-09-24. Deterministic all-orders cross-order theorem. This combines the
same-day strict principal-restriction identity with the elementary one-sided
extremum floor already isolated in the stable-pair note. The resulting
cross-order inequality was not previously stated in the repository.

Let A be any complete symmetric zero-diagonal signing of order N, let S be a
principal set of size n=N-t, let T=S^c, and put F=Phi(A). Assume t>=2.

For B=A[T], write

    P(B)=max_y Q_B(y),   N(B)=-min_y Q_B(y).

The strict restriction identity proved earlier today gives

    Phi(A[S]) <= F-min(P(B),N(B)).                           (1)

The key point is that every complete signing B of order t satisfies

    min(P(B),N(B)) >= floor(t/2).                           (2)

Proof of (2). If t is even, take a global maximizer y of Q_B. One-spin
optimality gives

    y_i(By)_i >= 0

for every i. Since t-1 is odd, each local field (By)_i is an odd integer, so
every y_i(By)_i is at least 1. Therefore

    2P(B)=sum_i y_i(By)_i >= t,

hence P(B)>=t/2. Applying the same argument to -B gives N(B)>=t/2.

If t is odd, restrict B to any t-1 vertices. The even case gives a positive
maximum at least (t-1)/2 on that principal block. Extend its maximizing spin
assignment to the last vertex, choosing the better of the two signs; the
linear incident term contributes in absolute value, so the full positive
maximum is at least (t-1)/2. Again apply this to -B for the negative side.
Thus (2) holds for all t>=2.

Combining (1) and (2) gives the universal host-free restriction theorem

    boxed:
    Phi(A[S]) <= Phi(A)-floor(t/2).                          (3)

This holds for EVERY principal S, with no optimizer assumption.

Applying (3) to an exact order-N minimizer yields

    boxed:
    m_N-m_(N-t) >= floor(t/2),   t>=2.                      (4)

Equivalently, for every N>=n+2,

    boxed:
    m_N-m_n >= floor((N-n)/2).                              (5)

Because each m_k has the parity of binom(k,2), (5) may be rounded further to
the least integer at least floor((N-n)/2) with the required endpoint parity.

Asymptotically, for any t=t(N) tending to infinity,

    m_N-m_(N-t) >= (1/2-o(1))t.                             (6)

This improves the same-day universal host-free slopes

    1/32 = 0.03125,
    (3-sqrt(7))/2 = 0.177124...,

to the exact coefficient 1/2.

The fourth-moment completion theorem remains stronger when its retained
cross-block term is large:

    Phi(A[S]) <= F-c4 sqrt(binomial(t,2)+||A[T,S]x||_2^2),

but after discarding the cross term, (3) is the canonical universal bound.

The significance is cross-order rather than asymptotic normalization: every
two deleted vertices now force at least one full raw unit of Boolean norm
loss, uniformly over all complete signings and all principal subsets.
