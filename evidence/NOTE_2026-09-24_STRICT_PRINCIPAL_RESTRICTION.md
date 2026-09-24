# Quantitative strict principal restriction from two-sided Boolean energy

2026-09-24. This is a deterministic cross-order theorem for the original
quadratic minimax problem. It sharpens ordinary restriction monotonicity.
The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N and put

    F=Phi(A).

Let S be ANY principal set of size n=N-t, with t>=2, and T=S^c.

For an order-t complete signing B define

    P(B)=max_x Q_B(x),   N(B)=-min_x Q_B(x).

The exact product argument from
NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md actually retained

    P(B) N(B) >= Delta_t,                                  (1)

where

    Delta_t = a^2 mu_b^2/4,
    a=ceil(t/2), b=floor(t/2),
    mu_b=E|eps_1+...+eps_b|.                               (2)

The previous note replaced this exact constant by t^3/96. Here the exact
quantity is kept.

The main result is

    THEOREM.
    For every S with |S|=N-t and t>=2,

        Phi(A[S]) <= F - Delta_t/F.                        (3)

In particular every proper restriction deleting at least two vertices has
strictly smaller Boolean norm.

## 1. One-sided extension identity

Write

    P_T=max_y Q_(A[T])(y),
    N_T=-min_y Q_(A[T])(y).

Fix x on S. For y on T,

    Q_A(x,y)=Q_(A[S])(x)+Q_(A[T])(y)+L_x(y),

where L_x is linear in y. Pair y with -y. Their T-internal energies agree
and their cross terms have opposite sign, so

    max{Q_A(x,y),Q_A(x,-y)}
      = Q_(A[S])(x)+Q_(A[T])(y)+|L_x(y)|
      >= Q_(A[S])(x)+Q_(A[T])(y).                          (4)

Taking y maximizing Q_(A[T]) and then x maximizing Q_(A[S]) gives

    max_z Q_A(z) >= max_x Q_(A[S])(x)+P_T.                 (5)

Apply the same argument to -A:

    -min_z Q_A(z) >= -min_x Q_(A[S])(x)+N_T.               (6)

Since both one-sided full extrema are at most F,

    Phi(A[S]) <= F-min(P_T,N_T).                           (7)

## 2. The complement cannot be one-sided-flat

For B=A[T], principal restriction and the trivial complete-graph ceiling give

    max(P_T,N_T)=Phi(B)<=min(F,binom(t,2)).                 (8)

The exact product bound (1) therefore implies the STRONGER estimate

    min(P_T,N_T)
      >= P_T N_T/max(P_T,N_T)
      >= Delta_t/min(F,binom(t,2)).                         (9)

Combining (7) and (9) proves

    Phi(A[S])
      <= F-Delta_t/min(F,binom(t,2))
      <= F-Delta_t/F.                                      (9a)

Thus (3) is only the host-scale corollary; (9a) is the canonical form.

This proof is finite, requires no optimizer assumption, and applies to every
choice of the retained principal set S.

## 3. Explicit form of Delta_t

For b>=1,

    mu_b
      = b 2^(-(b-1)) binom(b-1,floor((b-1)/2)).            (10)

Thus Delta_t is completely explicit.

The sharp p=1 Khintchine inequality for equal coefficients, or directly
Wallis' central-binomial estimate in (10), gives

    mu_b^2 >= b/2.

Since a^2 b>=t^3/8,

    Delta_t >= t^3/64.                                     (11)

Hence (9a) gives TWO clean finite corollaries:

    boxed: Phi(A[S]) <= F - t^3/(64F),                     (12)

and, because binom(t,2)<=t^2/2,

    boxed: Phi(A[S]) <= F - t^2/[32(t-1)]
                        <= F - t/32.                       (12a)

The second bound is independent of the host norm. Thus deleting any t>=2
vertices from a complete signing forces a raw Boolean-norm loss at least
t/32 (and in fact at least t^2/[32(t-1)]).

This improves the 1/96 constant implicit in the September 23 interior-gap
argument. More precisely, as t->infinity,

    mu_b^2=(2/pi)b+O(1),
    Delta_t=(1/(16pi)+o(1)) t^3.                           (13)

So the asymptotic restriction constant is 1/(16pi), almost twice 1/96.

## 4. New cross-order inequalities for m_n

Take A to be an exact order-N minimizer. Since every n-vertex principal
restriction is a legal order-n signing, (9a) gives, for t=N-n>=2,

    boxed:
    m_n <= m_N-Delta_t/min(m_N,binom(t,2))
        <= m_N-t^2/[32(t-1)]
        <= m_N-t/32.                                       (14)

The host-scale form retained for algebra below is

    m_n <= m_N-Delta_t/m_N.                                (14a)

Equivalently, (14a) gives

    m_N^2-m_n m_N-Delta_t >=0,

so

    boxed:
    m_N >= [m_n+sqrt(m_n^2+4 Delta_t)]/2.                  (15)

Using (11),

    m_N >= [m_n+sqrt(m_n^2+t^3/16)]/2.                    (16)

This is a quantitative strict-growth theorem across every gap of at least
two orders.

In particular,

    boxed: m_N-m_n >= (N-n)/32 whenever N>=n+2.            (17)

Since the difference is integral,

    boxed: m_N-m_n >= ceil((N-n)/32).                      (17a)

The parity of each endpoint is fixed by binom(k,2), so (17a) may be rounded
further to the least integer at least (N-n)/32 with the required endpoint
parity. In particular m_N>m_n whenever N>=n+2.

Using the elementary all-orders upper bound
m_N<=sqrt(log 2) N^(3/2), (14) also gives the explicit separation

    m_N-m_n
      >= Delta_t/[sqrt(log 2) N^(3/2)]
      >= t^3/[64 sqrt(log 2) N^(3/2)].                     (18)

For proportional deletion t=lambda N, the exact asymptotic form is

    m_N-m_((1-lambda)N)
      >= [lambda^3/(16 pi sqrt(log 2))+o(1)] N^(3/2).       (19)

No averaging over principal subsets is needed.

## 5. Scope

This is genuine cross-order control, but its gain is cubic in the deleted
fraction. For t=o(N) it is much smaller than the linear-scale control that
would directly force convergence of alpha_n. It therefore does not settle
the original limit by itself.

Its useful new content is that restriction is not merely monotone: every
two-or-more-vertex deletion has a quantified mandatory norm loss.  The
complement ceiling upgrades the loss to at least one thirty-second of the
number of deleted vertices at ALL scales, while the exact product constant
gives the stronger proportional-scale n^(3/2) separation when t is macroscopic.
