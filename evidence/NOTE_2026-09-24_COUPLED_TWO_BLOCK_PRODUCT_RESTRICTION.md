# Coupled two-block product restriction

2026-09-24. Deterministic cross-order strengthening. It uses BOTH principal
blocks simultaneously, rather than charging only the deleted block. The
original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N, partition
the vertices into S,T with |S|=n>=2, |T|=t>=2, n+t=N, and put

    F=Phi(A).

For U=S,T define

    p_U=max Q_(A[U]),   q_U=-min Q_(A[U]).

Let Delta_k^* be the optimized product constant from
NOTE_2026-09-24_OPTIMIZED_PRODUCT_SPLIT.md:

    Delta_k^*
      =(1/4) max_(1<=b<=k-1) (k-b)^2 mu_b^2,

so every order-k complete signing satisfies

    p_U q_U >= Delta_k^*.

Set

    delta_n=Delta_n^*,   delta_t=Delta_t^*,

    A0=(F^2+delta_t-delta_n)/F,

and

    R_(n,t)(F)
      = [A0-sqrt(A0^2-4 delta_t)]/2.                        (1)

Feasibility of the actual four one-sided extrema guarantees that the square
root is real.

Then

    boxed:
    min(p_T,q_T) >= R_(n,t)(F),                             (2)

and therefore

    boxed:
    Phi(A[S]) <= F-R_(n,t)(F).                              (3)

## Proof

The disjoint-phase inequalities give

    p_S+p_T<=F,
    q_S+q_T<=F.                                             (4)

Let c=min(p_T,q_T). Interchange the positive and negative labels if
necessary so c=p_T<=q_T. Since

    p_T q_T>=delta_t,

we have

    q_T>=delta_t/c.                                         (5)

The retained block satisfies

    p_S q_S>=delta_n,

while (4) gives

    p_S<=F-c,
    q_S<=F-q_T<=F-delta_t/c.

Hence

    (F-c)(F-delta_t/c)>=delta_n.                            (6)

After multiplying by c and rearranging,

    c^2-A0 c+delta_t<=0.                                    (7)

Therefore c lies between the two roots of this quadratic, proving (2).
Equation (3) follows from

    Phi(A[S])=max(p_S,q_S)
      <=F-min(p_T,q_T).

## Strict improvement over the one-block charge

The previous host-scale restriction used only T and gave

    min(p_T,q_T)>=delta_t/F.                                (8)

Here both product floors are positive. Also Cauchy and (4) give

    sqrt(delta_n)+sqrt(delta_t)
      <=sqrt((p_S+p_T)(q_S+q_T))
      <=F,

so delta_n+delta_t<F^2. Thus delta_t/F lies strictly to the left of the
quadratic vertex in (7). Evaluating the left side of (7) at delta_t/F gives

    delta_n delta_t/F^2 >0.

Consequently

    boxed:
    R_(n,t)(F) > delta_t/F.                                 (9)

So the coupled theorem strictly improves the earlier host-scale bound at
every nontrivial two-block split.

## Balanced form

For N=2k, n=t=k, put delta=Delta_k^*. Then A0=F and

    boxed:
    R_(k,k)(F)
      =[F-sqrt(F^2-4delta)]/2.                              (10)

Hence every k-vertex principal half satisfies

    boxed:
    Phi(A[S])
      <= [F+sqrt(F^2-4Delta_k^*)]/2.                        (11)

For an exact order-2k minimizer,

    boxed:
    m_(2k)-m_k
      >= [m_(2k)-sqrt(m_(2k)^2-4Delta_k^*)]/2.              (12)

This is a strict nonlinear improvement over
Delta_k^*/m_(2k).

## Asymptotic input

Since

    Delta_k^*=[2/(27pi)+o(1)]k^3,

the balanced improvement is on the full k^(3/2) scale. If
F=C k^(3/2), its normalized loss coefficient is

    [C-sqrt(C^2-8/(27pi))]/2,                               (13)

strictly larger than 2/(27pi C).

The theorem is finite and deterministic and uses no optimizer
classification, spectral assumption, Gaussian approximation, or census.
