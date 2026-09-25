# Boolean phase range as a bipartite cut norm

2026-09-25. Deterministic all-orders theorem. The original convergence
problem remains OPEN.

Let B be a complete symmetric zero-diagonal signing of order k. Write

    P(B)=max_x Q_B(x),   N(B)=-min_x Q_B(x),

and define the phase range

    R(B)=P(B)+N(B)=max_x Q_B(x)-min_x Q_B(x).

For a bipartition U,V of the vertices, let B[U,V] be the rectangular cross
block and define its Boolean bilinear norm

    beta(U,V)=max_(u in {+-1}^U, v in {+-1}^V) u^T B[U,V] v.

Then

    boxed:
    R(B)=2 max_(U disjoint union V=[k]) beta(U,V).          (1)

## Proof of the exact identity

For Boolean x,y define

    U={i:x_i!=y_i},   V={i:x_i=y_i}.

Write a_i=x_i on U and b_j=x_j=y_j on V. Internal pair products on U and on
V are unchanged when x is replaced by y, while every U--V product changes
sign. Hence

    Q_B(x)-Q_B(y)
      =2 a^T B[U,V] b.                                     (2)

Taking x to maximize Q_B and y to minimize Q_B gives

    R(B)<=2 max_(U,V) beta(U,V).

Conversely, given U,V and Boolean a,b, form x by x_U=a,x_V=b and y by
y_U=-a,y_V=b. Then (2) holds, so

    R(B)>=Q_B(x)-Q_B(y)=2a^T B[U,V]b.

Maximizing proves (1).

## A universal k^(3/2) range lower bound

Fix |U|=a, |V|=b with a+b=k. For uniform random v on V, choose

    u_i=sign((B[U,V]v)_i).

Then

    u^T B[U,V]v
      =sum_(i in U) |(B[U,V]v)_i|.

Every row is a length-b sign vector, so each row sum has the law of
eps_1+...+eps_b. Therefore, with

    mu_b=E|eps_1+...+eps_b|,

there exists v such that

    beta(U,V)>=a mu_b.

Using (1),

    boxed:
    P(B)+N(B)
      >= R_k^*:=2 max_(1<=b<=k-1) (k-b) mu_b.              (3)

Asymptotically mu_b=sqrt(2b/pi)+o(sqrt b). The profile
(1-lambda)sqrt(lambda) is maximized at lambda=1/3, so

    boxed:
    R_k^*
      =[4/3 sqrt(2/(3pi))+o(1)] k^(3/2)                    (4)

with coefficient

    4/3 sqrt(2/(3pi)) = 0.6142... .

Thus the optimal elementary phase-range split is again response:seed = 2:1,
but here the conclusion is a RANGE bound, not the product bound.

## New cross-order consequence

Let A be a complete signing of order N with F=Phi(A). Delete a principal
set T of size t>=2 and retain S. Put B=A[T]. The strict restriction identity
gives

    Phi(A[S])<=F-min(P(B),N(B)).                            (5)

Since max(P(B),N(B))=Phi(B)<=F, (3) gives

    min(P(B),N(B))
      =P(B)+N(B)-Phi(B)
      >=R_t^*-F.                                           (6)

Therefore

    boxed:
    Phi(A[S])<=2F-R_t^*.                                   (7)

whenever the right-hand range payment is positive; in general combine (7)
with the independent half-per-vertex and product payments.

Applying (7) to an exact order-N minimizer and using
m_(N-t)<=Phi(A[S]) gives

    boxed:
    2m_N >= m_(N-t)+R_t^*.                                 (8)

Equivalently,

    boxed:
    m_N >= (m_(N-t)+R_t^*)/2.                              (9)

At proportional deletion t=lambda N, this is a genuine N^(3/2)-scale
cross-order inequality:

    alpha_N
      >= (1/2)(1-lambda)^(3/2) alpha_(N-t)
       + [2/3 sqrt(2/(3pi))+o(1)] lambda^(3/2).             (10)

The coefficient in the second term is

    2/3 sqrt(2/(3pi)) = 0.3071... .

Equation (1) is the main structural gain: the positive/negative phase range
of a quadratic signing is exactly twice a maximum Boolean bilinear cut norm.
This creates a direct bridge from the original minimax problem to
rectangular discrepancy tools, while (8) is the corresponding unconditional
cross-order consequence.
