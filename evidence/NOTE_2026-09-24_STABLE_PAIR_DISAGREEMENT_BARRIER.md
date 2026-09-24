# Stable-pair disagreement energy barrier from strict restriction

2026-09-24. Deterministic structural theorem for the exact near-edge
skeleton. The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N and put

    F=Phi(A).

Fix one orientation sigma in {+1,-1}. For Boolean states x,y, write

    D_x = F-sigma Q_A(x),
    D_y = F-sigma Q_A(y).

Let

    T={i:x_i!=y_i},   t=|T|,   S=T^c,

and define the oriented internal disagreement energy

    q_T = sigma Q_(A[T])(x_T).

Because y=-x on T, the pair products inside T are unchanged, so the same
q_T is obtained from y_T.

For the induced signing A[T], write

    P_T=max_z Q_(A[T])(z),
    N_T=-min_z Q_(A[T])(z),
    ell_T=min(P_T,N_T).

## 1. Exact pair barrier

The two states have the same restriction to S and opposite cross terms.
Hence, if

    q_S=sigma Q_(A[S])(x_S),

then adding their two oriented energies gives

    q_S + q_T = F-(D_x+D_y)/2.                              (1)

The strict principal-restriction identity gives

    Phi(A[S]) <= F-ell_T.                                   (2)

Since q_S<=Phi(A[S]), combining (1)--(2) yields

    boxed:
    q_T + (D_x+D_y)/2 >= ell_T.                             (3)

Equivalently,

    boxed:
    q_T >= ell_T-(D_x+D_y)/2.                               (4)

This statement does NOT require x or y to be stable.

Thus two same-orientation states can both lie near the Boolean edge only if
their disagreement set carries correspondingly positive oriented internal
energy.

Using the exact two-sided product constant from the September-24 strict
restriction theorem,

    ell_T >= Delta_t/min(F,binom(t,2)),                     (5)

where

    Delta_t = a^2 mu_b^2/4,
    a=ceil(t/2), b=floor(t/2),

we obtain the explicit finite bound

    q_T >= Delta_t/min(F,binom(t,2))
           -(D_x+D_y)/2.                                    (6)

The simpler universal one-sided extremum bound
ell_T>=floor(t/2) also gives

    q_T >= floor(t/2)-(D_x+D_y)/2.                          (7)

## 2. Stable-pair local-field mass

Now assume x and y are both sigma-stable. Define

    w_i^sigma(x)=sigma x_i(Ax)_i,
    W_T(x)=sum_(i in T) w_i^sigma(x).

Flipping exactly T changes oriented energy by

    sigma Q_A(y)
      = sigma Q_A(x)-2W_T(x)+4q_T.                          (8)

Therefore

    W_T(x)=2q_T+(D_y-D_x)/2,                                (9)
    W_T(y)=2q_T+(D_x-D_y)/2.                               (10)

Insert (3). This gives the asymmetric exact lower bounds

    boxed:
    W_T(x) >= 2 ell_T-(3D_x+D_y)/2,                         (11)

    boxed:
    W_T(y) >= 2 ell_T-(D_x+3D_y)/2.                         (12)

In particular, if D_x,D_y<=d,

    boxed:
    W_T(x), W_T(y) >= 2 ell_T-2d.                           (13)

For two exact same-orientation maximizers,

    boxed:
    q_T>=ell_T,
    W_T(x)=W_T(y)=2q_T>=2ell_T.                             (14)

This removes the parity limitation of the earlier local-field argument:
the disagreement set of two exact maximizers carries positive internal
energy and positive local-field mass at every order.

## 3. Macroscopic consequence

Suppose t=lambda N with fixed 0<lambda<=1 and F=alpha N^(3/2).
The exact product constant has

    Delta_t=(1/(16pi)+o(1))t^3.

Whenever binom(t,2)>=F (true for fixed lambda and large N), (5) gives

    ell_T >=
      [lambda^3/(16pi alpha)+o(1)] N^(3/2).                 (15)

Hence two sigma-stable states with D_x,D_y=O(N) and Hamming distance
lambda N satisfy

    W_T(x)/t, W_T(y)/t
      >= [lambda^2/(8pi alpha)+o(1)] sqrt(N).               (16)

So a macroscopically separated pair of active near-edge stable states cannot
disagree mainly on light local-field coordinates: the AVERAGE stability
margin across their disagreement set is forced to be order sqrt(N).

This directly complements the earlier stable-skeleton light-set theorem,
which confined nearby stable states to light coordinates. The two results
now squeeze the skeleton from opposite directions:

- small Hamming moves can occur only through light coordinates;
- macroscopic Hamming separation of two near-edge stable states forces
  order-sqrt(N) average field mass on the disagreement set.

The theorem is finite and deterministic. No Gaussian approximation,
operator cap, or finite census is used.
