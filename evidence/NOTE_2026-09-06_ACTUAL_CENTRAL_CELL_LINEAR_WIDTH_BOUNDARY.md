# Actual central cells saturate the Gaussian linear-field width

2026-09-06. Analytic architecture boundary, not an impossibility theorem
for the original MO problem or for other Gaussian comparisons. No
mathematical program, scan, solver, signing census, or numerical run
was used. The states and fields below are actual, not formal profiles.

## 1. Actual objects and a large central set

Let K=[[A,B],[B^T,-A]] be any complete paired signing of order N=2n,
and let D>0 be diagonal with D+-K>=0. Write

    H_L=D_L^(-1/2) A D_L^(-1/2),
    H_R=D_R^(-1/2) A D_R^(-1/2),
    W=D_L^(-1/2) B D_R^(-1/2).

All three operator norms are at most one. For z=(x,y) in the Boolean
cube, retain ORIGINAL p=x^TAx, q=y^TAy, c=x^TBy, and weighted
p_D=x^TH_Lx, q_D=y^TH_Ry, c_D=x^TWy. Define the actual central set C by

    |p|,|q|,|c|<=n^(5/4),
    |p_D|,|q_D|,|c_D|<=n^(3/4).                          (1)

For uniform independent Boolean x,y, zero diagonals give

    E p^2=E q^2=2n(n-1),       E c^2=n^2,
    E p_D^2=2tr(H_L^2)<=2n,
    E q_D^2=2tr(H_R^2)<=2n,    E c_D^2=||W||_F^2<=n.

Chebyshev and a union bound consequently prove

                         Pr(z not in C)<=10n^(-1/2).    (2)

Thus C is nonempty for all sufficiently large n. It is invariant under
the simultaneous global sign change z->-z. These facts require no
trace cap, diagonal near-flatness, or source/conditional optimality.

## 2. The actual central union has asymptotically full iid field width

Fix 0<=k<=w<=1, with w>0 for now. Let g0~N(0,wnI_N). Its maximizing
Boolean state is z*=sign(g0), which is uniform. Its signs are independent
of all coordinate magnitudes, hence of L0=sum_i |(g0)_i|. On z* in C,
max_C g0^Tz=L0. On the complementary event max_C g0^Tz>=0, because
C is nonempty and globally sign-symmetric. Therefore

    E max_C g0^Tz >= E[L0 1_(z* in C)]
        >=(1-10n^(-1/2)) 2sqrt(kappa w)n^(3/2),          (3)
    kappa=2/pi.

The reverse upper E max_C g0^Tz<=E L0=2sqrt(kappa w)n^(3/2) is exact.
Thus the central restrictions on all SIX actual energies do not remove
leading Gaussian linear-field width.

## 3. At least one genuine refined cell retains that full width

Intersect C with the weighted-shell theorem's partition by ORIGINAL
integer triple and weighted bins of side 1/n. There are at most

                         m<=(2n^2+1)^6                 (4)

nonempty resulting cells C_j. Choose each representative within its
FINAL intersection with C. For the common iid g0, every cell maximum
has Gaussian Lipschitz constant sqrt(wnN)=n sqrt(2w). Concentration
and the finite maximum bound give

    E max_C g0^Tz
       <=max_j E max_(C_j) g0^Tz+2n sqrt(w log m).        (5)

No independence between these maxima is assumed.
The representative theta_j=(p_D,j,q_D,j,c_D,j) supplies the actual
positive field covariance of the weighted-shell theorem:

    M_j=wnI+k[[q_D,j H_L,-c_D,j W],[-c_D,j W^T,p_D,j H_R]],
    M_j>=0,          (M_j)_ii=wn.

Because the representative belongs to C,

                         ||M_j-wnI||op<=2k n^(3/4).     (6)

Couple g_j=M_j^(1/2)G and g0=sqrt(wn)G with the same standard G.
Since scalar I commutes with M_j, the scalar square-root inequality
implies ||M_j^(1/2)-sqrt(wn)I||op<=sqrt(||M_j-wnI||op).
Cauchy--Schwarz then gives, uniformly over every cell,

    |E max_(C_j) g_j^Tz-E max_(C_j) g0^Tz|
                     <=2sqrt(2k)n^(11/8).               (7)

Indeed the expectation of the full-cube supremum of the coupling
difference is at most N times its operator square-root difference.
On the other hand, the actual constant diagonal of every M_j gives

    E max_(C_j) g_j^Tz<=sum_i E|(g_j)_i|
                                      =2sqrt(kappa w)n^(3/2).

Equations (3)-(7) prove the exact leading conclusion

    max_j E max_(C_j) g_j^Tz/(2n^(3/2)) ->sqrt(kappa w).  (8)

In particular the same limit holds when the maximum also includes
all other actual original/weighted cells: the constant-diagonal upper
holds for those cells too. For w=0, k=0 and (8) is immediate.

## 4. Consequence for the current linear-field upper, and its limits

For any fixed drift parameter |s|<=1, the ORIGINAL cell offset is
a_j=(p_j-q_j)/2+s c_j. On these central cells, |a_j|<=2n^(5/4).
Each C_j is globally sign-symmetric, so even pointwise

    max_(C_j) |a_j+g_j^Tz|=|a_j|+max_(C_j) g_j^Tz.

Consequently the maximum of these actual field-plus-offset widths
over the central cells still has normalized limit sqrt(kappa w).
Any valid field-width upper W_j>=E max_(C_j) g_j^Tz, including an
improved ellipsoid metric, must respect (8).

For the centered sign law, w=1, so the current all-cell linear-field
upper cannot be evaluated below 2sqrt(kappa)n^(3/2)+o(n^(3/2)).
More generally a fixed shifted-sign threshold has w=1-s^2. Along an
actual source subsequence with Phi(A)/n^(3/2)->alpha, evaluating this
field-upper architecture at the target 2sqrt(2)Phi(A) necessarily
requires

                         w<=2alpha^2/kappa.             (9)

This is a necessary condition for that UPPER ARGUMENT, not a new
constraint on actual sources. In particular suppressing the central
cells with bias still leaves the separate competition from active and
near-active original cells; conditional optimality does not evaluate it.

Crucially, (8) is NOT a lower bound on the original Gaussian cross
process X_z=x^TZy. That process is upper-compared to the linear fields,
not identified with them. This result does not preclude a sharper
cross-process comparison, a biased-threshold argument satisfying (9),
or another proof of the original convergence problem. It identifies
why another centered all-cell linear-field metric cannot by itself
finish the current approach. No generic Jtheta metric is packaged here.

## 5. Provenance

The complete 381-line weighted-shell prerequisite was read:
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
It supplies the actual positive M_j and the polynomial partition.
The proof author derived the actual central-set, sign/magnitude,
selection, and covariance-coupling argument; root requested the
classified diagnostic and specified the scope distinctions. A targeted
search found no existing proof of this exact actual-cell saturation.
Only this separate /tmp source was written. No mathematical execution,
canonical change, publication, or backup was performed by its author.
