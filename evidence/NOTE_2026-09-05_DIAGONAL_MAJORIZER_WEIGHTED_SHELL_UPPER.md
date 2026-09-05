# Weighted linear fields and the original cross-energy ellipsoid metric

2026-09-05. Exact analytic comparison and a polynomial-cell reduction.
The original internal source and cross drift are unchanged. No numerical
evaluation, mathematical computation, signing census, or optimization
was run. The final upper bound is not evaluated at the original target.

The new covariance removes the old scalar tensor normalization, but its
linear field is controlled by WEIGHTED source energies. This note keeps
those energies distinct from the original integer shell values and
identifies the remaining compatibility term explicitly.

## 1. Covariance and the two distinct sets of energies

Let n>=2, N=2n, and

    K=[[A,B],[B^T,-A]],        D=diag(d_1,...,d_N)>0,
    D-K>=0,   D+K>=0,         d_0=tr D<=C_D N^(3/2),

where A is an actual symmetric zero-diagonal signing and B is an n by n
sign matrix. All parameters are deterministic before the Gaussian.
Write

    U=diag(d_1^(-1/2),...,d_n^(-1/2)),
    V=diag(d_(n+1)^(-1/2),...,d_(2n)^(-1/2)),
    A_L=UAU,      A_R=VAV,      W_D=UBV,
    T=D^(-1/2)KD^(-1/2)=[[A_L,W_D],[W_D^T,-A_R]].

Thus ||T||op<=1. The covariance construction and Gaussian reduction are
in `original_mo_diagonal_majorizer_cross_covariance.md`, SHA256
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
The only covariance facts needed here are

    R_0=I+A_L tensor A_R-S_(W_D),       0<=R_0<=2I,
    R_D=R_0+U^2 tensor V^2,            diag R_D=1.       (1.1)

Let k,v>=0, w=k+v<=1. Let Z_0 be a centered Gaussian cross matrix
with covariance kR_0+vI, and put X_(x,y)=x^T Z_0 y. For two states,

    Cov(X_(x,y),X_(x',y'))
      =w(x^T x')(y^T y')
       +k[(x^T A_L x')(y^T A_R y')
                    -(x^T W_D y')(x'^T W_D y)].         (1.2)

For z=(x,y) in the Boolean cube, distinguish the ORIGINAL energies

    p=x^T A x,       q=y^T A y,       c=x^T B y          (1.3)

from the weighted parameters

    p_D=x^T A_L x,   q_D=y^T A_R y,   c_D=x^T W_D y,
    theta(z)=(p_D,q_D,c_D) in [-n,n]^3.                 (1.4)

The last inclusion follows by compressing the contraction T. In
general neither triple determines the other.

For later use, Gaussian variance padding costs at most

    2 sqrt(k n log(2) (tr U^2)(tr V^2))=O_(C_D)(n)      (1.5)

in an expected Boolean maximum, including with arbitrary fixed offsets.
Indeed the independent padding entries have variances
k/(d_i d_(n+j)); every Boolean pairing has variance
k(tr U^2)(tr V^2), and there are 2^(2n) such pairings.
The identity tr D^(-1)<=d_0/(N-1), proved from diag(T^2)<=1,
bounds their product by d_0^2/[4(N-1)^2]. This padding is removed
only from the Gaussian model, never from the actual unit-variance
rounded-sign law.

## 2. Exact positive linear field on a weighted shell

Suppose theta=(p_D,q_D,c_D) is attained by an actual Boolean state.
Define

    M_theta=wn I_(2n)
       +k[[q_D A_L,-c_D W_D],[-c_D W_D^T,p_D A_R]].      (2.1)

Then

    0<=M_theta<=2wn I_(2n).                             (2.2)

In particular no rank correction or Gaussian with an indefinite
covariance is required. To prove this, compress T to the orthonormal
vectors (x,0)/sqrt(n),(0,y)/sqrt(n) at a state attaining theta. The
matrix [[p_D,c_D],[c_D,-q_D]] has operator norm at most n. Its signed
swap conjugate

    H_theta=[[q_D,-c_D],[-c_D,-p_D]]

also has norm at most n. The nonconstant block in (2.1) is the
compression of H_theta tensor T onto the first internal block of the
first tensor copy and the second internal block of the second copy.
Its norm is therefore at most n. Adding wnI proves (2.2), since w>=k.

Let g_theta=(xi,eta) be centered Gaussian with covariance M_theta.
On any nonempty set of states with this EXACT weighted triple,

    E max X_(x,y) <= E max (xi^T x+eta^T y).             (2.3)

Here is the full increment identity. For two states in that set put

    r_x=x^T x',  r_y=y^T y',
    d=x^T W_D y',  e=x'^T W_D y,
    delta x=x-x', delta y=y-y'.

The excess of the linear-field increment variance over that of X,
divided by two, is

    k/4 [<delta x tensor delta y,
                  R_0(delta x tensor delta y)>+(d-e)^2]
             +v(n-r_x)(n-r_y) >=0.                     (2.4)

This follows by expanding (1.2) and (2.1), retaining the exchange
term -de. Positive semidefiniteness in (1.1) then proves the sign.
The finite Gaussian increment comparison gives (2.3), including
singular covariances. The same comparison is valid with identical
arbitrary deterministic offsets attached to the states.

If the full Gaussian covariance kR_D+vI is retained instead, add

    k diag((tr V^2)U^2,(tr U^2)V^2)

to M_theta. It remains positive semidefinite, and (2.4) holds with
R_D in place of R_0. The simplified base form (2.1) suffices below
because of the explicit padding bound (1.5).

## 3. Real weighted values: a controlled polynomial-cell reduction

Exact weighted triples need not have polynomially many distinct values.
Fix 0<delta<=4n and partition [-n,n]^3 into half-open cubes of side
delta, assigning the final boundary to its last cell. In each nonempty
cell of ACTUAL states choose a deterministic representative z_0 and
write theta_0=theta(z_0). Its field covariance M_(theta_0) is genuinely
positive by Section 2; no possibly unattained cell center is substituted.

For any two states z,z' in that cell let

    bar theta=(theta(z)+theta(z'))/2,
    Delta p=p_D(z)-p_D(z'),
    Delta q=q_D(z)-q_D(z'),
    Delta c=c_D(z)-c_D(z').

The same expansion as (2.4), using M_(bar theta) algebraically, gives
the increment excess divided by two exactly as

    k/4 [<delta x tensor delta y,
                  R_0(delta x tensor delta y)>+(d-e)^2
                            -Delta p Delta q+(Delta c)^2]
             +v(n-r_x)(n-r_y).                          (3.1)

No Gaussian with M_(bar theta) is needed for this algebraic identity.
Since every parameter difference is at most delta,
the expression in (3.1) is at least -k delta^2/4.
Moreover

    ||M_(theta_0)-M_(bar theta)||op<=2k delta,
    ||z-z'||^2<=8n.

Replacing the midpoint covariance by the representative therefore
leaves increment excess divided by two at least

                       -8kn delta-k delta^2/4.          (3.2)

Introduce independent centered Gaussians epsilon_z, one for each
state in the cell, independent also of g_(theta_0), with common
variance

                       sigma_delta^2=8kn delta+k delta^2/4.

For distinct states their increment variance is 2sigma_delta^2, so
(3.2) proves valid increment domination by
`g_(theta_0)^T z+epsilon_z`. The equal-state assertion is trivial.
It follows that, with arbitrary matching deterministic offsets,

    E max_cell (a_z+X_z)
      <= E max_cell (a_z+g_(theta_0)^T z)
                           +6n sqrt(k delta log(2)).     (3.3)

Indeed the cell has at most 2^(2n) states and
sigma_delta^2<=9kn delta; the Gaussian finite-maximum bound gives
the last term. This is an exact comparison error, not a claim that
all weighted energies in the cell are equal.

Take delta=1/n. There are at most (2n^2+1)^3 weighted cells. Further
partition by the original integer triple (p,q,c), which lies in
[-n^2,n^2]^3. The total number m of nonempty cells is therefore

                         m<=(2n^2+1)^6.                (3.4)

The comparison error in (3.3) is now at most
`6 sqrt(k n log(2))`, uniformly over all cells.

Selecting cells also has a controlled cost. A supremum of the base
process X over any fixed subset is a Lipschitz function of its
underlying standard Gaussian with Lipschitz constant at most
`n sqrt(2k+v)`: each coefficient state has norm n and the base
covariance has norm at most 2k+v. Gaussian concentration and the
finite exponential bound consequently give

    E max_(j<=m) [a_j+max_(z in cell j) X_z]
      <= max_(j<=m) [a_j+E max_(z in cell j) X_z]
                          +n sqrt(2(2k+v)log(m)).        (3.5)

Arbitrary deterministic cell offsets do not change this concentration
argument. For a maximum of absolute values include both signs of
the objective and replace log(m) by log(2m). No independence between
the original cell suprema is assumed. Thus the complete cell-selection
error is O(n sqrt(log n))=o(n^(3/2)); the exact weighted-value count
has never been asserted to be polynomial.

## 4. The same D majorizes the ORIGINAL cross-energy form

Put

    H_B=[[0,B],[B^T,0]],      J=diag(I_n,-I_n).

Conjugating D+-K by J preserves D. Since
`H_B=(K-JKJ)/2`, averaging `D-K` with `D+JKJ`, and then reversing
the signs, proves

                             D-H_B>=0, D+H_B>=0.        (4.1)

Define the contraction and the metric

    L_D=D^(-1/2)H_BD^(-1/2)=[[0,W_D],[W_D^T,0]],
    P_eta=D-eta H_B,     E_eta=(1-|eta|)D,   |eta|<1.   (4.2)

Then P_eta>0 and 0<=E_eta<=P_eta. On any cell with ORIGINAL cross
energy c, the metric radius is exactly

    z^T P_eta z=d_0-2eta c=d_0(1-eta u),
    u=2c/d_0 in [-1,1].                                 (4.3)

It is c, not c_D, that appears in this shell radius. Thus this
metric uses precisely the diagonal majorizer that supplied the
weighted covariance while preserving the original cross constraint.

For a representative theta define

    cal H_theta=[[q_D A,-c_D B],[-c_D B^T,p_D A]],
    Mhat_theta=D^(-1/2)M_theta D^(-1/2)
              =wn D^(-1)+kD^(-1)cal H_theta D^(-1),
    F_eta=(I-eta L_D)^(-1),
    T_eta=tr(Mhat_theta F_eta),
    R_eta=tr(Mhat_theta F_eta^2).                         (4.4)

These are the exact two field traces:
`T_eta=tr(M_theta P_eta^(-1))` and
`R_eta=tr(D P_eta^(-1)M_theta P_eta^(-1))`.
In particular no simultaneous diagonalization of Mhat_theta and
L_D is assumed. The cyclic trace identities suffice.

The Boolean ellipsoid remainder theorem of
`NOTE_2026-09-05_BOOLEAN_ELLIPSOID_SHELL_UPPER.md`, Sections 1 and 4,
now gives the genuine upper
bound on the representative field over the actual cell,

    W(theta,c):=inf_(|eta|<1) sqrt(d_0) {
       sqrt((|eta|-eta u)[T_eta-(1-|eta|)R_eta])
                    +sqrt(kappa)(1-|eta|)sqrt(R_eta)},   (4.5)

where kappa=2/pi. In other words
`E max_cell g_theta^T z<=W(theta,c)`.
The expression `T_eta-(1-|eta|)R_eta` is nonnegative: it is
`tr(M_theta P_eta^(-1)(P_eta-E_eta)P_eta^(-1))`.
The cited exact Boolean remainder theorem applies because (4.3)
holds for every state in the cell. Refining a shell does not affect
that hypothesis. Each eta is fixed in this finite-dimensional
inequality; no uncontrolled endpoint exchange of limits is used.

For clarity, its short completion-square mechanism is as follows.
For a Gaussian g with covariance M, shell radius z^TPz=q, and
0<=E<=P with E diagonal, completing the square and retaining
`sum_i E_ii(|(P^(-1)g)_i|/t-1)^2` gives, after optimizing t>0,

    E max g^Tz <= sqrt((q-tr E)
       [tr(MP^(-1))-tr(E P^(-1)MP^(-1))])
           +sqrt(kappa) sum_i E_ii sqrt((P^(-1)MP^(-1))_ii).

Cauchy--Schwarz on the final sum with E=(1-|eta|)D yields (4.5).
No independence of field coordinates is used.

Combining (3.3)--(3.5) with (4.5) gives, for any deterministic drift
constant a_j on each original triple, the fully justified upper

    E max_z |a_z+X_z|
      <= max_(nonempty actual cells j) [|a_j|+W(theta_j,c_j)]
         +6 sqrt(k n log(2))
         +n sqrt(2(2k+v)log(2m)).                        (4.6)

For the original block objective the drift is exactly
`a_j=(p_j-q_j)/2+s c_j`; neither A nor B has been reweighted in it.
For the full Gaussian covariance kR_D+vI, add (1.5) to (4.6).
The tuples in the maximum are ACTUAL nonempty cells and their actual
representatives, not an independently relaxed box of energy values.

## 5. Exact first feedback and the explicit compatibility remainder

If w=0 all fields vanish. Otherwise set

    r_L=sum_(i<=n) d_i^(-3/2),
    r_R=sum_(j>n) d_j^(-3/2),
    T_0=wn tr D^(-1)>0.

The zero diagonals of A and H_B and the sign identity B_ij^2=1 give

    tr(Mhat_theta L_D)=-2k c_D r_L r_R,
    v_0:=tr(Mhat_theta L_D)/T_0
                            =-2k c_D r_L r_R/T_0.       (5.1)

Consequently the simpler affine ellipsoid consequence is also valid:

    E max_cell g_theta^T z <= sqrt(d_0 T_0)
       inf_(|eta|<1) sqrt(
        [1-eta u-(1-kappa)(1-|eta|)](1+eta v_0)
                                      /(1-eta^2)).       (5.2)

One has |v_0|<=1 by positivity of Mhat_theta and contraction of L_D.
This affine consequence uses the resolvent upper bound, not an
equality replacing the stronger two-trace expression in (4.5).

To state the mismatch without hiding it in notation, choose any fixed
positive reference scales alpha,beta, and define at an actual state

    Delta_L=x^T(UAU-alpha^2 A)x,
    Delta_R=y^T(VAV-beta^2 A)y,
    Delta_B=x^T(UBV-alpha beta B)y.                       (5.3)

Then, exactly,

    p_D=alpha^2 p+Delta_L,
    q_D=beta^2 q+Delta_R,
    c_D=alpha beta c+Delta_B.                            (5.4)

For example alpha=tr(U)/n and beta=tr(V)/n are deterministic choices;
no particular choice is needed for the identities. Equation (5.1)
becomes

    v_0=-lambda_D u-(2k r_L r_R/T_0)Delta_B,
    lambda_D=k alpha beta d_0 r_L r_R/T_0.               (5.5)

This is the explicit original-cross/weighted-feedback remainder.
Replacing v_0 by a multiple of -u requires a justified estimate on
Delta_B. The original source energy constraints do not identify it.

The stronger trace bound also retains internal compatibility terms.
Let E=U^2, F=V^2 and write the block diagonal resolvents as

    F_(eta,LL)=(I-eta^2 W_D W_D^T)^(-1),
    F_(eta,RR)=(I-eta^2 W_D^T W_D)^(-1).

Their exact contribution to T_eta is

    k q_D tr(E A E F_(eta,LL))
                 +k p_D tr(F A F F_(eta,RR)).            (5.6)

These vanish at eta=0 but have no prescribed sign in general.
The analogous contribution to R_eta uses the respective diagonal
blocks of F_eta^2. In particular p=q=0 does not make (5.6) vanish:
the weighted parameters then equal Delta_L,Delta_R.

Equations (5.3)--(5.6) are not asserted to be small under the trace
cap alone. Nor is an impossibility theorem for controlling them
claimed. They specify the additional estimate needed to turn this
same-source weighted upper into an evaluated original-order result.

## 6. Status and exact remaining implication

The local-scale Gaussian covariance, positive weighted reference
fields, polynomial-cell reduction, and same-D original-cross metric
are valid. They remove neither the unweighted internal drift nor the
independent Hermite cushion. Their auxiliary Gaussian comparison
errors are o(n^(3/2)) under the stated fixed trace cap.

What remains is an upper evaluation of (4.6) on the ACTUAL coupled
original/weighted cells, using any available conditional optimality
and source constraints. One may use either the exact traces (4.5)
or the weaker affine consequence (5.2), but must retain or bound the
compatibility remainders just displayed. No scalar weak-feedback
counterexample is promoted to an actual source, and no conclusion
about original MO convergence is claimed.
