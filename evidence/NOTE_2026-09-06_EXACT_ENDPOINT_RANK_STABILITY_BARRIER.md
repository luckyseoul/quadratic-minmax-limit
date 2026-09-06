# Exact endpoint rank obstruction is not asymptotically stable

2026-09-06. Bounded analytic barrier, not a milestone or an actual
counterexample to the source-cap problem. No computation was run.

## 1. Exact obstruction

For a complete symmetric zero-diagonal sign matrix A of order n,
reduction modulo two gives A=J-I over F_2. If n is even this matrix
is invertible; if n is odd its kernel is the span of the all-one
vector. A nonzero minor modulo two is nonzero over the reals, so

    n even: ker(A)=0;
    n odd:  dim ker(A)<=1.                                    (1.1)

Suppose a scalar paired contraction has H=A/d and W=B/d, where
all nonzero singular values of B equal d. Put
P_L=WW^T and P_R=W^TW; both are projections. The two block row
inequalities give H^2+P_L<=I and H^2+P_R<=I. Hence H annihilates
both cross ranges, and

                             AB=BA=0.

Thus rank(B)<=1. Even n is impossible for nonzero complete B;
for odd n a complete rank-one B has its single singular value n,
so d=n. In particular an exactly flat endpoint cross block of
positive limiting rank fraction is impossible at scale d=O(sqrt(n)).

This is only an exact rank obstruction. It supplies no positive
lower bound on the nonzero singular values of A after division
by sqrt(n). The next family shows why that distinction matters.

## 2. Explicit actual near-annihilation family

Let ell>=1, t=4^ell, n=2t, and define

    F=2I_4-J_4,       H=F^{tensor ell},
    C_2=[[1,-1],[-1,1]],
    A=H tensor J_2-I_n,       B=H tensor C_2.

Directly F^2=4I_4 and diag(F)=1. Thus H is symmetric,
H^2=tI_t, and diag(H)=1. Consequently A has zero diagonal and
every off-diagonal entry is a sign, while every entry of B is a sign.
Since J_2 C_2=C_2 J_2=0,

                             AB=BA=-B.                         (2.1)

On the within-pair constant subspace, A has eigenvalues
2lambda(H)-1=+/-2sqrt(t)-1 and B vanishes. On the within-pair
zero-sum subspace, A=-I and B has eigenvalues +/-2sqrt(t).
The latter subspace has dimension t=n/2. Thus B has rank n/2
and every nonzero singular value is exactly 2sqrt(t).

For the actual paired matrix K=[[A,B],[B,-A]], commutation in
(2.1) gives K^2=diag(A^2+B^2,A^2+B^2). Therefore

    ||K||op=d:=2sqrt(t)+1,
    D=d I_(2n) is feasible,       delta=0.

Both signs of sqrt(t) occur in H, so the claimed maximum is attained.
Writing N=2n=4t, the trace cap is explicitly

    tr(D)/N^(3/2)=1+1/(2sqrt(t))<=5/4.                        (2.2)

The FULL actual W=B/d has squared-singular-value law

    nu_t=(1/2)delta_0+(1/2)delta_(4t/d^2)
                       ->(1/2)delta_0+(1/2)delta_1.             (2.3)

Nevertheless A remains invertible, with the n/2-dimensional
eigenspace A=-I supporting the entire range of B. Its normalized
annihilation error tends to zero even in operator norm:

    ||(A/d)(B/d)||op=2sqrt(t)/d^2 ->0.                         (2.4)

Thus exact scalar feasibility, a fixed trace cap, full endpoint-law
convergence, and approximate annihilation do NOT upgrade (1.1)
to an asymptotic rank exclusion. Complete B can indeed occupy a
macroscopic small-eigenvalue subspace of a complete A.

## 3. Why this does not realize the strengthened formal profile

This family fails the crucial small original-source norm cap.
The vector 1_4 is an eigenvector of F with eigenvalue -2, and
(1,1,-1,-1) is a sign eigenvector with eigenvalue +2. Their tensor
products give a Boolean vector v with H v=-sqrt(t)v. Hence
v tensor 1_2 is a Boolean A-eigenvector with eigenvalue -d,
so the spectral upper bound on its original quadratic norm is attained:

    Phi(A)=n d/2,
    Phi(A)/n^(3/2)->1/sqrt(2), not 2/5.                       (3.1)

Similarly a Boolean tensor eigenvector attains beta(B)=2n sqrt(t).
The paired norm satisfies beta(B)<=Phi(K)<=n d, and therefore

                    Phi(K)/n^(3/2)->sqrt(2), not 4/3.          (3.2)

The family is a countermechanism to the GENERIC rank-stability
argument only. It does not settle a rigidity statement that also
uses the original source cap Phi(A)<= (2/5+o(1))n^(3/2), and it
does not realize the specific m=9/25 formal profile.

The proof worker derived the exact parity obstruction and this
explicit family after the root requested a stability check and
suggested the replicated-row mechanism. Root checked the displayed
algebra. No canonical source was edited; this note records a bounded
barrier so the exact-rank argument is not mistaken for asymptotic
progress on the original small-source-norm branch.
