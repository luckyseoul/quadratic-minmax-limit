# Retaining both initial phases improves the disagreement floor eightfold

2026-09-24. Supporting analytic result, author-reviewed in this session.
The original convergence problem remains OPEN. No independent human or
proof-assistant review is claimed.

The earlier September 24 note selected one orientation per coordinate and
discarded the other one's nonnegative contribution. A new square-completion
identity bounds the input/field residual variance in BOTH orientations.
Convexity then improves the initial mean changed-coordinate fraction from
1/10 to 1/5, and the usable eta floor from 1/1000 to 1/125.

## 1. Statement

Use the initial tilted paired sources with fixed t=993/1000 and

    kappa=2/pi, M=A/sqrt(n), q=(n-1)/n, d=1+t^2 q,
    R_s=(I+s t M)^2/d, G_s ~ N(0,R_s), X_s=sign(G_s),
    F_(s,i)=(s M X_s)_i, s in {-1,+1}.

Keep the old sign at a zero field. Let k_s count changed coordinates in
the synchronous best response and let

    eta_n=(1/(2n^3))sum_s E[k_s^3-1_(k_s=1)].

Uniformly over each fixed cap ||M||op<=L, the new conclusions are

    (1/(2n))sum_s E k_s >= H(t)-o_L(1),   H(t)>1/5,       (1)
    eta_n >= 1/125-o_L(1).                               (2)

Define the same limiting constants used by the tilted theorem:

    z=t^2/(1+t^2), a=1-kappa asin(z)+kappa z,
    e=kappa t/(1+t^2), f=sqrt(kappa a),
    C=2sqrt(kappa)t/(1+t^2), beta=a-kappa>0.

One explicit bound in (1) is

    H(t)=min{1/4, atan(sqrt(beta)/C)/pi}.                 (3)

The constants are independent of the fixed cap. The errors and the
finite-order threshold can depend on it. These are estimates for the
initial source laws, not arbitrary laws after an update.

## 2. Reused hypotheses and the new residual-variance identity

Reuse the tilted covariance expansion, the complete-entry moment
inequality, and the distinguished-input joint Gaussianization from
[the previous source-disagreement note](NOTE_2026-09-24_PAIRED_DISAGREEMENT_FLOOR.md),
Sections 1--3. Put

    z_q=t^2 q/d, b=kappa t^2/d,
    a_q=q[1-kappa asin(z_q)+kappa z_q],
    r_i=(M^3)_ii, w_i=(M^4)_ii.

For v_(s,i)=E F_(s,i)^2 and c_(s,i)=E[G_(s,i)F_(s,i)], those
results give

    v_(s,i)>=a_q+b(w_i-q^2)+s(2kappa t/d)r_i+delta_(s,i),
    (1/n)sum_i |delta_(s,i)|=o_L(1),
    w_i-q^2>=r_i^2/q,
    c_(s,i)=sqrt(kappa)(2tq+s t^2 r_i)/d.                 (4)

After replacing w_i-q^2 by r_i^2/q, the exact new identity is

    a_q+(b/q)r^2+s(2kappa t/d)r
      - kappa(2tq+s t^2 r)^2/d^2
    = a_q-kappa q
      + kappa[t r+s q(1-t^2 q)]^2/(q d^2).               (5)

Thus BOTH phases have the residual-variance floor

    v_(s,i)-c_(s,i)^2 >= beta_q+delta_(s,i),
    beta_q=a_q-kappa q >= (2/3-kappa)q > 0.              (6)

Here a_q/q>=2/3 follows from 0<=z_q<=1/2. In fact for n>=2,
kappa<16/25 and q>=1/2 give beta_q>1/75. The exact remainder in
(5) is a nonnegative square, not a bound obtained by assigning unrelated
worst signs to the two phases.

Because Var(G_(s,i))=1, the quantity v_(s,i)-c_(s,i)^2 is exactly
E[(F_(s,i)-c_(s,i)G_(s,i))^2]. It becomes the conditional variance
in the limiting Gaussian pair. No finite-n Gaussian conditional law of
the actual field is asserted.

## 3. Pair the disagreement angles

Let eps_n tend to zero slowly enough that |delta_(s,i)|<=eps_n in
both phases outside o(n) exceptional coordinates. On the remaining
coordinates put b_n=beta_q-eps_n>0. Then

    v_(s,i)>=c_(s,i)^2+b_n.                              (7)

The existing fixed-cap joint Gaussianization and the positive variance
floor give, uniformly on these coordinates,

    P(Y_(s,i)!=X_(s,i))
       =arccos(c_(s,i)/sqrt(v_(s,i)))/pi+o_L(1).          (8)

Ties at zero local field have probability o_L(1), as in the previous
note. Only one Gaussian input and one field are used in each application.

For b>0, define on c>=0

    phi_b(c)=atan(sqrt(b)/c)/pi, phi_b(0)=1/2.

This decreasing function is convex, since for c>0

    phi_b''(c)=2sqrt(b)c/[pi(b+c^2)^2]>=0.                (9)

The two exact input/field covariances have the fixed mean

    [c_(+,i)+c_(-,i)]/2=C_q=2sqrt(kappa)tq/d>0.          (10)

If both covariances are nonnegative, (7)--(9) and Jensen give

    (1/2)sum_s P(Y_(s,i)!=X_(s,i))
       >= phi_(b_n)(C_q)-o_L(1).                         (11)

If either covariance is negative, its angle in (8) is at least pi/2,
so the phase average is at least 1/4-o_L(1). The other phase contributes
a nonnegative probability. These cases include all coordinates and yield

    (1/2)sum_s P(Y_(s,i)!=X_(s,i))
       >= min{1/4,phi_(b_n)(C_q)}-o_L(1).                (12)

Average over i, discard only the o(n) exceptional coordinates, and take
n to infinity. Now beta_q->beta, C_q->C, and (12) proves (1)--(3).
There is no assumption of independence between changed coordinates or
between phases, and no hybrid phase is constructed.

## 4. A rational certificate for H(t)>1/5

The angle phi_beta(C) has squared cosine C^2/(C^2+beta). Since

    C^2=4kappa t^2/(1+t^2)^2<=kappa,
    beta=a-kappa>0,

monotonicity in C^2 gives

    C^2/(C^2+beta)<=kappa/a
       < (16/25)/[2/3+(7/11)(49/100)]
       =2112/3229 <1309/2000 <cos(pi/5)^2.               (13)

This reuses the exact rational chain from the prior note, including
z>49/100. Hence phi_beta(C)>1/5. The other branch 1/4 in (3) is
also strictly greater than 1/5, proving the strict bound on H(t).

Jensen on the phase/state mixture, with the k=1 correction retained,
now gives

    eta_n >= [(1/(2n))sum_s E k_s]^3-1/n^3
           >=1/125-o_L(1).                              (14)

This improves the conservative initial eta floor by a factor of eight.

## 5. A strict unconditional lower-bound improvement

Keep both the same tilt and the same already-admissible optimizing weight:

    p=(f-2e)/(f+sqrt(f^2+(f-2e)^2)), 0<p<1/2,
    B=B_tilt=e-f/2+sqrt(f^2+(f-2e)^2)/2.

Insert (14) into the finite interior-gap mean-update inequality, using
the existing source inputs e_n>=e-o_L(1), f_n>=f-o_L(1):

    alpha^2-B alpha-p^2/[3000(1+p^2)] >=-o_L(1).

Taking the positive root and removing the fixed operator cap in the
same order of limits as before gives

    liminf alpha_n >= B_both,
    B_both=[B+sqrt(B^2+p^2/[750(1+p^2)])]/2.              (15)

The previous floor 1/1000 used 6000 in place of 750. Therefore

    B_both > B_int > B_tilt                              (16)

without any numerical parameter search. The lower-bound theorem is the
exact expression (15). Using the earlier reported scalar values gives a
rounded estimate B_both approximately 0.3258626; the new remote checker
will evaluate the expression at 60-digit working precision.

For cap removal, apply the fixed-L conclusion at L=K+8 to the
same-order regularized minimizers. It gives
liminf alpha_n>=B_both-2sqrt(2 Gamma/K), with the existing
Gamma=4pi/log(1+sqrt(2)). Let K tend to infinity after n. No error
uniform in a growing cap is assumed. No finite-order cutoff is given.

## 6. Scope, duplication check, and verification

The new contribution is the residual-variance identity (5) followed by
the paired angle estimate (12). The initial one-phase mismatch proof,
source variance estimates, Gaussianization, interior correction, Jensen
step, and cap-removal mechanism are reused explicitly.

The duplicate-work check covered the active checkout, all already-fetched
remote branch heads, the current source and update notes, solution.md,
and the September 6 archive's UNPUBLISHED_WORKING_FRONTIER.md. Only one
local worktree is accessible here. A read of GitHub's main reference
confirmed 889a4c88c8aa0634693e9b7db54de6db4fa7c18d during this session;
the new work continues local commit fe61649aaff6f9882d26e9adf3cc17205474af7f.
No archived family search or finite signing census was restarted.

Analytic review checked both signs in (5), the positivity of beta_q,
the negative-covariance case in (12), the strict rational angle bound,
the phase/state Jensen normalization, and the fixed-cap order of limits.
The focused symbolic/numerical checker is
`both_phase_disagreement_20260924/check.py`. It has not been executed
in this session: the remote SSH preflight failed with a socket permission
error. No remote check or independent proof review is claimed.

The earlier user-reported Orin exit code and its limited scope are recorded
in `paired_disagreement_floor_20260924/USER_REPORTED_CHECK.md`; they are
not evidence that the new two-phase checker passed.

The present variance floor cannot simply be reapplied to non-Gaussian
successor signs. The subsequent [independent flip-noise note](NOTE_2026-09-24_SUCCESSOR_FLIP_NOISE.md)
instead proves an actual-successor estimate by a different moment argument,
and improves (15) further. Long-time trajectory control with a useful
energy conclusion, or an independent cross-order comparison, remains
missing. The original convergence problem remains OPEN.
