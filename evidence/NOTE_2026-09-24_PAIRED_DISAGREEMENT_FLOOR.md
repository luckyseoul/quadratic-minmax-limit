# A positive disagreement floor for the initial paired Gaussian phases

2026-09-24. Supporting analytic result, reviewed by its author in this
session. The original convergence problem remains OPEN. There is no
independent human review or proof-assistant verification of this note.

The new implication is a cap-independent positive lower bound on the
changed-coordinate statistic in the September 23 interior-gap inequality,
for the INITIAL tilted paired Gaussian laws. It makes that correction
strictly positive and gives a strict improvement of the exact tilted lower
bound. It makes no Gaussian assumption about a post-update law.

## 1. Statement and reused inputs

Fix t=993/1000 and kappa=2/pi. For a complete signing A of order n put

    M=A/sqrt(n), q=(n-1)/n, alpha=Phi(A)/n^(3/2),
    d=1+t^2 q, R_s=(I+s t M)^2/d,
    G_s ~ N(0,R_s), X_s=sign(G_s), F_(s,i)=(s M X_s)_i.

Let Y_s be the synchronous best response, retaining X_(s,i) at a zero
local field. Define

    k_s=#{i:Y_(s,i)!=X_(s,i)},
    g(k)=k^3-1_(k=1),
    eta_n=(1/(2n^3)) sum_s E g(k_s).

For every fixed operator cap ||M||op<=L, uniformly over this class,

    liminf (1/(2n)) sum_s E k_s >= h(t)>1/10,             (1)
    liminf eta_n >= h(t)^3 > 1/1000,                      (2)

where, using the existing tilted-source constants,

    z=t^2/(1+t^2), a=1-kappa asin(z)+kappa z,
    e=kappa t/(1+t^2), f=sqrt(kappa a),
    rho=2e/f, h(t)=arccos(rho)/(2pi).

The lower constants do not depend on L; the vanishing errors may depend
on L. The order of limits is still n first, L second.

Reused inputs, with their original scopes:

* [Tilted paired fields](NOTE_2026-09-23_TILTED_PAIRED_FIELD_LOWER.md),
  Sections 1--3: actual phase covariance and variance estimates.
* [Distinguished-coordinate Gaussianization](NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md),
  Section 6: joint convergence of ONE Gaussian input and ONE local field
  under a fixed covariance cap. The lemma itself does not require the
  near-flat hypotheses of that note's main theorem.
* [Finite interior gap](NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md),
  equation (9): the strengthened mean update for arbitrary actual laws.
* [Same-order regularization](NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md),
  Section 5: removal of a fixed source operator cap at the minimum.

## 2. Joint Gaussianization with either covariance sign

Section 6 of the cited Gaussianization proof establishes convergence of
(G_i,F) to a Gaussian pair with its actual covariance, uniformly by a
compactness argument. The stated sign-probability corollary assumes
E[G_i F]>=c_*>0 only to ensure a nonzero limiting field variance.

The identical proof gives the following slightly more convenient form:
if Var(F)>=v_*>0 and the same fixed covariance-cap and delocalized-row
hypotheses hold, then, uniformly, with c=E[G_i F] and v=E F^2,

    P(sign(F)!=sign(G_i))=arccos(c/sqrt(v))/pi+o(1).       (3)

Indeed all convergent covariance subsequences have Var(G_i)=1 and
Var(F)>=v_*. The boundary of the sign event is the union of the two
coordinate axes, of zero limiting Gaussian probability. This works for
negative c and for degenerate pairs with correlation +/-1 too. It also
gives P(F=0)=o(1), so retaining the old sign at zero field has no effect
on (3). No new growing-dimensional central limit theorem is invoked.

## 3. Select one orientation per coordinate

Write

    r_i=(M^3)_ii, w_i=(M^4)_ii,
    b=kappa t^2/d, c0=1-kappa asin(t^2 q/d),
    a_q=c0 q+b q^2, H=2kappa t/d.

The tilted covariance estimate gives

    v_(s,i):=E F_(s,i)^2
       >= a_q+b(w_i-q^2)+s H r_i+delta_(s,i),             (4)
    (1/n)sum_i |delta_(s,i)|=o_L(1).

The complete-entry moment inequality is r_i^2<=q(w_i-q^2).
Also the Gaussian sign identity E[G_i sign(G_j)]=sqrt(kappa) R_ij
gives the EXACT input/field covariance

    c_(s,i):=E[G_(s,i) F_(s,i)]
       =sqrt(kappa)(2tq+s t^2 r_i)/d.                     (5)

Choose s_i so that s_i r_i=-|r_i|, with either choice at r_i=0,
and put u=|r_i|. For this chosen orientation, (4)--(5) imply

    v_(s_i,i) >= U(u)+delta_(s_i,i),
    U(u)=a_q+(b/q)u^2-Hu,
    c_(s_i,i)=C_q(1-tu/(2q)), C_q=2sqrt(kappa)tq/d.       (6)

For n>=20, the fixed tilt has d>48/25. Since 7/11<kappa<16/25
and 0<=t^2 q/d<=1/2, elementary arcsine bounds give

    2/3 <= a_q/q=1-kappa asin(t^2 q/d)+kappa t^2 q/d <=1,
    2kappa/d < 2/3,      4kappa/d > 1.                    (7)

The relevant exact scalar identity is

    U(u)-a_q(1-tu/(2q))^2
      = (t/q)(a_q-2kappa q/d)u
        + (t^2/(4q^2))(4kappa q/d-a_q)u^2 >=0.            (8)

Both coefficients are positive by (7). In addition, completing the
square gives the uniform positive floor

    U(u)>=a_q-kappa q/d>q/3>=19/60.                       (9)

Discard the o(n) coordinates on which either covariance-expansion error
has magnitude larger than a vanishing tolerance. On all remaining
coordinates (3) applies uniformly; their variance is bounded below by,
for example, 19/120 for sufficiently large n.

If c_(s_i,i)>=0, equations (6)--(9) yield

    c_(s_i,i)/sqrt(v_(s_i,i)) <= C_q/sqrt(a_q)+o_L(1).

If c_(s_i,i)<0, (3) instead gives disagreement probability at least
1/2-o_L(1). Thus in either case

    P(Y_(s_i,i)!=X_(s_i,i))
       >= arccos(C_q/sqrt(a_q))/pi-o_L(1).                (10)

Here C_q^2<=kappa q and a_q>=2q/3, so the arccos argument is uniformly
less than one. As n tends to infinity it tends to rho=2e/f.

Summing (10), retaining the nonnegative contribution of the OTHER phase,
and dividing by 2n proves (1) with h(t)=arccos(rho)/(2pi). The selected
orientation can differ by coordinate: this is only a lower bound on the
sum of the two actual marginal probabilities, not a new hybrid phase.
No independence of changed coordinates is assumed.

## 4. An exact conservative floor for eta

For our rational t, z>49/100. Since 4t^2/(1+t^2)^2<=1,

    rho^2 = 4kappa t^2/[(1+t^2)^2 a]
       < (16/25)/[2/3+(7/11)(49/100)]
       = 2112/3229
       < 1309/2000
       < (3+sqrt(5))/8 = cos(pi/5)^2.                    (11)

The middle rational comparison is 4224000<4226761. The last uses
sqrt(5)>559/250, whose squared comparison is 312481<312500.
Consequently rho<cos(pi/5), hence h(t)>1/10.

Apply Jensen to k/n on the equally weighted mixture of the two phase
laws. Because g(k)=k^3-1_(k=1),

    eta_n >= [(1/(2n))sum_s E k_s]^3-1/n^3.               (12)

Equations (1), (11), and (12) prove (2). In particular the usable
cap-independent bound is eta_n>=1/1000-o_L(1).

## 5. Strict lower-bound consequence without another optimization

Keep the PREVIOUS optimizing weight, rather than retune any parameter:

    p=(f-2e)/(f+sqrt(f^2+(f-2e)^2)),
    B=e-f/2+sqrt(f^2+(f-2e)^2)/2 = B_tilt.

Equation (11) implies f>2e>0, so 0<p<1/2. Thus the finite interior-gap
mean update applies at this fixed weight:

    (1+p^2)alpha
       >= (1-p)^2 e_n+p(1-p)f_n+p^2 eta_n/(24alpha).

For each fixed L, use the already-proved bounds e_n>=e-o_L(1),
f_n>=f-o_L(1), and the new eta_n>=1/1000-o_L(1). Multiplying by alpha,
which is positive and bounded at fixed L, gives

    alpha^2-B alpha-p^2/[24000(1+p^2)] >= -o_L(1).

The positive root is

    B_int = [B+sqrt(B^2+p^2/[6000(1+p^2)])]/2 > B.        (13)

Hence liminf alpha_n^(L)>=B_int for every fixed eligible cap. Applying
same-order regularization at L=K+8 gives

    liminf alpha_n >= B_int-2sqrt(2 Gamma/K),
    Gamma=4pi/log(1+sqrt(2)).

Let K tend to infinity AFTER n. The resulting unconditional conclusion is

    liminf alpha_n >= B_int > B_tilt.                     (14)

This is an exact expression, not an unverified new decimal evaluation.
The stronger h(t)^3 in (2) could be retained in (13), but is unnecessary
for a strict advance. There is no new finite-order cutoff.

## 6. Remaining gap, provenance, and verification limits

This controls eta for the initial tilted Gaussian phases. It does NOT
control eta, field covariances, or the energy gain for arbitrarily many
actual updates. In particular (3) cannot be reapplied to sign(M X_s)
merely because it was valid for X_s=sign(G_s). Nor does this prove a
one-vertex row controlling the stable-state deficit shells.

The same-order bounds now have a quantitatively positive interior term;
the original convergence problem remains OPEN. The next substantial
implication would be a quantitative update estimate for the actual
non-Gaussian successor law, or the independent neutral-increment
extension estimate with summable normalized excess.

Baseline: main at 889a4c88c8aa0634693e9b7db54de6db4fa7c18d.
The duplicate-work check covered the local checkout, all fetched remote
branch heads, the active source/update notes, and the analytic scratch
member UNPUBLISHED_WORKING_FRONTIER.md in the September 6 campaign archive.
The archived cubic estimate and near-flat mismatch estimate are reused
context, not claimed new. This environment exposes only the new checkout;
unpublished changes on the user's other machines were not inspected.

Verification here is analytic author review of (5), (8), (9), (11), (12),
the positive-root algebra, and the reuse of the joint-Gaussianization proof.
No matrix census, simulation, numerical constant evaluation, or test run
was performed. Repository instructions require offloaded mathematical
checks; the NUKA preflight failed because its hostname was not resolvable
from this environment. No mesh receipt or independent review is claimed.
This is a supporting lemma and constant correction, not a major structural
milestone; no new large-drive backup was taken.
