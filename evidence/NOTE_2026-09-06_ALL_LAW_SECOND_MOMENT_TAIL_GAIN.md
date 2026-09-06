# All-law original-source gain with second-moment spectral-tail control

2026-09-06. Analytic extension for ACTUAL complete symmetric signings.
No mathematical program, checker, solver, spectral scan, numerical
integral, optimization run, or construction search was executed.

## 1. The fixed-cutoff asymptotic theorem

Let A be an actual complete real symmetric zero-diagonal signing of
order n, with n tending to infinity, and write

    M=A/sqrt(n),       q=(n-1)/n,
    Q_M(x)=x^T M x/2,
    Phi(A)=max_(x in {+1,-1}^n)|x^T A x/2|,
    alpha(A)=Phi(A)/n^(3/2)=max_(x in {+1,-1}^n)|Q_M(x)|/n,
    ell=tr|M|/n,       kappa=2/pi.

Fix C>=1. Assume only the normalized SECOND-MOMENT tail condition

    t_n=(1/n)tr[M^2 1_(|M|>C)] ->0.                       (1.1)

For a,w>=0 let

    Psi(a,w)=[(a+w)arctan(sqrt(w/a))-sqrt(a w)]/pi,
                          when a>0,
    Psi(0,w)=w/2,       Psi(a,0)=0,
    D_C(ell)=(1-exp(-1))/(2C)
                         Psi(kappa C^2,(1-kappa)ell^2),
    F_C(ell)=kappa/(2ell)+D_C(ell).

Then the EXACT same gain function as in the bounded all-law theorem
satisfies the asymptotic conclusion

    liminf [alpha(A)-F_C(ell)]>=0.                         (1.2)

The nuclear moment ell need not converge. No limiting spectral law,
diagonal homogeneity, flat spectrum, or actual operator cap is assumed.
Without a tail rate in (1.1), no error depending only on C and n is
promised uniformly over all sources. Every energy below is the full
original same-source energy; no source coordinates are removed.

## 2. Bounded bulk and one common good-coordinate set

Put

    M_b=M 1_(|M|<=C),       M_o=M-M_b,       L_b=|M_b|.

Spectral orthogonality and completeness give

    ||M_b||op<=C,       M_b M_o=0,
    M^2=M_b^2+M_o^2,       (M^2)_ii=q,
    tr M_o^2=nt_n=o(n),       ||M_o||_1<=n sqrt(t_n)=o(n),
    ell_b=tr L_b/n=ell-o(1),
    (q-t_n)/C<=ell_b<=ell<=sqrt(q).                        (2.1)

For matrices ||.||_1 denotes the nuclear norm; for vectors it denotes
the usual sum of absolute values. The nuclear difference in (2.1)
uses |M|=L_b+|M_o|. The lower bound follows from M_b^2<=C L_b.
In particular ell stays bounded away from zero.

Define

    u_i=(M_o^2)_ii,       h_i=(L_b)_ii,
    delta_i=(M_b)_ii=-(M_o)_ii,
    eta_n=max(n^(-1/8),t_n^(1/4)),
    I={i:u_i<=eta_n^2},       m=|I|.

Then eta_n->0, (n-m)/n<=t_n/eta_n^2->0, and on I,

    |delta_i|<=eta_n,
    (q-eta_n^2)/C<=h_i<=sqrt(q).                           (2.2)

The lower bound again uses M_b^2<=C L_b; the upper bound is
coordinatewise spectral Cauchy--Schwarz. This same I is used for both
phases. All coordinates outside I remain present as genuine variables.

## 3. Each truncated phase uses its own actual diagonal

For s in {+1,-1} set

    P_s=L_b+sM_b>=0,       ||P_s||op<=2C,
    d_(s,i)=(P_s)_ii=h_i+s delta_i.

On I, d_(s,i)>0 for all sufficiently large n. Define R_s[I,I] by
normalizing P_s[I,I] with these OWN diagonals. Set R_s[I^c,I^c]=I
and the two off-blocks to zero. Thus

    diag R_s=1,       R_s>=0,
    ||R_s||op<=B_n=2C^2/(q-eta_n^2-C eta_n)
                         =2C^2/q+o(1).                   (3.1)

The displayed denominator is positive eventually and B_n>=1. The
normalizations are not assumed equal. Let J be diagonal with
J_ii=h_i^(-1/2) on I and zero outside I. On I the diagonal ratio
sqrt(h_i/d_(s,i)) is 1+O_C(eta_n), uniformly in i and s. Consequently

    ||R_s-(J P_s J+I_(I^c))||op=O_C(eta_n)=o(1),
    R_+-R_-=2J M_b J+E_n,       ||E_n||op=o(1).            (3.2)

This follows by conjugating J P_s J by that diagonal ratio; the
middle matrix has bounded operator norm by (2.2). The independent
identity blocks cancel in the second formula.

Let G_s be centered Gaussian with covariance R_s and X_s=sign(G_s),
using sign(0)=+1. Each coordinate of G_s is standard Gaussian, and
singular joint covariances are allowed. Schur multiplication by a
correlation matrix preserves PSD order and fixes I, so for every k>=1,

    0<=R_s^(circ k)<=B_n I.

The normalized Hermite expansion of sign gives

    C_(X_s)=E[X_s X_s^T]=kappa R_s+C_(tail,s)<=B_n I,
    C_(tail,s)=sum_(odd k>=3)c_k^2 R_s^(circ k),
    sum_(odd k>=3)c_k^2=1-kappa.                           (3.3)

The series converges in operator norm. Actual complete entries imply,
uniformly in odd k>=3 and in s,

    |tr(M R_s^(circ k))|
       <=n^(-1/2)sum_(i!=j)|R_(s,ij)|^k
       <=n^(-1/2)tr R_s^2<=B_n sqrt(n).                  (3.4)

## 4. The same original nuclear baseline survives for these phases

Write e_s=E Q_(sM)(X_s)/n for the two oriented actual phase energies.
Equations (3.3)--(3.4) show

    e_s=(s kappa/(2n))tr(M R_s)+o(1).

Use (3.2) and ||M||_1<=sqrt(n)||M||_F<=n to obtain

    (e_++e_-)/2=(kappa/(2n))tr(M J M_b J)+o(1).

Replacing M_b by M in this trace costs at most

    ||J||op^2 ||M||_F ||M_o||_F=o(n),

because J is uniformly bounded. The complete off-diagonal squares are
M_ij^2=1/n. Therefore

    (e_++e_-)/2
       =(kappa/(2n^2))sum_(i,j in I, i!=j)1/sqrt(h_i h_j)+o(1).

The elementary inequalities 1/sqrt(ab)>=2/(a+b) and Jensen on the
unordered pairs give

    sum_(i,j in I, i!=j)1/sqrt(h_i h_j)
                              >=m^2(m-1)/(sum_(i in I)h_i).

Since sum_(i in I)h_i<=n ell_b<=n ell and m/n->1, we conclude

    (e_++e_-)/2>=kappa/(2ell)-o(1).                        (4.1)

This is a lower on the average of the SAME two actual phase energies
which will be improved below. It is not merely an unrelated already
maximized nuclear inequality. No global Lipschitz bound for arcsine,
common actual phase diagonal, or actual operator cap was used.

## 5. Full actual higher-chaos noise, uniformly over the Hermite tail

For a fixed phase put K_s=diag(sqrt(d_(s,i))) on ALL coordinates.
Some of its bad-coordinate entries may be zero; it is never inverted.
On I the matrix K_s R_s K_s agrees with P_s, and on I^c it is the
diagonal of P_s. Hence

    Z_s=K_s R_s K_s-P_s,
    rank Z_s<=2|I^c|,       ||Z_s||op<=4C,
    ||Z_s||_1=o(n),
    L_b=K_s R_s K_s-sM_b-Z_s.                             (5.1)

For odd k>=3, pairing the last identity with R_s^(circ k) gives

    tr(L_b R_s^(circ k))
       =sum_(i,j)sqrt(d_(s,i)d_(s,j))R_(s,ij)^(k+1)
                     -s tr(M_b R_s^(circ k))-tr(Z_s R_s^(circ k)).

The sum is at least sum_i d_(s,i)=tr P_s because k+1 is even.
The last trace is o(n) uniformly in k by (3.1) and (5.1). Also

    |tr(M_o R_s^(circ k))|<=B_n ||M_o||_1=o(n),
    tr P_s=tr L_b+s tr M_b=n ell+o(n),

where tr M_b=-tr M_o and (2.1) were used. Combining with the ACTUAL
complete-source bound (3.4), we obtain a single a_n->0 such that

    (1/n)tr(L_b R_s^(circ k))>=ell-a_n
                        for both phases and every odd k>=3.          (5.2)

Since |M|-L_b=|M_o|>=0 and R_s^(circ k)>=0, (5.2) also holds with
the FULL |M| in place of L_b. For every t>=0, spectral calculus gives

    M^2=|M|^2>=t|M|-(t^2/4)I.

Pair with R_s^(circ k), whose trace is n, and take
t=2(ell-a_n)_+. Thus

    (1/n)tr(M^2 R_s^(circ k))>=(ell-a_n)_+^2=ell^2-o(1),
                         uniformly for both phases and odd k>=3.

For the actual local fields define

    F_(s,i)=(sM X_s)_i,
    v_(s,i)=(M C_(tail,s) M)_ii,
    b_(s,i)=kappa(M R_s M)_ii,
    c_(s,i)=E[(G_s)_i F_(s,i)],
    sigma_(s,i)^2=E F_(s,i)^2=b_(s,i)+v_(s,i),
    w_(s,i)=sigma_(s,i)^2-c_(s,i)^2>=v_(s,i).

Orthogonality of different Gaussian chaoses and first-chaos
Cauchy--Schwarz give the displayed relations. Summing the UNIFORM
trace bound with the nonnegative Hermite weights proves, for each s,

    (1/n)sum_i w_(s,i)>=(1-kappa)ell^2-o(1).               (5.3)

The remaining covariance controls also refer to actual fields:

    sigma_(s,i)^2<=B_n q=2C^2+o(1),       uniformly in i,s,
    (1/n)sum_i c_(s,i)^2
       <=(kappa/n)tr(M^2 R_s)
       <=kappa[C^2+B_n t_n]=kappa C^2+o(1).               (5.4)

Here M^2=M_b^2+M_o^2, M_b^2<=C^2 I, and tr R_s=n justify the last
bound. No first-chaos alignment or homogeneity of any variance is needed.

## 6. Importing the unchanged Gaussian gain through full-source smoothing

Fix an auxiliary update cutoff U=C+zeta with zeta>0. Equations
(5.4) ensure, for all sufficiently large n,

    sigma_(s,i)^2<=2U^2,       (1/n)sum_i c_(s,i)^2<=kappa U^2.

Use the continuous, nonnegative, globally 2-Lipschitz function

    g_U(r)=r^2/(2U)-(r-2U)_+^2/(2U),       r>=0.

Sections 6, 8, and 9 of the frozen 553-line all-law source give the
following gain for each phase, with r_(s,i)=(-X_(s,i)F_(s,i))_+:

    (1/n)sum_i E g_U(r_(s,i))
       >=(1-exp(-1))/(2U)
                         Psi(kappa U^2,(1-kappa)ell^2)-o(1)
        =D_U(ell)-o(1).                                   (6.1)

Its joint marginal Gaussianization applies to the ACTUAL coefficients
d_j=sM_ij: max|d_j|<=n^(-1/2) and sum_j d_j^2=q<=1. The frame norms
are bounded by (3.1). Rows of M_b are NOT used for this CLT. The radial
clipping argument uses the literal variance bound 2U^2. Convexity and
monotonicity of Psi, with 0<=partial_w Psi<=1/2, then use (5.3)--(5.4)
exactly as in that source. No new CLT, scalar integral, or Psi estimate
is being asserted, and the local parameters need not converge.

The actual operator penalty is bypassed by Section 6 of the frozen
325-line cap-free source. Here are its hypotheses and mean in the
present setting. For each phase and fixed epsilon>0 put

    X_(s,epsilon)=clip(G_s/epsilon,[-1,1]^n),
    Y_(s,epsilon)=clip((I+sM_b/U)X_(s,epsilon),[-1,1]^n).

The deterministic bulk clipping identity permits a nonzero diagonal
of M_b and has gain sum_i g_U((-X_(s,i)(sM_b X_s)_i)_+). Replacing
those bulk fields by the ACTUAL F_(s,i) costs o(n), since
E||M_o X_s||^2<=B_n tr M_o^2=o(n) and g_U is 2-Lipschitz.
The bulk phase baseline differs from its full baseline by o(n), since
|E Q_(sM_o)(X_s)|<=B_n ||M_o||_1/2=o(n).

The smoothed output map is odd and has Euclidean Gaussian-input
Lipschitz constant at most 2sqrt(B)/epsilon for a fixed eventual bound
B on B_n. Gaussian Poincare therefore gives

    Cov(Y_(s,epsilon))<=4B epsilon^(-2)I,
    |E Q_(sM_o)(Y_(s,epsilon))|
                            <=2B epsilon^(-2)||M_o||_1=o_epsilon(n).

The bounded bulk energy changes by at most O_U(n sqrt(epsilon))
between this smoothed mean and its unsmoothed mean. These are exactly
the estimates proved in Section 6 of that source; none uses flatness.
Finally independent Boolean rounding of the full original mean has
expected Q_(sM) equal to its mean energy because diag M=0. Thus

    alpha(A)>=e_s+(1/n)sum_i E g_U(r_(s,i))
                                  -O_U(sqrt(epsilon))-o_epsilon(1).

Take n->infinity for each fixed epsilon and then epsilon->0. Use
(6.1), and average the two valid oriented bounds with (4.1), to obtain

    liminf[alpha(A)-kappa/(2ell)-D_(C+zeta)(ell)]>=0.

By (2.1), ell lies in a compact subinterval of (0,1] eventually.
Uniform continuity of D_U there allows zeta->0 and proves (1.2).
This extra cutoff also avoids treating 2C^2+o(1) in (5.4) as a literal
finite bound of 2C^2. No covariance bound on the unsmoothed output is used.

## 7. Consequential region and the bulk-cutoff boundary

Suppose instead of (1.1) at one fixed cutoff that

    for EVERY fixed C'>5/3,
         (1/n)tr[M^2 1_(|M|>C')] ->0,
    limsup ell<=4/5.                                      (7.1)

Apply (1.2) at each C' and then let C' decrease to 5/3. The same
compactness and continuity argument gives the F_(5/3) lower at every
accumulation point of ell, which lies in [3/5,4/5] by (2.1).
The unchanged scalar region in Section 10 of the 553-line source gives

    liminf alpha(A)>=35/88+3/1250
                    =2/5+7/55000.                        (7.2)

Thus (7.1) excludes alpha(A)->2/5 for an entire actual all-law region.
No convergence of the full empirical spectral law is required in (7.1).

A sufficient spectral formulation of the TAIL premise is weak convergence
of that full law to a measure supported on [-5/3,5/3] whose second moment
is exactly 1.
For every fixed C'>5/3 its bounded-cutoff second moment tends to 1;
subtracting from tr M^2/n=q->1 proves (7.1)'s tail premise. The same
argument works for any bulk endpoint C>=1. The cutoff must be taken
STRICTLY above the endpoint before passing to the limit: atoms at C
approached from outside can prevent tail vanishing at C itself.

Rank alone does not directly control second-moment tails. For the actual
all-positive complete signing A=1_n 1_n^T-I, where 1_n is the all-ones
vector, M has one eigenvalue (n-1)/sqrt(n) and n-1 eigenvalues -1/sqrt(n).
Its tail rank is one at every fixed C>=1, but
tr[M^2 1_(|M|>C)]/n->1, not zero. Its alpha grows without bound.
The existing source norm bootstrap supplies the missing additional
control on bounded-alpha subsequences, as follows.

### Stronger bulk-rank corollary via the existing source bootstrap

The conclusion (7.2) in fact holds under the weaker hypotheses

    for EVERY fixed C'>5/3,
         rank 1_(|M|>C')/n ->0,
    limsup ell<=4/5.                                      (7.3)

No alpha cap is assumed in this corollary. To prove it, suppose there
were a subsequence violating (7.2) by a fixed positive amount. On that
subsequence alpha<=H for a fixed constant H. Section 2, equation (8),
of the existing norm-only phase-moment bootstrap identified in Section 9
then gives, without an actual operator cap,

    tr|A|^3<=(4/kappa+o_H(1))(n-1)Phi(A),
    (1/n)tr|M|^3<=(4/kappa+o_H(1))q alpha=O_H(1).          (7.4)

For a spectral projector P of rank r, finite Holder gives

    (1/n)tr(M^2 P)
       <=[(1/n)tr|M|^3]^(2/3)(r/n)^(1/3).

Thus (7.3) implies every tail premise in (7.1) on the violating
subsequence. The already proved (7.2) contradicts that violation.
This proves the stronger corollary. In particular a full weak limiting
law supported on [-5/3,5/3], together with limsup ell<=4/5, suffices;
no separate second-moment-saturation assumption is needed in this
stronger formulation. The support threshold itself is not inferred
from bounded alpha. Rank alone is still not a general matrix tail
estimate; (7.4) is the essential extra source input in this argument.

## 8. A global strict improvement of the original 1/pi lower bound

Let m_n be the minimum of Phi(A) over all complete symmetric
zero-diagonal signings A of order n. Then

    liminf_(n->infinity) m_n/n^(3/2)>1/pi.                 (8.1)

Equivalently, there exist epsilon0>0 and n0 such that every such
actual A with n>=n0 satisfies alpha(A)>=1/pi+epsilon0. This corollary
does NOT supply a numerical value of epsilon0 or prove convergence.

The previously proved exact nuclear baseline, recorded in Section 2
of the 553-line source, is alpha(A)>=kappa q/(2ell). Since
ell<=sqrt(q), it already gives liminf m_n/n^(3/2)>=kappa/2=1/pi.
If (8.1) failed, choose actual signings with n->infinity and
alpha(A)->kappa/2. The same nuclear bound forces ell->1. Hence

    (1/n)tr[(|M|-I)^2]=q-2ell+1 ->0.

For EVERY fixed C'>1, the scalar inequality

    x^2 1_(x>C')<=[C'/(C'-1)]^2(x-1)^2,       x>=0,

shows that (1/n)tr[M^2 1_(|M|>C')]->0 on this hypothetical sequence.
Apply the fixed-cutoff theorem at each C'>1, using ell->1, then let
C' decrease to 1. It gives

    liminf alpha(A)>=F_1(1)
       =kappa/2+(1-exp(-1))Psi(kappa,1-kappa)/2
       >kappa/2,

because 0<kappa<1 and the displayed gain is strictly positive. This
contradicts saturation and proves (8.1). No moment bootstrap is needed
for this argument; the displayed second-moment tail estimate is direct.

The value F_1(1) is a bound on the HYPOTHETICAL saturation sequence,
not an unconditional lower for all actual sources. In particular the
proof does NOT identify epsilon0 with F_1(1)-1/pi. Other source regimes,
the paired all-active-cell implication, convergence of the minima, and
the global original MO target remain OPEN.

## 9. Frozen imports, authorship, and full-read boundary

The unchanged joint Gaussianization, radial clipping, convex Psi
argument, and explicit scalar region are imported from

    original_mo_all_law_adaptive_nuclear_gain.md, 553 lines
    SHA256 0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51.

The full-original-source smoothed adaptive-mean transfer is imported
from Section 6 of

    original_mo_original_source_near_flat_cap_free_gain.md, 325 lines
    SHA256 0dfa5f62baaa57850a661bbc98d33d32440c783cccb11eaf5446feffbd81f7d4.

The norm-only cubic-moment estimate (7.4) is imported from the existing
canonical note, Section 2, equation (8),

    evidence/NOTE_2026-09-05_ORIGINAL_PHASE_MOMENT_BOOTSTRAP.md, 176 lines
    SHA256 3736db69d904b5a63ade46b32f6fddcc0505019f45ef483110c3ee67b24c8915.

All three complete frozen sources were directly read and their hashes
rechecked; their previously disclosed prerequisite provenance remains.
Root and the exact worker independently derived and checked the new
own-diagonal bulk frames, actual nuclear baseline, full higher-chaos
noise transfer, and actual covariance bounds before writing. Root
specified the endpoint-tail and asymptotic-error scope, the existing
bootstrap application, and the global strict-1/pi contradiction. The
docs-gate worker independently checked the bootstrap and global links;
the exact worker also checked them and authored this full extension,
including its complete author self-read.
These are contributions, not independent whole-new-source reviews.
Any new-source reviewer must disclose older prerequisite authorship.

This is one new /tmp source, with no separate author receipt. It and
all frozen imports remain unchanged after the announced final hash;
full-source independent review is a separate step. No canonical file,
actual source, covariance, paired block, or active-state definition was
changed by this author. Publication and documentation gates remain
root's workflow; no mathematical job, commit, or backup was performed here.
