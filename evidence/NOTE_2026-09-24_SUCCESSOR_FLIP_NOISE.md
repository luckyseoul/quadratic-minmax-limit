# Independent flip noise forces an actual successor-field gap

2026-09-24. Supporting analytic results, author-reviewed in this session.
The original convergence problem remains OPEN. No independent human or
proof-assistant review is claimed. The new scalar checker has not run.

This note supplies a quantitative implication for an ACTUAL non-Gaussian
successor law. It uses only independent partial-flip coins and moments,
not a new Gaussian approximation. Combined with today's improved initial
changed-coordinate fraction, it also gives a further strict unconditional
lower-bound improvement.

## 1. Definitions and finite statement

Let A be a complete symmetric zero-diagonal signing of order n>=2,
M=A/sqrt(n), and alpha=Phi(A)/n^(3/2). Fix 0<p<=1/2. In each
orientation s, take ANY Boolean law X_s and put

    F_s=s M X_s,
    S_s={i:X_(s,i) F_(s,i)<0}, k_s=|S_s|,
    R_s=sum_(i in S_s) |F_(s,i)|.

Zero-field ties keep the old sign. Independently flip every coordinate of
S_s with probability p, leaving the other coordinates unchanged. Call the
actual Boolean successor X'_s. All coins are independent conditional on
the source state. Define phase-averaged quantities

    mu=(1/(2n))sum_s E k_s,
    e=(1/(2n))sum_s E Q_(sM)(X_s),
    f=(1/(2n))sum_s E ||M X_s||_1,
    g=f-2e=(1/n)sum_s E R_s,
    eta=(1/(2n^3))sum_s E[k_s^3-1_(k_s=1)].

Use primes for the successor. With

    d_p=2sqrt(p(1-p)/3), E_p=d_p+4/3,

the following finite, cap-free estimates hold:

    mu' >= p^2 mu,                                           (1)
    g' >= 2p d_p mu^(3/2)-2p E_p/sqrt(n),                     (2)
    e' >= (1-p^2)e+p(1-p)g-p^2 alpha,                        (3)
    eta' >= (mu')^3-1/n^3.                                  (4)

For even n, (1) improves to mu'>=p mu. The energy inequality (3)
and Jensen inequality (4) are reused; (1)--(2) are the new successor
inputs. Negative right sides in (2) cause no issue since g'>=0.

## 2. Two coins preserve a positive disagreement fraction

Condition on one source state x and its bad set S. For i in S, the new
field F'_i does not depend on i's own flip coin, because M_ii=0. Thus,
conditional on all other coins,

    P(x'_i F'_i<0 | other coins)>=p 1_(F'_i!=0),
    E[(-x'_i F'_i)_+ | other coins]>=p |F'_i|.               (5)

The second assertion follows from
(1-p)(-a)_++p a_+>=p|a| for a=x_i F'_i and p<=1/2.

If k=|S|>=2, choose j in S distinct from i and condition on every coin
except i,j. Changing j's coin changes F'_i by a nonzero number of
magnitude 2/sqrt(n). Its two possible values cannot both be zero.
Each coin value has probability at least p, so P(F'_i!=0)>=p.
The first inequality in (5) therefore gives P(i in S')>=p^2.

If k=1, the only original bad coordinate's field is unchanged by its
own coin. It stays bad with probability 1-p>=p^2. If k=0 there is
nothing to prove. Summing only over the original bad coordinates gives
E[k'|x]>=p^2 k; averaging proves (1). For even n every Boolean field
is a sum of n-1 odd-in-number signs divided by sqrt(n), so it is never
zero; (5) directly gives the stronger factor p.

## 3. An elementary shifted-noise moment lemma

Let Z=sum_j Z_j, with independent centered variables |Z_j|<=b, and
sigma^2=E Z^2. For any real shift a, set W=a+Z and
v=E W^2=a^2+sigma^2. Then

    E|W| >= sigma/sqrt(3)-2b/3.                             (6)

Indeed,

    |E Z^3|<=b sigma^2,
    E Z^4<=3sigma^4+b^2 sigma^2.

The second bound follows by expanding the fourth moment of independent
centered summands and using E Z_j^4<=b^2 E Z_j^2. Consequently,

    E W^4
      <=a^4+6a^2 sigma^2+3sigma^4+4|a|b sigma^2+b^2 sigma^2
      <=3v^2+4b v^(3/2)+b^2 v
      <=[sqrt(3)v+2b sqrt(v)]^2.                            (7)

Here |a|sigma^2<=v^(3/2). Holder interpolation gives
E|W|>=(E W^2)^(3/2)/(E W^4)^(1/2). If v>0, (7) implies

    E|W|>=v/[sqrt(3)sqrt(v)+2b]
          >=sqrt(v)/sqrt(3)-2b/3
          >=sigma/sqrt(3)-2b/3.

For the middle inequality use u^2/(u+c)>=u-c with u=sqrt(v),
c=2b/sqrt(3). If v=0, (6) is immediate. This argument is uniform
in the shift a; no centered-field or Gaussian-field assumption is made.

## 4. Apply the moment lemma to actual updated fields

Condition on x. For i in S, write the field exactly as

    F'_i=a_i+sum_(j in S, j!=i) u_ij (xi_j-p),
    u_ij=-2s M_ij x_j in {-2/sqrt(n),+2/sqrt(n)},
    xi_j independently Bernoulli(p).

The shift a_i is arbitrary. The centered summands have absolute value
at most b=2/sqrt(n), and

    sigma_i^2=4p(1-p)(k-1)/n.

Equations (6) and sqrt(k-1)>=sqrt(k)-1 give, for k>=1,

    E[|F'_i| | x]
      >=d_p sqrt((k-1)/n)-4/(3sqrt(n))
      >=d_p sqrt(k/n)-E_p/sqrt(n).                         (8)

Apply the second inequality in (5) and sum over i in S only:

    E[R'|x]/n
      >=p d_p (k/n)^(3/2)-p E_p (k/n)/sqrt(n).              (9)

Jensen over the phase/state mixture and k/n<=1 prove (2). This is a
finite estimate for the actual successor law, not for a replacement law.

For completeness, multilinearity and independent flips give

    E[Q_(sM)(X')|x]
      =(1-p)^2 Q_(sM)(x)+p(1-p)||Mx||_1+p^2 Q_(sM)(Y),

where Y is the full best response. Since Q_(sM)(Y)>=-alpha n,
averaging gives (3). Jensen, retaining the k'=1 correction, gives (4).

## 5. A further strict unconditional lower bound

Start from the tilted paired Gaussian source at t=993/1000. Reuse
[the both-phase theorem](NOTE_2026-09-24_BOTH_PHASE_DISAGREEMENT.md):
mu_0>=1/5-o_L(1) for each fixed cap ||M||op<=L. As before set

    kappa=2/pi, z=t^2/(1+t^2), a=1-kappa asin(z)+kappa z,
    e0=kappa t/(1+t^2), f0=sqrt(kappa a),
    p=(f0-2e0)/(f0+sqrt(f0^2+(f0-2e0)^2)),
    B=e0-f0/2+sqrt(f0^2+(f0-2e0)^2)/2.

Here e0,f0 denote limiting lower inputs, not exact finite-source means.
Use this same p for ONE actual partial-flip update. Define

    G=2p d_p/5^(3/2),  G^2=16p^3(1-p)/375>0.

The new field estimate and the reused energy estimate give

    g_1>=G-o_L(1),
    e_1>=(1+p^2)B-p^2 alpha-o_L(1).                       (10)

The [sharp deficit-disagreement law](NOTE_2026-09-23_POINTWISE_DEFICIT_DISAGREEMENT.md)
applies to this actual successor, without a Gaussian hypothesis. For ANY
law with D=alpha-e and g=f-2e it implies

    D>=alpha+g/2-sqrt(alpha^2+alpha g)
      =g^2/[4(alpha+g/2+sqrt(alpha^2+alpha g))]
      >=g^2/(8alpha+4g).                                  (11)

The first bound is the inverse of g/2<=D+sqrt(2alpha D).
The last expression increases with g>=0. Applying it to (10) yields

    (1+p^2)(alpha-B)>=G^2/(8alpha+4G)-o_L(1),
    (alpha-B)(alpha+G/2)>=G^2/[8(1+p^2)]-o_L(1).           (12)

The positive root is

    B_noise=
      [B-G/2+sqrt((B+G/2)^2+G^2/[2(1+p^2)])]/2.            (13)

Thus liminf alpha_n>=B_noise after the same-order cap removal. Explicitly,
apply the fixed-cap proof to the regularized minimizers with L=K+8,
take n to infinity first, and obtain
liminf alpha_n>=B_noise-2sqrt(2 Gamma/K), where
Gamma=4pi/log(1+sqrt(2)). Then take K to infinity. No growing-cap
uniform Gaussianization is assumed. The successor estimates themselves
are cap-free; only their initial Gaussian input needs regularization.

### Exact comparison with the new interior-only bound

Let x=B_both from the preceding note. It satisfies

    x(x-B)=p^2/[3000(1+p^2)],  x>1/4.

For the last inequality, B>=e0>(7/11)(993/1000)/2>1/4.

The polynomial in the first line of (12), after clearing its positive
denominator, is negative at x precisely when

    [32p(1-p)-2]x>G.                                      (14)

This comparison can be certified without decimal optimization. The reused
rational source bound is

    (2e0/f0)^2<2112/3229<(119/143)^2.

Since (f0-2e0)/f0=2p/(1-p^2), this proves p>1/12. Also p<1/2.
Hence p(1-p)>11/144 and the left side of (14) is greater than
(4/9)(1/4)=1/9. Meanwhile p^3(1-p)<=1/16 gives
G^2<=1/375<1/81, so G<1/9. This proves (14), and therefore

    B_noise > B_both > B_int > B_tilt.                     (15)

We compare two valid gaps; we do NOT add the interior gap to the successor
gap. Both constrain the same successor energy, so adding them without a
further argument would double count. Expression (13), rather than an
unverified rounded decimal, is the canonical new lower constant.

## 6. Fixed-time repeated-update consequences and what is still missing

Iterate the same state-independent p with fresh independent coins at each
step. No Gaussianization of any successor is needed. Equations (1)--(4)
give explicit one-step inequalities for the actual sequence of laws.
In particular,

    mu_j>=p^(2j) mu_0,
    eta_j>=p^(6j) mu_0^3-1/n^3,
    g_(j+1)>=2p d_p p^(3j) mu_0^(3/2)-2p E_p/sqrt(n).       (16)

For the initial tilted source and every FIXED j,

    eta_j>=p^(6j)/125-o_L(1),
    g_(j+1)>=G p^(3j)-o_L(1).                             (17)

These floors decay exponentially with j. The lower energy recursion (3)
does not prove energy monotonicity or approach to alpha; a fixed p can
lose energy on some states. The proof does not cover a state-dependent
choice of p without additional lower control on those probabilities.
The finite recursions (16) hold at every j, but their shrinking floors
do not provide a useful nondecaying estimate or energy-convergence theorem
as j diverges. No cross-order comparison has been proved. The original
limit therefore remains OPEN.

The remaining target is sharper trajectory control with a useful
long-time energy conclusion, or an independent cross-order argument;
it is no longer accurate to say there is no quantitative actual-successor
estimate at all.

## 7. Provenance and verification

Duplicate checks covered the active checkout, all already-fetched remote
branch heads, the active source/update proofs, and the archived
UNPUBLISHED_WORKING_FRONTIER.md indexed by ARTIFACTS.md. Only one local
worktree is exposed. The preserved archive's cubic posterior-overlap
estimate concerns a different comparison and is not this flip-noise law.
No family scan, signing census, or old computation was restarted.

Reused results: the mean-update energy identity, sharp pointwise deficit
closure, phase/state Jensen inequality, source bounds, and same-order
cap removal. New implications: (1), the shifted-noise estimate (6) as
applied in (2), and their actual-successor and lower-bound consequences.

Analytic author review checked the zero-field cases, absence of a diagonal
term, conditional coin independence, uniformity in the field shift, fourth
moment bounds, phase normalizations, positive-root algebra, comparison
(14), and the fixed-cap order of limits. The focused checker is
`successor_flip_noise_20260924/check.py`; it remains UNEXECUTED. SSH
preflight to NUKA failed with a socket permission error. No mathematical
test was run on the controller, and no independent review is claimed.
These are supporting same-order results, not a completed convergence
milestone; no new major-milestone backup is claimed.
