# Counting entropy detects the normalized signing minimum

2026-09-13. An unconditional abundance bound and a sufficient convergence
criterion. Convergence of the counting entropy is still unproved.

This continues Section 1 of
`NOTE_2026-09-02_DIRECT_LIMIT_LITERATURE_AUDIT.md`. The probability bound
used below is already proved in Section 1 of
`NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md`; the new step is
to count the rounded signings and retain a bounded entropy at empty events.
The randomness is over edge signings, not over the spin states of one host.

## 1. A finite abundance bound around every signing

Let `n>=2`, `E=binom(n,2)`, and let A be any complete symmetric
zero-diagonal signing. All logarithms are natural. For `0<q<=1/2`, put

    h(q) = -q log q - (1-q) log(1-q),
    L_n = (n+1) log 2 + log n,
    T_n(q) = sqrt(8q(1-q) E L_n) + 4L_n/3.

Then the number of complete signings B satisfying

    Phi(B) <= (1-2q) Phi(A) + T_n(q)                         (1)

is at least

    exp(E h(q) - E log(2)/n - log(2)).                       (2)

To prove this, independently flip each edge of A with probability q.
The centered edge variables of `B-(1-2q)A` have magnitude at most two
and total variance `4q(1-q)E`. The Bernstein and union-bound argument in
the cited variational note, now with L_n as above, gives

    P(Phi(B-(1-2q)A) > T_n(q))
        <= 2^(n+1) exp(-L_n) = 1/n.                         (3)

Since `1-2q>=0`, every outcome outside this exceptional event satisfies
(1). Let G be the set in (1), let M be its cardinality, and let
`eta=P(B not in G)<=1/n`. In particular M is nonzero. The distribution
of B has Shannon entropy exactly `E h(q)`: its edge coordinates are
independent, and multiplying by the fixed signs of A is a bijection.
Condition on membership in G. The entropy chain rule and the bound by
the logarithm of the support size give

    E h(q)
      <= log 2 + (1-eta) log M + eta E log 2
      <= log 2 + log M + E log(2)/n.                        (4)

This proves (2). The conditional term with zero probability is taken as
zero. The bound uses the entropy of the whole product distribution; it
does not mistake a single high-probability atom for exponentially many
outcomes. No independence of the different spin energies is assumed.

## 2. Explicit positive entropy within any fixed tolerance of the minimum

Write `alpha_n=m_n/n^(3/2)`. For `0<delta<=1`, define

    q_delta = delta^2/64,
    N_delta = ceil(max(2, 256/(9 delta^2),
                       6 log(2)/h(q_delta))).              (5)

For every `n>=N_delta`,

    #{B: Phi(B) <= m_n + delta n^(3/2)}
        >= exp(E h(q_delta)/2)
        >= exp(n^2 h(q_delta)/8).                           (6)

Indeed, `log 2<=1`, `log n<=n-1`, and `E<=n^2/2` give `L_n<=2n` and

    T_n(q)/n^(3/2) <= sqrt(8q) + 8/(3 sqrt n).

At (5) this is at most

    delta/(2 sqrt 2) + delta/2 < delta.                     (7)

Choose an actual norm minimizer for A in (1). Its mean part has norm
`(1-2q_delta)m_n<=m_n`, so (7) puts every signing counted by (2) within
the tolerance in (6). Moreover, since `n>=2`,

    log(2)/n + log(2)/E <= 3 log(2)/n <= h(q_delta)/2.

This proves the first inequality in (6); `E>=n^2/4` proves the second.
The threshold N_delta is independent of the minimizing signing and of
any matrix-family hypothesis. It can be large; no computational density
claim at small n follows from this asymptotic application.

## 3. A bounded counting-entropy criterion for convergence

Define the actual labeled counts and their bounded entropy by

    C_n(t) = #{A: Phi(A)/n^(3/2) <= t},
    s_n(t) = log(1+C_n(t))/E_n,     E_n=binom(n,2).          (8)

The addition of one retains the empty event as `s_n(t)=0`. One always
has `0<=s_n(t)<=log 2 + log(2)/E_n`. Put

    a = liminf_n alpha_n,       b = limsup_n alpha_n.

The boundedness of alpha_n is supplied by CORE. Then

    a = inf{t: limsup_n s_n(t)>0},
    b = inf{t: liminf_n s_n(t)>0}.                          (9)

For `t<a`, the event in (8) is eventually empty. For `t>a`, choose
`0<delta<=1` with `a<t-delta`. Infinitely many n have
`alpha_n<=t-delta`, and (6) gives `s_n(t)>=h(q_delta)/2` on that
subsequence once n is large. These statements prove the first formula.
For `t<b`, infinitely many n have `alpha_n>t`, hence `s_n(t)=0` there.
For `t>b`, choose delta with `b<t-delta` and apply (6) for every
sufficiently large n. This proves the second formula. No assertion at
the boundary points t=a,b is needed.

Consequently, if `s_n(t)` converges for every t in a dense subset of
the real line, then alpha_n converges. Otherwise a<b, and a point t of
that dense set in `(a,b)` would have

    liminf_n s_n(t)=0 < limsup_n s_n(t),

contradicting the assumed convergence. This is a sufficient condition;
convergence of alpha_n does not by itself prove convergence of these
entire entropy curves.

## 4. Relation to the original probability-rate target

Under the uniform law on all `2^E_n` signings, put

    R_n(t) = -log P(Phi(A)/n^(3/2)<=t)/E_n,

with `R_n(t)=infinity` when the event is empty. Exactly,

    C_n(t)>0 => 0<=R_n(t)<=log 2,
    R_n(t)=log 2-log C_n(t)/E_n   when C_n(t)>0.            (10)

Thus convergence of R_n(t) on a dense set already suffices for
convergence of alpha_n, without the separate unique-threshold assumption
in the September 2 audit: if a<b, any t between them has an infinite-rate
subsequence and a subsequence bounded by log 2. They cannot have a common
extended-real limit.

The abundance theorem proves more. Uniformly for `n>=N_delta`,

    t>=alpha_n+delta => R_n(t)<=log 2-h(q_delta)/2.         (11)

It is (6), rather than the atom bound alone, that permits replacing the
infinite-valued rate by the bounded quantity s_n in Section 3. For an
arbitrary sequence of finite sets, alternating an empty accepted set
with a singleton makes `log(1+C_n)/E_n` tend to zero while nonemptiness
oscillates. The complete-signing noise bound rules out that isolated
optimum phenomenon at every fixed positive normalized tolerance.

## 5. Remaining implication and verification scope

The result supplies no comparison between s_n and s_N at different
orders. Establishing convergence of those bounded entropy curves is the
remaining sufficient implication on this route. Ordinary graphon limits
do not supply it, as already documented in
`NOTE_2026-09-02_THERMODYNAMIC_INTERPOLATION_GATE.md`.

Even monotonicity in n in either direction is false for the full finite
sequence: the previously certified values m_5=4,m_6=5,m_7=9 give
`alpha_6<7/20<alpha_5` and `alpha_6<2/5<alpha_7`. Therefore s_5(7/20)=0
and s_6(7/20)>0, whereas s_6(2/5)>0 and s_7(2/5)=0. No new enumeration
is needed for these exclusions. They do not refute eventual convergence
or an approximate comparison with a controlled error.

The source Bernstein argument was reused, not rerun. Root derived and
reviewed the entropy chain rule, the uniform quantifiers in (5)--(6), and
both directions of (9). Two symbolic scalar checks, recorded in
`bernoulli_counting_entropy_20260913/`, corroborate the new finite cutoff
arithmetic for all admissible real parameters. They do not formalize
the probability proof or the convergence criterion and are not an
independent human review. The original convergence problem remains OPEN.

The duplication check covered main, the nine other preserved local branch
heads, the three linked worktrees, and the unpublished analytic scratch
indexed in the September 6 campaign. Prior inputs are the September 2
rate formulation and the September 5 Bernstein rounding bound. No earlier
bounded counting-entropy threshold proof was found in that check.
