# Scope freeze: one-vertex extension is the only live target

2026-09-25. Project-management decision for the original convergence problem.

The body of work has become too broad. From this point, a new theorem,
computation, construction, or structural lemma is in scope only if it
directly advances the one-vertex extension inequality below.

For an exact order-n minimizer A with F=m_n, define

    r_n = m_n [ (1+1/n)^(3/2) - 1 ].

Using the stable-state reduction, for an added row a in {+-1}^n,

    Phi(A extended by a)-F
      = max_[x in T_A] ( |a.x| - D_A([x]) ),

where T_A is the antipodal stable skeleton and
D_A([x])=F-E_A([x]) is the oriented deficit.

The sole live objective is to prove that there exists a row a such that

    max_[x in T_A] ( |a.x| - D_A([x]) )
      <= r_n + rho_n,                                       (1)

with an error satisfying

    sum_n rho_n/(n+1)^(3/2) < infinity.                     (2)

Then

    m_(n+1) <= m_n + r_n + rho_n

implies

    alpha_(n+1) <= alpha_n + rho_n/(n+1)^(3/2),

so alpha_n=m_n/n^(3/2) converges.

The strongest target is eventual neutral extension,

    boxed:
    m_(n+1) <= m_n + r_n,                                  (3)

which would make alpha_n eventually nonincreasing.

## What is now out of scope

The following are frozen unless they can be connected immediately and
quantitatively to (1):

- new lower bounds for liminf alpha_n;
- new principal-restriction inequalities by themselves;
- new two-block or doubling identities by themselves;
- new stable-pair geometry by itself;
- new Paley/conference structure by itself;
- new pressure/free-energy reformulations by themselves;
- optimizer classifications not yielding an extension row;
- generic discrepancy theorems not converted into the thresholded stable
  constraints in (1);
- computational searches that do not distinguish a concrete live mechanism
  for constructing or certifying the row a.

The value L=1/2 is also deferred. First prove that the limit exists. Only
after convergence is established should ratio-dense conference constructions
or other subsequence information be used to identify the limiting value.

## Acceptance rule for future work

A research turn counts as progress only if it does at least one of:

1. constructs a row a satisfying a quantitatively stronger version of (1);
2. proves a new bound on the thresholded stable family that reduces the left
   side of (1);
3. proves a summable rho_n;
4. converts an existing theorem into one of the previous three items;
5. produces a computation whose output is a certificate or discriminator for
   one of those items.

Duplicate structure, dead-end audits, and unrelated constant improvements do
not count as completed research turns.

This note does not claim a proof of convergence. It freezes the project onto
the shortest known route to one.
