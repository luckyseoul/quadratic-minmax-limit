# Core-signature pairing for one-vertex extension

2026-09-25. Direct one-vertex extension theorem. This is in scope under
`SCOPE_FREEZE_2026-09-25_ONE_VERTEX_ONLY.md`: it constructs an actual
incident sign row from a strictly more correlation-sensitive certificate
than the reciprocal-square mass `Xi_A(r)`.

DUPLICATION AUDIT: the core-signature pairing itself is NOT new. Section 7
of `NOTE_2026-09-01_ONE_TWO_BIBALANCED_HADAMARD.md` already partitions
coordinates by complete relative-sign signatures, pairs inside each class,
annihilates every fixed anchor on the paired core, and records the
`2^(R-1)` leftover bound. The new content here is specifically the
one-vertex conversion: after that existing annihilation step, apply the
Guo--Fang--Lu Komlos theorem to the residual pair-difference vectors with
state-dependent weights `1/(D_s+r-u)`, yielding the bottleneck residual
mass criterion (2)--(4). Do not count the signature-pairing mechanism itself
as a new result.

Let `A` be a complete symmetric zero-diagonal signing of order `n`,
put `F=Phi(A)`, and let `T_A` be the antipodal stable skeleton from
`NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md`. For each
`s in T_A`, choose a Boolean representative `x^s`, write its stable
deficit as `D_s`, and for a desired increment `r>0` put

    t_s = D_s + r.

As usual, states with `t_s>n` are automatically harmless and may be
discarded. Let `H` denote the active family `t_s<=n`.

## Theorem: annihilate a core, balance only residual pair disagreements

Choose any core `C subset H`. For coordinate `i`, define its core
signature

    sigma_i = (x_i^c)_(c in C) in {+-1}^C.

Declare coordinates `i,j` projectively equivalent when

    sigma_j = lambda sigma_i

for some `lambda in {+-1}`. Inside every projective signature class,
pair coordinates arbitrarily except for at most one leftover coordinate.
Let `P` be the resulting set of disjoint pairs and `U` the leftover
coordinates, with

    u = |U|.

For a pair `p=(i,j)`, let `lambda_p` be the sign satisfying
`sigma_j=lambda_p sigma_i`. Define the residual weighted disagreement mass

    M_p(r,u)
      = sum_(s in H\C : x_i^s != lambda_p x_j^s)
          1/(D_s+r-u)^2.

If

    u < r                                                     (1)

and the pairs can be chosen so that

    4 max_(p in P) M_p(r,u) <= 1/(18 pi),                    (2)

then there exists an actual incident row `a in {+-1}^n` such that

    Phi(A extended by a) < F+r.                              (3)

Equivalently, it is enough that every selected pair obey

    M_p(r,u) <= 1/(72 pi).                                   (4)

The theorem is finite and deterministic.

## Proof

For every pair `p=(i,j)`, impose the relation

    a_j = -lambda_p a_i.

Write the remaining free pair sign as `eps_p=a_i`. Give the leftover
coordinates in `U` arbitrary fixed signs.

For a core state `c in C`,

    x_j^c = lambda_p x_i^c,

so every paired contribution cancels exactly:

    a_i x_i^c + a_j x_j^c
      = eps_p (x_i^c-lambda_p x_j^c)
      = 0.

Therefore

    |a.x^c| <= u < r <= D_c+r.                              (5)

Now fix a residual active state `s in H\C`. The leftover coordinates
contribute an offset `b_s` with `|b_s|<=u`. Put

    tau_s = D_s+r-u >0.

For every pair `p=(i,j)`, form the vector indexed by residual states

    v_p(s)
      = (x_i^s-lambda_p x_j^s)/tau_s.

Its entries are `0,+-2/tau_s`, and

    ||v_p||_2^2
      = 4 sum_(s in H\C : x_i^s != lambda_p x_j^s)
          1/tau_s^2
      = 4 M_p(r,u).                                         (6)

By (2), every `v_p` has Euclidean norm at most
`1/(3 sqrt(2 pi))`. Apply the Guo--Fang--Lu Komlos theorem to the
scaled vectors `3 sqrt(2 pi) v_p`. There are pair signs `eps_p` such
that simultaneously for every residual state `s`,

    |sum_p eps_p v_p(s)| < 1.                               (7)

Hence the paired part of `a.x^s` has magnitude strictly below
`tau_s`. Adding the leftover offset gives

    |a.x^s|
      < tau_s+u
      = D_s+r.                                              (8)

Inactive states have `D_s+r>n>=|a.x^s|`, so they are harmless as well.
The exact stable-state extension identity then yields (3). QED.

## Why this is different from the old Xi certificate

The reciprocal-square theorem used one vector for every OLD coordinate with

    ||v_i||_2^2 = Xi_A(r),

so every active stable state charged every coordinate equally, regardless
of correlation among the constraints.

The present theorem first quotients the coordinate set by a chosen core.
Every selected pair is then identically zero on the entire core, and on the
residual family it pays only for states that actually separate that pair.
Thus the relevant quantity is no longer

    sum_s 1/(D_s+r)^2,

but the maximum, over the chosen coordinate pairs, of a weighted
DISAGREEMENT sum.

This can remain small even when `Xi_A(r)` is huge because many dangerous
stable states are mutually redundant as constraints.

The cost of killing the core is only the number `u` of odd projective
signature classes, not `|C|`. In particular

    u <= number of nonempty projective signatures <= 2^(|C|-1)

for nonempty `C`.

## Neutral and summable-error consequences

For an exact order-`n` minimizer, take

    r_n = m_n[(1+1/n)^(3/2)-1].

If some core and projective-signature pairing satisfy (1)--(2) with
`r=r_n`, then

    alpha_(n+1) < alpha_n.

More generally replace `r_n` by `r_n+rho_n`. If the theorem's
certificate holds and

    sum_n rho_n/(n+1)^(3/2) < infinity,

then the existing one-vertex argument gives convergence of `alpha_n`.

## Exact finite discriminator

This theorem suggests a concrete search that is directly relevant to the
frozen target:

1. enumerate only the active stable classes and their deficits;
2. choose a small core `C`;
3. build projective coordinate signatures on `C`;
4. pair coordinates within each signature class;
5. minimize the bottleneck residual mass `max_p M_p`.

A successful output is itself an extension-row certificate, not merely a
structural statistic.

For the recorded exact order-15 minimizer, the old neutral Komlos quantity
is extremely loose: there are 301 active antipodal stable classes and
`Xi_A(r_15) approximately 14.1414`, while `1/(18 pi) approximately
0.01768`. This explains why a correlation-sensitive certificate is needed;
the known exact one-vertex optimum nevertheless increases `Phi` only from
27 to 30.

This note does not claim that (1)--(2) hold asymptotically for all exact
minimizers. It converts stable-state redundancy into a quantitative
one-vertex mechanism and gives a finite certificate to hunt next.
