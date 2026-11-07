# Exact mixed-center, mixed-temperature BSC bridge accounting

2026-11-07. New scoped algebraic lemma extending the November 4
same-center BSC identities. The original quadratic-minmax limit remains OPEN.
This is **not** a cross-order comparison for the actual minima.

## Exact theorem

Let E be a finite coordinate set. Fix k>=2, temperatures beta_r>=0
(r=1,...,k), centers x^(r) in {+-1}^E, and one shared output
y in {+-1}^E. Write

    T = sum_r beta_r,
    K_b(y|x) = exp(b*x*y)/(2 cosh b),
    R_b(y|x) = (1/2)/K_b(y|x) = cosh(b) exp(-b*x*y),
    Gamma = log cosh(T) - sum_r log cosh(beta_r) >= 0.

At each e choose a **weighted majority center** z_e in {+-1} maximizing
z_e sum_r beta_r x^(r)_e; ties may be resolved arbitrarily. Define

    w_e = sum_{r: x^(r)_e != z_e} beta_r,
    W = sum_e w_e.

Then for every output y the following identity is EXACT:

    sum_e log R_T(y_e|z_e)
      - sum_r sum_e log R_(beta_r)(y_e|x^(r)_e)
      = |E| Gamma - 2 sum_e w_e y_e z_e.                 (1)

Consequently,

    sum_e log R_T(y_e|z_e)
      >= sum_r sum_e log R_(beta_r)(y_e|x^(r)_e)
         + |E| Gamma - 2W.                               (2)

The lower bound is sharp over output vectors: equality holds when y=z
on every coordinate with w_e>0. Therefore the pointwise reverse-budget
domination, for **every** shared output y, holds if and only if

    2W <= |E| Gamma.                                      (3)

This permits unequal centers and unequal temperatures; it does not
require exact synchronization.

## Proof

For each e,

    sum_r beta_r x^(r)_e = (T-2w_e) z_e.

Using the definition of R and subtracting its log expressions gives

    log R_T(y_e|z_e)
      - sum_r log R_(beta_r)(y_e|x^(r)_e)
      = Gamma - y_e [T z_e - sum_r beta_r x^(r)_e]
      = Gamma - 2w_e y_e z_e.

Sum over e for (1). Since y_e z_e<=1, (2) follows, with equality
for y=z. Finally cosh(a+b)>=cosh(a)cosh(b) for a,b>=0,
so Gamma>=0. This proves (3) and all claims.

## Equal-temperature corollary and sharp tolerance

If beta_r=beta, let h_e be the number of centers disagreeing with
the majority center z_e, H=sum_e h_e. Then W=beta H and

    Gamma = log cosh(k beta) - k log cosh beta.

The exact worst-output reverse-budget margin is

    |E| Gamma - 2 beta H.                                 (4)

For beta>0, full pointwise domination holds exactly when

    H <= |E| Gamma/(2 beta).                              (5)

As beta -> infinity, Gamma -> (k-1) log 2. In particular, at
beta=log n, the sharp tolerated mismatch fraction is asymptotic to

    (k-1) log 2/(2 log n).                                (6)

This is a genuine relaxation of the exact same-center hypothesis in
NOTE_2026-11-04_REVERSE_DINF_BSC_TENSORIZATION.md. A constant positive
fraction of mismatched shared coordinates, however, fails the uniform
pointwise comparison at sufficiently high temperature.

## Forward-likelihood caution

For forward surprisal S_b=-log K_b, the corresponding exact identity is

    sum_e S_T(y_e|z_e) - sum_r sum_e S_(beta_r)(y_e|x^(r)_e)
      = -|E| s - 2 sum_e w_e y_e z_e,                    (7)
    s = (k-1) log 2 - Gamma >= 0.

Thus the worst-output forward budget loses precisely |E|s+2W.
For equal beta=log n and |E|=O(n^2), |E|s=O(1) for each
fixed k, while the mismatch penalty is 2(log n)H. To retain an o(n)
forward loss by this pointwise certificate, it suffices that
H=o(n/log n). This is a sufficient condition for this certificate,
not a necessary condition for convergence or for every other bridge.

## Scope and deduplication

The closest existing notes are:

- NOTE_2026-11-04_EXACT_BSC_SYNCHRONIZATION.md (same center, forward);
- NOTE_2026-11-04_REVERSE_DINF_BSC_TENSORIZATION.md (same center, reverse);
- NOTE_2026-11-04_HYPERCONTRACTIVE_ROW_COLUMN_BRIDGE.md (noise-smoothed
  row/column factors).

Searches of this repository and Robby955/mo-413935-research found no
matching explicit weighted-majority mismatch budget. This statement is
proved by direct algebra, not by a finite computation.

The source definitions of MOP-DINF-CLOSURE have **not** been identified
with the BSC factors in (1). Shared output coordinates are still assumed,
even though equal centers and temperatures are not. Equations (1)-(7)
do not construct a signing of order m from one of order n, do not prove
alpha_m <= alpha_n+o(1), and do not alter src/original_mo_status.py.
