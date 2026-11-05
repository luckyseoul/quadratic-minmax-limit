# Experiment 1 — one-sided all-future modulus (finite diagnostic)

## Question

Does \(\alpha_n = F(n)/n^{3/2}\) converge? (Ivanisvili, MO 413935; still open.)

## Source of the attack

`luckyseoul/quadratic-minmax-limit` HANDOFF 2026-11-04: a one-sided all-future
upper comparison is sufficient for convergence. If there exists
\(\varepsilon_n\to 0\) such that for all large \(n\) and all \(m\ge n\),

\[
\alpha_m \le \alpha_n + \varepsilon_n,
\]

then \(\liminf\alpha_n=\limsup\alpha_n\). Consecutive-gap control
\(O(n^{-1/2})\) alone is not enough unless the errors sum to a uniform
all-future bound.

## Working hypothesis (H1)

On the exact table \(F(1),\dots,F(15)\) (OEIS A399208 / both research
archives), the finite one-sided excess

\[
E_n(N) := \max_{n\le m\le N}\alpha_m - \alpha_n
\]

with \(N=15\) is eventually decreasing in \(n\) for \(n\ge 8\), and is
dominated by the candidate schedule \(\varepsilon_n = n^{-1/2}\)
(i.e. \(E_n(15)\le n^{-1/2}\) for all \(n\ge 8\)).

Extending with published upper envelopes
\(F(16)\le 30\), \(F(17)\le 32\), and the analytic ceiling \(\alpha_m\le 1/2\),
the mixed excess

\[
\overline E_n := \max\bigl(E_n(15),\;
\tfrac{30}{16^{3/2}}-\alpha_n,\;
\tfrac{32}{17^{3/2}}-\alpha_n,\;
\tfrac12-\alpha_n\bigr)
\]

still satisfies \(\overline E_n \le n^{-1/2}\) for \(n\ge 8\).

A failure of the \(\varepsilon_n=n^{-1/2}\) domination on exact data would
show that this particular schedule is already too tight at finite \(n\);
it would not prove nonconvergence. Success is only a consistency check.

## Method

1. Load exact \(F(n)\) for \(n\le 15\).
2. Verify \(F(n)\) exhaustively for \(n\le 6\) (Seidel switching).
3. Compute \(\alpha_n\), \(E_n(15)\), \(\overline E_n\), and
   \(E_n(15)/n^{-1/2}\).
4. Record pass/fail for each clause of H1.
5. Write `results/experiment1.json` and `results/experiment1_series.csv`.

## Expected result

- Exact verification matches A399208.
- \(E_n(15)\) is small for large \(n\) but **not** monotone in \(n\)
  (because \(\alpha_n\) oscillates; \(m_n\) is monotone, \(\alpha_n\) is not —
  CORE.md §2).
- The schedule \(\varepsilon_n=n^{-1/2}\) likely covers exact \(E_n(15)\) for
  \(n\ge 8\); the analytic \(1/2\) ceiling may dominate \(\overline E_n\) and
  keep it from shrinking like \(n^{-1/2}\) until better uniform upper bounds
  exist.

## Analysis

Tabulate \(E_n\), compare to \(n^{-1/2}\) and to the conference gap
\(1/2-\alpha_n\). Report which H1 clauses hold.

## Limitations

- Finite \(N=15\) cannot prove \(\varepsilon_n\to 0\) for all future \(m\).
- Upper envelopes in \(\overline E_n\) can be loose; a large \(\overline E_n\)
  may reflect weak constructions rather than true \(\alpha_m\).
- Does not settle existence of \(\lim\alpha_n\).
