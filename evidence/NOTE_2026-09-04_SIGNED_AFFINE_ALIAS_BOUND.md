# Proposition 15.763: signed affine-alias incidence bound

**Status:** proved conditional theorem for every admissible odd affine
parameter. It strengthens Proposition 15.755 on its affine-alias branch,
including the first shell at `r=1`. It does not classify the full first defect
shell and does not close residual (ii), E1, `L=1/2`, or the original
MathOverflow limit.

## Setup

Use Proposition 15.755's dangerous shared-maximizer normalization. Let `G`
have even cardinality, let `e` be outside `G`, put

\[
 H=G\cup\{e\},\qquad A=C\mathbin\triangle G,
 \qquad B=C\mathbin\triangle H,
\]

and assume

\[
 \Phi(A)=\Phi-2,\qquad \Phi(B)\le\Phi-4.
\]

Choose a signed maximizer `x` of `A`, with phase
\(\epsilon\in\{+1,-1\}\). The one-edge comparison gives

\[
 \epsilon f_e(x)=1,
 \qquad \epsilon S_H(x)=2-\frac{\delta_\epsilon(x)}2.       \tag{1}
\]

Suppose this particular `x` is an odd-parameter affine alias. Thus, for a
positive odd integer `r`,

\[
 \delta_\epsilon(x)=2pr^2,
 \qquad m={p+1\over2}+r\le p,                               \tag{2}
\]

and `x` is constant on `m` positive parallel fibres. Flipping the union
\(T_J\) of any `r` of those fibres gives an
\(\epsilon p\)-Boolean eigenvector

\[
 w_J=x^{T_J},\qquad |J|=r.                                  \tag{3}
\]

There are

\[
 N={m\choose r}                                             \tag{4}
\]

such aliases. This is a hypothesis on the selected active cube point, not a
classification of every point at defect `2pr^2` (or even every first-shell
point when `r=1`).

## Signed alias cuts

For an edge \(h=\{u,v\}\in H\), put

\[
 a_h=\epsilon C_{uv}x_ux_v\in\{+1,-1\}.
\]

Equations (1)--(2) give the exact signed row sum

\[
 \sum_{h\in H}a_h=2-pr^2.                                  \tag{5}
\]

For an alias cut define

\[
 L_J=\sum_{\substack{h\in H\\h\text{ crosses }T_J}}a_h.
\]

Switching `x` on \(T_J\) changes the sign of exactly those terms, so

\[
 \epsilon S_H(w_J)=2-pr^2-2L_J.                            \tag{6}
\]

Because \(w_J\) is on the Paley eigenshell,

\[
 \epsilon Q_B(w_J)=\Phi-2\epsilon S_H(w_J)\le\Phi-4.
\]

Therefore \(\epsilon S_H(w_J)\ge2\). The set `H` is odd, so its signed
score is odd and in fact

\[
 \epsilon S_H(w_J)\ge3,
 \qquad -L_J\ge{pr^2+1\over2}.                             \tag{7}
\]

## Retaining the signs

Let \(\mu_h\) count the alias cuts crossed by `h`. An edge between two
distinct positive fibres crosses

\[
 M=2{m-2\choose r-1}                                       \tag{8}
\]

alias cuts. An edge between a positive fibre and the outside crosses
\({m-1\choose r-1}\le M\), and the other edge types cross no more. Hence
\(\mu_h\le M\) for every edge.

If \(N_-\) is the number of edges with \(a_h=-1\), (5) gives

\[
 N_-={|H|+pr^2-2\over2}.                                   \tag{9}
\]

Summing (7) and retaining the signs, rather than replacing every summand by
its absolute value, gives

\[
 {N(pr^2+1)\over2}
 \le \sum_{h\in H}\mu_h(-a_h)
 \le M N_-.                                                \tag{10}
\]

The positive edges can only decrease the middle expression. Substitute
(4), (8), and (9), and use

\[
 {{m\choose r}\over {m-2\choose r-1}}
 ={m(m-1)\over r(m-r)}.
\]

This proves

\[
 \boxed{
 |H|\ge
 { (pr^2+1)m(m-1)\over 2r(m-r)}-pr^2+2.}                  \tag{11}
\]

The right side is rational, while `H` is odd. The exact integral statement
is that `|H|` is at least the least odd integer greater than or equal to the
right side. The separate trivial consequence of (5),
\(|H|\ge pr^2-2\), and the old Proposition 15.755 bound remain available;
the effective lower bound is the maximum of all three.

Proposition 15.755 discarded signs in (10) and obtained

\[
 B_{15.755}={ (pr^2+1)m(m-1)\over4r(m-r)}.
\]

The exact comparison is

\[
 B_{\rm signed}=2B_{15.755}-pr^2+2.                        \tag{12}
\]

The signed expression need not dominate for every large admissible `r`, so
the implementation retains both bounds. At `r=1`, however,

\[
 B_{15.755}={(p+1)(p+3)\over8},\qquad
 \boxed{B_{\rm signed}={p^2+11\over4}}.                    \tag{13}
\]

For odd `p`, \(p^2\equiv1\pmod8\), so the latter is already an odd integer.
It is strictly larger than the parity-adjusted old bound for every `p>=5`.

## Critical-alias alternative

The active edge has \(a_e=1\). If `e` does not cross \(T_J\), then
\(\epsilon f_e(w_J)=1\); if it crosses, that sign is `-1`. Consequently an
internal alias with \(\epsilon S_H(w_J)=3\) gives

\[
 \epsilon S_{H\setminus\{e\}}(w_J)=2.                     \tag{14}
\]

Suppose no alias has both properties in (14). Every noncrossing alias then
has odd H-score at least five, improving (7) by one. If \(\mu_e\) aliases
cross `e`, the summed lower bound gains \(N-\mu_e\), while the positive
edge `e` contributes exactly \(-\mu_e\) to the signed upper bound. These
two occurrences cancel, leaving

\[
 {N(pr^2+3)\over2}\le M N_-.
\]

Thus either (14) occurs, or

\[
 \boxed{
 |H|\ge
 { (pr^2+3)m(m-1)\over 2r(m-r)}-pr^2+2,}                  \tag{15}
\]

again rounded up to an odd integer. At `r=1` the rational expression in
(15) is \((p^2+2p+17)/4\). This dichotomy identifies exactly what the affine
family contributes to the missing eigenshell bridge; it is not an exclusion
of residual (ii).

## Scope and replay

The bad-sign equality case of Proposition 15.755 yields an integral
eigenvector with one coordinate of magnitude three. The affine fibre family
is an explicit subfamily of that case, not a proved classification. In an
all-deletions minimal graph, the active maximizer and its affine coordinates
may also vary with the deleted edge. Equations (11) and (15) therefore cannot
be aggregated across deletions without a new common-coordinate or
minimum-shell classification theorem.

Exact replay:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q tests/test_prop15763.py
PYTHONPATH=src /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_prop15763.py
```

The source uses `fractions.Fraction` throughout. It performs no finite-prime,
graph, orbit, solver, or eigenshell census. Residual (ii), E1, `L=1/2`, and
the original MathOverflow limit remain **OPEN**.
