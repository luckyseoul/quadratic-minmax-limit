# The first post-third even candidate shell is empty

Date: 2026-08-25. Proposition 15.638. Retain the Paley dual lattice
\(L^*=P\mathbb Z^n\), the \(R=(p+1)/2\) square-direction profiles
\(a_{j,s}\), their common sum \(t\), and the scaled norm

\[
s=2p\lVert x\rVert^2.
\]

Proposition 15.637 excludes the zero-common-sum branch at the first
post-third even candidate value \(s=2(p+3)\). This note excludes every remaining
common sum and proves

\[
\boxed{\{x\in L^*:2p\lVert x\rVert^2=2(p+3)\}=\varnothing
       \qquad(p\ge11).}
\]

It does not identify the next nonempty shell or control the full harmonic
theta tail required by R1. In particular, the separate odd-phase floor
\(3p-6\) is not being reclassified here; at \(p=11\) it lies below
\(2(p+3)\).

## Balancing leaves three nonzero sums

For \(|t|=ap+b\), \(0\le b<p\), let

\[
f_p(t)=(p-b)a^2+b(a+1)^2.
\]

The circle-frame identity and profile balancing give

\[
\sum_{j,s}a_{j,s}^2={s+t^2\over2},\qquad
g_p(t)=(p+1)f_p(t)-t^2
      =pa^2+2ab+b(p+1-b)\le s.                   \tag{1}
\]

Here \(s\) is even, so \(t\) is even. If \(a\ge2\), then
\(g_p(t)\ge4p>2(p+3)\). Checking \(0\le |t|<2p\) in (1) leaves exactly

\[
|t|\in\{0,2,p-1,p+1\}.                           \tag{2}
\]

The \(t=0\) case is Proposition 15.637. Replacing \(x\) by \(-x\) handles
negative \(t\), so consider the three positive values.

Put \(k=(p-1)/2\), so \(R=k+1\). Relative to balanced profiles, the total
excess squared energies are

\[
\begin{array}{c|c|c}
t&\text{balanced profile}&\text{total excess}\\ \hline
2&\delta_a+\delta_b&4\\
p-1&\mathbf1-\delta_a&4\\
p+1&\mathbf1+\delta_a&2.
\end{array}                                      \tag{3}
\]

For \(1\le d\le4\), profile glue supplies a degree-\(d\) binary form
\(q_d\) whose values on the selected projective directions are the
\(d\)-th power sums. These four degrees are available because
\(R-2\ge4\) when \(p\ge11\).

## The \(t=2\) branch

An ordinary profile is an unordered pair
\(\delta_a+\delta_b\) with \(a\ne b\). Its moments satisfy

\[
\begin{aligned}
2q_3-3q_1q_2+q_1^3&=0,                           \tag{4}\\
2q_4-q_2^2-2q_1^2q_2+q_1^4&=0.                  \tag{5}
\end{aligned}
\]

The local integral profiles of sum two and energies four and six are,
up to permuting roots,

\[
\begin{array}{c|c}
4&2\delta_a,\quad
  \delta_a+\delta_b+\delta_c-\delta_d\\
6&2\delta_a+\delta_b-\delta_c,\quad
  \delta_a+\delta_b+\delta_c+\delta_d-\delta_e-\delta_f.
\end{array}                                      \tag{6}
\]

All roots displayed with opposite signs are disjoint, and roots carrying
the same unit sign are distinct.

### Two energy-four exceptions

There are \(R-2=k-1\ge4\) ordinary directions. The left side of (4) is a
binary cubic, so its ordinary zeros force it to vanish identically. On the
second energy-four pattern in (6), its value factors as

\[
6(a-d)(b-d)(c-d)\ne0.                            \tag{7}
\]

Thus both exceptions must be doubled points \(2\delta_a\).

Define the binary quadratic

\[
D=2q_2-q_1^2.
\]

It is a nonzero square on every ordinary selected direction and vanishes
at the two distinct exceptional directions. Let \(N\) be the anisotropic
binary norm form defining the selected half

\[
T=\{z\in\mathbb P^1(\mathbb F_p):\eta(N(z))=\varepsilon\},
\qquad |T|=R.
\]

The split nondegenerate quadratic \(D\) has projective character sum zero.
Since its two roots lie in \(T\) and \(\eta(D)=1\) on the rest of \(T\),

\[
{p-3\over2}
=\sum_{z\in T}\eta(D(z))
={\varepsilon\over2}
  \sum_{z\in\mathbb P^1(\mathbb F_p)}\eta(N(z)D(z)).
\]

Consequently

\[
\left|\sum_{\mathbb P^1(\mathbb F_p)}\eta(ND)\right|=p-3. \tag{8}
\]

The quartic \(ND\) is squarefree: \(N\) has two nonrational conjugate roots
and \(D\) has two distinct rational roots. Hence \(Y^2=ND\) is a smooth
genus-one double cover of \(\mathbb P^1\). Hasse gives

\[
\left|\sum_{\mathbb P^1(\mathbb F_p)}\eta(ND)\right|
\le2\sqrt p,
\]

contradicting \(p-3>2\sqrt p\) for \(p\ge11\).

### One energy-six exception

Now \(R-1=k\ge5\) directions are ordinary, so both (4) and (5) are binary
form identities. On \(2\delta_a+\delta_b-\delta_c\), the cubic defect is

\[
6(a-c)^2(b-c)\ne0.                               \tag{9}
\]

For the remaining four-positive/two-negative pattern, choose
\(u,v\) over the algebraic closure with

\[
u+v=q_1,\qquad u^2+v^2=q_2
\]

at the exceptional direction. Equations (4)--(5) give
\(u^d+v^d=q_d\) also for \(d=3,4\). Moving the two negative roots to the
other side produces two four-element multisets with equal first four power
sums. Since \(p>4\), Newton identities make the multisets equal. That is
impossible because the original positive and negative roots are disjoint.
This finishes \(t=2\).

## The \(t=p-1\) branch

Set \(b_{j,s}=1-a_{j,s}\). Each transformed profile has sum one. For
\(1\le d\le4\),

\[
\sum_{s\in\mathbb F_p}s^d=0,
\]

so its moment vectors are again values of degree-\(d\) binary forms, with
an irrelevant overall sign. Equation (3) says that the total transformed
energy is \(R+4\). Ordinary profiles are single deltas. The only
possibilities are:

- two energy-three profiles
  \(\delta_a+\delta_b-\delta_c\); or
- one energy-five profile, either
  \(2\delta_a-\delta_b\) or three positive and two negative unit roots.

For an ordinary delta, \(q_2-q_1^2=0\). In the first case this binary
quadratic has \(R-2=k-1\ge4\) ordinary zeros and is therefore an identity,
but its exceptional defect is

\[
-2(a-c)(b-c)\ne0.                                \tag{10}
\]

In the one-exception case there are \(R-1=k\ge5\) ordinary zeros. The
doubled pattern has defect

\[
-2(a-b)^2\ne0.                                   \tag{11}
\]

For the three-positive/two-negative pattern, the ordinary directions force
\(q_d=q_1^d\) as binary-form identities for \(d=1,2,3\). If
\(u=q_1\) at the exception, moving \(u\) to the negative side gives two
three-element multisets with equal first three power sums. Newton identities
again force equality, contradicting disjoint support.

## The \(t=p+1\) branch

Set \(b_{j,s}=a_{j,s}-1\). These profiles also have common sum one, but
(3) leaves only two excess units. Thus exactly one profile has energy
three and is \(\delta_a+\delta_b-\delta_c\); all other \(R-1=k\ge5\)
profiles are deltas. The binary quadratic \(q_2-q_1^2\) is forced to be an
identity and (10) gives the final contradiction.

Combining these three exclusions with Proposition 15.637 proves the boxed
empty-shell statement.

## Independent audits and search record

The certificate module exhausts the balancing values and checks every
defect identity. Its test suite also exhausts all integer coefficient
multisets of sums one and two through energy six. A separate finite-field
enumeration found no selected-half binary quadratic with the forbidden
two-root/nonnegative pattern for every tested prime \(11\le p\le43\).
These are audits; equations (1)--(11) are the proof.

Targeted searches found the general conference-ETF lattice literature and
standard Hasse/character-sum applications, but no Paley-dual theorem
excluding this shell. The candidate scaled values
\(28,32,40,44,52,64,68,\ldots\) collide with unrelated OEIS sequences
(for example A354810 and A092259); they are merely \(2(p+3)\), and no new
integer-sequence claim is made.

Evidence:

- src/e1_gmin_m4_prop15638.py;
- evidence/e1_gmin_m4_prop15638.json;
- scripts/r1_next_shell_half_conic_audit.py;
- evidence/r1_next_shell_half_conic_11_43.json;
- tests/test_prop15638.py.
