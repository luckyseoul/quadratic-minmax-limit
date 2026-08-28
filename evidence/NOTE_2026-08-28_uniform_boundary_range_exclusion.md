# Uniform non-Walsh boundary-range exclusion

Date: 2026-08-28. This is Proposition 15.669. It extends the exact
hypergeometric parity floors of Propositions 15.652 and 15.657 from at most
six odd fibres to the full middle range. Combining those floors with the
type-split affine budget and the finite-plane pair-deficit bound gives a
uniform exclusion range for every odd prime at least 17, plus new exact
cases at 11 and 13.

This is not a closure of residual (ii). The first count profiles beyond the
stated ranges survive the floor-and-pair relaxation; no actual residual
graph is constructed by such a profile.

## 1. The full middle parity floor

Let \(X\) be uniform on the middle slice
\(J(p,(p+1)/2)\), let \(B\subseteq\mathbb F_p\) have size \(b\), and put
\(t=|X\cap B|\). Its first two moments are

\[
 \mu=\mathbb E t={b(p+1)\over2p},\qquad
 \mathbb E t^2={b(b+1)(p+1)\over4p},
\]

\[
 \sigma^2=\operatorname{Var}(t)
 ={(p+1)b(p-b)\over4p^2}.                         \tag{15.669.1}
\]

For phase \(\eta\), a feasible averaged quadratic majorant obeys

\[
 q(t)\ge (t+\eta\bmod2).
\]

The constant polynomial \(q=1\) is feasible. To certify it as optimal, it
suffices to put a positive measure with moments (15.669.1) on the contact
nodes \(t+\eta\) odd.

Complementing \(B\) sends \(t\) to \((p+1)/2-t\), so it suffices to take

\[
 5\le b\le{p-1\over2}.
\]

Let \(j=1-\eta\), let \(r=j\), and let \(R\) be the largest integer at most
\(b\) congruent to \(j\) modulo two. At the prescribed mean, the largest
variance of a measure on the contact nodes is

\[
 (\mu-r)(R-\mu).
\]

The smallest variance is
\((\mu-a)(a+2-\mu)\), where \(a,a+2\) are the adjacent contact nodes
bracketing \(\mu\); it is zero when \(\mu\) itself is a contact.

After multiplication by \(4p\), the four possible upper-envelope margins
factor as

\[
\begin{array}{c|c}
(j,R)&4p\bigl((\mu-r)(R-\mu)-\sigma^2\bigr)\\ \hline
(0,b)&b(b-1)(p+1)\\
(1,b)&b\{b(p+1)-3p+1\}\\
(0,b-1)&b(b-3)(p+1)\\
(1,b-1)&(b-1)\{b(p+1)-4p\}.
\end{array}                                           \tag{15.669.2}
\]

Depending on \(b\bmod4\) and \(j\), the four lower-envelope margins are

\[
 {b(p-b-3)\over4p},\quad
 {b(p+1-b)-4p\over4p},\quad
 {b(p-1-b)-3p\over4p},\quad
 {(p-b)(b-3)\over4p}.                                \tag{15.669.3}
\]

All expressions in (15.669.2) are positive for \(b\ge5\). In
(15.669.3), the first and fourth are immediate. The two middle numerators
are concave in \(b\). On their relevant intervals their endpoint lower
bounds are

\[
 2p-30,\qquad {p^2-14p-3\over4},\qquad
 {p^2-14p+1\over4},
\]

which are positive for \(p\ge17\).

Thus \(\sigma^2\) lies between the two envelope variances. Mix the
two-point lower-envelope measure with the two-point endpoint measure in the
unique rational proportion giving \(\sigma^2\). This is an explicit
positive degree-two quadrature on parity-one contacts. It has the same
moments through degree two as the hypergeometric law, so every feasible
quadratic has expectation at least one. Therefore

\[
 \boxed{M(p,b,\eta)=1,\qquad
        f_\eta(b)=2p\quad(5\le b\le p-5,\ p\ge17).}    \tag{15.669.4}
\]

The source constructs the rational weights, transports them back across
complementation, and verifies positivity, support, contact parity, and all
three moments.

## 2. No infinity in the boundary

Let the all-finite odd-degree boundary have even size \(s\). There are
\(m=(p+1)/2\) directions of each quadratic type; their phases are opposite.
Each type has budget \(m(p+1)\). Starting from middle cost \(2p\) in every
direction, each type must save

\[
 R=m(p-1)={p^2-1\over2}.                            \tag{15.669.5}
\]

For phase zero the only savings are

\[
\begin{array}{c|ccc}
b&0&2&4\\ \hline
\text{saving}&2p&p-1&6\\
\text{deficit }s-b&s&s-2&s-4.
\end{array}
\]

For \(p\ge17,s\ge6\), \(b=0\) has the best saving per deficit, and
\(b=2\) is best among the remaining choices:

\[
 s(p+1)>4p,\qquad s(p-7)>4(p-4).                    \tag{15.669.6}
\]

Put \(a=\lfloor R/(2p)\rfloor\) and \(r=R-2pa\). If at least \(a+1\)
zero-fibre directions occur, their deficit is at least \((a+1)s\). If at
most \(a\) occur, (15.669.6) gives deficit at least

\[
 as+{r(s-2)\over p-1}.
\]

For \(p\equiv1\pmod4\), \(a=(p-1)/4,r=(p-1)/2\). For
\(p\equiv3\pmod4\), \(a=(p-3)/4,r=(3p-1)/2\). Both branches therefore give

\[
 D_0\ge{(p+1)s\over4}-1.                            \tag{15.669.7}
\]

In phase one, only \(b=2\) saves anything, namely \(p+1\). Hence at least
\((p-1)/2\) such directions are required and

\[
 D_1\ge{(p-1)(s-2)\over2}.                          \tag{15.669.8}
\]

Every pair of finite boundary points collides in one projective direction,
so Proposition 15.657 gives

\[
 \sum_d(s-b_d)\le s(s-1).                           \tag{15.669.9}
\]

Subtracting (15.669.9) from (15.669.7)--(15.669.8), four times the
contradiction gap is

\[
 h_p(s)=s(3p+3-4s)-4p.                              \tag{15.669.10}
\]

This is concave in \(s\). At the endpoints of
\([6,3(p-1)/4]\), its values are

\[
 h_p(6)=14p-126>0,\qquad
 h_p(3(p-1)/4)={p-9\over2}>0.
\]

Consequently

\[
 \boxed{\text{every all-finite even boundary with }
 6\le s\le {3(p-1)\over4}\text{ is impossible for }p\ge17.} \tag{15.669.11}
\]

## 3. Infinity in the boundary

Now \(s\) is the odd number of finite boundary points and the phase is the
same in both quadratic types.

In phase zero, \(b=1\) saves \(p-1\), while every other allowed value of
\(b\) saves at most six. Meeting (15.669.5) therefore forces every direction
of both types to have \(b=1\), requiring total deficit

\[
 (p+1)(s-1).
\]

In phase one, \(b=1\) saves \(p+1\). Every other allowed value saves at most
six. The only extra endpoint possibility is \(b=p-4\) when
\(p\equiv1\pmod4\), obtained by complementing the phase-zero \(b=4\)
floor. If a type had at most \(m-2\) directions with \(b=1\), its total
saving would be at most

\[
 (m-2)(p+1)+12,
\]

which falls short of \(m(p-1)\) by \(p-11>0\). Thus each type needs at
least \(m-1=(p-1)/2\) directions with \(b=1\), requiring total deficit

\[
 (p-1)(s-1).
\]

For odd \(5\le s\le p-4\), both quantities strictly exceed the geometric
budget \(s(s-1)\). Hence

\[
 \boxed{\text{every infinity-present boundary with }
 5\le s\le p-4\text{ finite points is impossible for }p\ge17.} \tag{15.669.12}
\]

## 4. Exact extensions at p=11 and p=13

For the two primes below 17, the source evaluates the original rational
parity-majorant LP for every allowed \(b\), then performs an exact dynamic
program over the \((p+1)/2\) directions of one type. The minimum required
pair-deficit gaps are:

| prime and boundary | phase-zero gap | phase-one gap |
|---|---:|---:|
| \(p=11\), infinity plus 7 | 30 | 18 |
| \(p=13\), infinity plus 7 | 42 | 30 |
| \(p=13\), infinity plus 9 | 40 | 24 |

For \(p=13\), eight finite points have opposite-phase total gap \(4\).
All are positive, so these cases are excluded exactly. The next
floor-plus-pair survivors are eight finite points at \(p=11\), ten finite
points at \(p=13\), and infinity plus 9 or 11 respectively.

## 5. Sharp scope of this relaxation

The exact count-profile dynamic program gives the first all-finite
floor-plus-pair survivors

| \(p\) | 17 | 19 | 23 | 29 | 31 | 37 | 47 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| first surviving even \(s\) | 14 | 14 | 18 | 22 | 24 | 28 | 36 |

These are exactly the first even integers strictly above
\(3(p-1)/4\). With infinity, the first survivor has \(s=p-2\) finite
points. A survivor here is only a multiset of directional odd-fibre counts
compatible with the scalar floors and pair budget. It is not an affine
boundary, a slack catalog, or an edge set. Closing those profiles requires
new incidence information, higher moments, or exact finite models.

## 6. Literature and OEIS check

Targeted searches for parity-constrained integer-valued quadratics on
Johnson slices, hypergeometric quadratic majorants, and odd-fibre incidence
bounds in affine planes found only adjacent Johnson association-scheme and
finite-incidence literature. In particular, Filmus--Kindler--Mossel--Wimmer
concerns low-degree analysis on slices, Blekherman--Gouveia--Pfeiffer
concerns nonnegative polynomials on the cube, and Bernard--Crampé--Vinet
[arXiv:2306.01882](https://arxiv.org/abs/2306.01882) concerns the
non-binary Johnson scheme. None states (15.669.4) or combines it with the
Paley type budget and (15.669.9). This is a duplicate/context check, not an
unqualified priority claim.

Individual OEIS searches were made for the larger endpoint values 758,
7,506, 75,002, and 750,002. The last returned no entry; the others occur in
unrelated partition, recurrence, grid, and trajectory sequences. The
endpoint formula is elementary and no sequence submission or novelty claim
is proposed.

## 7. Reproduction

~~~bash
python src/e1_gmin_m4_prop15669.py
python -m pytest -q tests/test_prop15669.py
~~~

The generated exact record is
evidence/e1_gmin_m4_prop15669.json.
