# The complete shell at the first nonminimal odd scaled norm

Date: 2026-08-25. Proposition 15.639. Let

\[
n=p^2+1,\qquad P={1\over2}(I+C/p),\qquad L^*=P\mathbb Z^n.
\]

For every odd prime \(p\ge11\), Proposition 15.635 gives the first possible
nonminimal odd scaled norm

\[
s_{\rm odd}=2p\lVert x\rVert^2=3p-6.             \tag{1}
\]

This note classifies the complete shell at (1). It is the disjoint union of
two explicit signed families.

First, let \(\mathcal T\) consist of \(Pz\), where \(z\) has three
coordinates in \(\{\pm1\}\), all others zero, and

\[
z_i z_j C_{ij}=-1
\]

on all three support edges. Second, for each square
\(\mathbb F_p\)-subline \(S\), orient its signed complement \(w_S\) so
that \(Cw_S=pw_S\), and let

\[
\mathcal O=
\left\{\mathord\pm\left(Pe_i+{w_S\over p}\right):
             (w_S)_i=-1\right\}.                 \tag{2}
\]

The orientation in (2) is unique for each point outside \(S\). The theorem
is

\[
\boxed{\{x\in L^*:2p\lVert x\rVert^2=3p-6\}
       =\mathcal T\mathbin{\dot\cup}\mathcal O.}  \tag{3}
\]

Its signed size is

\[
\boxed{N_{\rm odd}(p)
=\binom{p^2+1}{3}+p^2(p-1)(p^2+1)
={p^2(p-1)(p+7)(p^2+1)\over6}.}                  \tag{4}
\]

For \(p=11,13\), this is the fourth nonempty shell. For \(p\ge17\), further
even candidates can lie between the third shell and (1), so (3) is not
asserted to be the fourth shell globally.

## Every vector at scaled norm 3p-6 has a unit scaled coordinate

Write \(x=Pz\), \(z\in\mathbb Z^n\), and

\[
r=2px=(pI+C)z\in\mathbb Z^n.
\]

Modulo two, every entry of \(pI+C\) is one. Also

\[
s=p\lVert z\rVert^2+z^TCz\equiv\sum_i z_i\pmod2.
\]

Thus \(s\) odd implies that every \(r_i\) is odd. On the shell (1),

\[
\sum_i r_i^2=2ps=6p^2-12p<9(p^2+1).              \tag{5}
\]

If every \(|r_i|\ge3\), the left side of (5) is at least \(9n\), a
contradiction. Hence some \(|r_i|=1\). A signed Paley automorphism is an
integral monomial matrix \(M\) with \(MCM^T=C\); it preserves
\(P\mathbb Z^n\), both families in (3), and the scaled norm. Signed
\(\operatorname{PSL}(2,p^2)\) is transitive on the coordinates. Thus move
that coordinate to the profile base point and negate \(x\) if necessary.
Since the common profile sum at that base point is the corresponding
coordinate of \(r\), it is now

\[
t=1.                                              \tag{6}
\]

## Equality in the \(t=1\) profile bound

Let \(R=(p+1)/2\), let \(m=Pe_u\) be the minimum vector determined by the
degree-one profile moment, and put \(y=x-m\). Its profiles

\[
b_j=a_j-\delta_{\mu_j}
\]

have sum zero. Suppose exactly \(h\) are active and let \(M\) be their
total positive mass. With

\[
\Delta={s-p\over4}=R-2,
\]

the proof of Proposition 15.635 gives

\[
\Delta\ge h,\qquad
\Delta\ge M-h,\qquad
M\ge h(R-h).                                     \tag{7}
\]

For \(1\le h<R\),

\[
h(R-h-1)-(R-2)=(h-1)(R-h-2).
\]

Combining this with (7) leaves only

\[
\boxed{h=1\quad\hbox{or}\quad h=R-2.}             \tag{8}
\]

## The one-active-profile case

Here \(M=R-1=(p-1)/2\), and equality holds in every local energy bound. If
the distinguished entry is \(b_{\mu}=-c\), equality says

\[
\sum_s |b_s|(|b_s|-1)=2(c-1).
\]

The left side contains \(c(c-1)\), so \(c\in\{1,2\}\). For \(c=2\), all
other nonzero entries are units and \(p\lVert y\rVert^2=p+1\). The positive
and negative degree-
\((p-1)/2\) root multisets then have equal power sums through one degree
less, cover all but two field elements, and one side has a repeated root.
This is exactly the profile excluded uniformly by Proposition 15.636.

Thus \(c=1\), every nonzero entry is a unit, and the one-profile equality
classification from Proposition 15.633 makes

\[
y={w_S\over p}
\]

for an oriented square-circle complement. Since
\(\langle Pe_u,y\rangle=y_u\), the norm (1) forces \(y_u=-1/p\).
Therefore \(x\) belongs to \(\mathcal O\).

## The \(R-2\)-active-profile case

Equality in (7) gives positive mass two in every active profile and local
excess one. An inactive original profile is \(\delta_\mu\); an active one
has three distinct support points and the form

\[
a_j=\delta_\alpha+\delta_\beta-\delta_\gamma,
\qquad
\alpha+\beta-\gamma=\mu.                         \tag{9}
\]

Put \(A=\alpha-\mu\), \(B=\beta-\mu\), so
\(\gamma=\mu+A+B\). If \(q_d\) is the degree-\(d\) moment form, define

\[
Q_2=q_2-\mu^2,\qquad Q_3=q_3-\mu^3.
\]

Here \(A\) and \(B\) are nonzero on every active direction: otherwise two
of \(\alpha,\beta,\gamma\) would coincide. On active directions,

\[
Q_2=-2AB,\qquad
Q_3=-3AB(2\mu+A+B).                              \tag{10}
\]

Both forms vanish at the two inactive directions, while \(Q_2\) is nonzero
at every active direction. The binary quadratic \(Q_2\) therefore has
exactly those two distinct projective roots, hence divides the binary cubic
\(Q_3\). Consequently

\[
S={2\over3}{Q_3\over Q_2}-2\mu                  \tag{11}
\]

is a binary linear form whose active values are \(A+B\). The linear form
\(\nu=\mu+S\) is the line-coordinate profile of one point \(v\). Set

\[
Y=x+Pe_v.
\]

On every active direction, the profiles of \(Y\) are
\(\delta_\alpha+\delta_\beta\). On an inactive direction they are

\[
\delta_\mu+\delta_{\mu+S}.
\]

Let \(r\) be the number of the two inactive directions at which \(S=0\).
Since \(S\) is linear, \(r\in\{0,1,2\}\), and \(Y\) has exactly \(r\)
doubled profiles. The circle-frame identity gives

\[
2p\lVert Y\rVert^2=2(p-1)+4r.                   \tag{12}
\]

For \(r=2\), (12) is the empty shell \(2(p+3)\) excluded by Proposition
15.638. For \(r=1\), it is the complete third shell \(2(p+1)\). Proposition
15.636 writes

\[
Y=\mathord\pm P(e_i+C_{ij}e_j).
\]

Subtracting \(Pe_v\) either cancels one support point, doubles it, or leaves
three points with one signed conference edge positive. Directly in
\(2p\lVert Pz\rVert^2=p\lVert z\rVert^2+z^TCz\), these cases have scaled
norm \(p\), at least \(3p-2\), or at least \(5p+4\), never \(3p-6\).
Thus \(r=1\) is impossible.

It follows that \(r=0\), so \(Y\) lies on the complete second shell from
Proposition 15.633. If \(Y\) is a projected point-pair vector, overlap of
\(v\) with its two support points gives scaled norm \(p\) or \(5p-4\), so
it cannot give (1). Otherwise \(x=Y-Pe_v\) has three unit coordinates. Its
norm is (1) exactly when all three signed conference edges are negative,
giving \(\mathcal T\). If
\(Y=w_S/p\), then

\[
\lVert Y-Pe_v\rVert^2={p-1\over p}+{1\over2}-2Y_v
\]

equals (1) exactly when \(Y_v=1/p\). Up to global sign this is precisely
\(\mathcal O\). Undoing the initial signed Paley automorphism preserves both
descriptions. This proves exhaustion in (3).

## Converse, disjointness, and count

For a negative signed triple \(z\),

\[
2p\lVert Pz\rVert^2=3p-6.
\]

For a point-circle vector in (2),

\[
\left\lVert Pe_i+{w_S\over p}\right\rVert^2
={1\over2}+{p-1\over p}+{2(w_S)_i\over p}
={3p-6\over2p}.
\]

The two families are visibly distinct in scaled coordinates. A negative
triple has exactly three coordinates of magnitude \(p-2\); every other
coordinate has magnitude at most three. A point-circle vector has exactly
one coordinate of magnitude \(p-2\), and all others have magnitude at most
three. Since \(p\ge11\), these signatures also make both parametrizations
injective.

Finally,

\[
\operatorname{tr}(C^3)=0
=6\sum_{\{i,j,k\}}C_{ij}C_{ik}C_{jk},
\]

because \(C\) has equally many eigenvalues \(p\) and \(-p\). Hence exactly
half of all coordinate triples are negative, and each has one antipodal
sign pair. This gives \(|\mathcal T|=\binom n3\).

There are \(p(p^2+1)/2\) unoriented square circles and
\(p^2-p\) points off each. For every off-circle point exactly one
orientation has value \(-1\), and global signs double the result. Therefore

\[
|\mathcal O|=p^2(p-1)(p^2+1).
\]

This proves (4).

## Boundary and search record

An independent exact PARI/GP `qfminim` enumeration on NUKA used the
saturated dual Gram form through scaled bound \(28\) at \(p=11\). It
returned

```
P=11
BOUND=28
SIGNED_COUNT=473970
MAXNORM=27
ELAPSED_MS=2033141
```

The cumulative count through the complete third shell is \(31,110\). The
exact residual is therefore \(473,970-31,110=442,860\), exactly (4), and
there is no vector at the excluded bound \(28\). This computation is an
audit, not an input to the uniform proof.

The theorem classifies the full scaled-norm \(3p-6\) shell. Proposition
15.640 subsequently computes its degree-four harmonic operator. Neither
result excludes the additional even candidates below \(3p-6\) when
\(p\ge17\); those and the later tail are real boundaries before any global
R1 claim.

Targeted literature and OEIS searches found no Paley-dual shell theorem or
entry for the count formula (4). The bare negative-triangle half-count is
an immediate conference-matrix trace identity; the content here is the
profile exhaustion and the second-shell reduction. No broad novelty claim
is made without external review.

Evidence:

- src/e1_gmin_m4_prop15639.py;
- evidence/e1_gmin_m4_prop15639.json;
- evidence/r1_dual_shell_count_p11_28.json;
- tests/test_prop15639.py.
