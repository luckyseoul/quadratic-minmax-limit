# Two-half lower endpoint: uniform near-pencil reduction

**Date:** 2026-09-02
**Status:** proved reduction; residual (ii) remains open

This note treats only the lower endpoint of the Proposition 15.758 branch B.
It is not a small-prime census and it does not assert that the surviving
near-pencil graphs satisfy the global undercutter premise.

Put

\[
 p=4r+1,\qquad m=2r+1,\qquad t=2r^2-5r,
 \qquad |H|=4p+2t+1={p^2-1\over4}+p+4.
\]

There are \(m\) Paley-hard directions. In each hard direction let
\(\ell_L\) be the exceptional affine fibre, and let

\[
 \nu_x=\#\{L:x\in\ell_L\},\qquad
 K=\sum_x\left\lfloor{(\nu_x-1)^2\over4}\right\rfloor .
\tag{1}
\]

The directional boundary parities give

\[
 |D|=p^2-m^2-4K.                                    \tag{2}
\]

Since every graph with boundary \(D\) has \(2|H|\ge |D|\), (2) implies

\[
 r^2-2r-2\le K\le r^2.                              \tag{3}
\]

For completeness, the upper bound is also exact. The ratio

\[
 {\lfloor(s-1)^2/4\rfloor\over {s\choose2}}
\]

is nondecreasing in \(s\). Charging every concurrence to its point-pairs
therefore gives

\[
 K\le {\lfloor(m-1)^2/4\rfloor\over {m\choose2}}{m\choose2}=r^2.
\]

A full pencil attains equality. The point of this note is the following
stability statement.

## Theorem

For every integer \(r\ge13\), (3) forces at least \(m-2\) of the exceptional
lines \(\ell_L\) to be concurrent.

### Proof

Dualize the \(m\) exceptional lines to a set \(S\) of \(m\) points. For
each line determined by \(S\), write \(s_\ell=|S\cap\ell|\), and put

\[
 W=\sum_\ell\left\lfloor{s_\ell\over2}\right\rfloor .
\]

For every integer \(s\ge2\),

\[
 {s\choose2}-2\left\lfloor{(s-1)^2\over4}\right\rfloor
 =\left\lfloor{s\over2}\right\rfloor .             \tag{4}
\]

Every pair of points of \(S\) determines one line, so (3)--(4) give

\[
 W={m\choose2}-2K=r+2(r^2-K)\le5r+4.                \tag{5}
\]

Let \(k=\max_\ell s_\ell\), let \(d=m-k\), and fix a \(k\)-secant \(B\).
The connector lines through each point of \(B\) partition the \(d\) outside
points. Therefore

\[
 W\ge\left\lfloor{k\over2}\right\rfloor
       +k\left\lceil{d\over2}\right\rceil .         \tag{6}
\]

There is also an exact pair-efficiency bound. Namely

\[
 {s\choose2}=\left\lfloor{s\over2}\right\rfloor e_s,
 \qquad
 e_s=\begin{cases}s,&s\text{ odd},\\s-1,&s\text{ even}.
 \end{cases}
\]

If \(q_k\) is the largest odd integer at most \(k\), then \(e_s\le q_k\)
for \(s\le k\), and hence

\[
 {m\choose2}\le q_k W.                              \tag{7}
\]

Suppose first that \(k\le r\). Then \(d\ge r+1\), and (5)--(6) imply

\[
 k\le {2(5r+4)\over r+1}=10-{2\over r+1}<10.
\]

For \(k\le6\), (7) gives

\[
 r(2r+1)\le5(5r+4),
\]

contrary to \(2r^2-24r-20>0\) for \(r\ge13\). For \(k=7,8,9\), the
right side of (6) is respectively

\[
 7r-18,\qquad8r-20,\qquad9r-32,
\]

and each exceeds \(5r+4\). Thus \(k\ge r+1\), so \(d\le r\).

If \(7\le d\le r\), (6) and (5) give

\[
 {d(2r+1-d)\over2}\le5r+4.                         \tag{8}
\]

The left side is concave in \(d\). At the two endpoints it equals
\(7r-21\) and \(r(r+1)/2\), both larger than \(5r+4\) for \(r\ge13\).
For \(d=5,6\), the exact bound (6) is \(7r-14,7r-18\), again too large.
Consequently \(d\le4\).

It remains to exclude \(d=3,4\). Every determined line other than \(B\)
contains at most one point of \(B\). Charging its contribution to its
pairs of outside points gives

\[
 K\le\left\lfloor{(k-1)^2\over4}\right\rfloor+{d\choose2}. \tag{9}
\]

For \(d=3\), (9) gives \(r^2-K\ge3r-5\); for \(d=4\), it gives
\(r^2-K\ge4r-10\). Both contradict (3), which says
\(r^2-K\le2r+2\). Hence \(d\le2\), as claimed. \(\square\)

## Exact surviving geometries

Let \(h_0=(p^2-1)/4\), let \(a\) be the common pencil point, and let

\[
 D_0=\mathbf F_p^2\setminus
 \bigcup_{L\ {\rm hard}}\{x:L(x)=L(a)\}.
\]

This is the Paley-hemisphere boundary, of size \(2h_0\). If an outlier
\(\ell_i\) replaces the pencil fibre \(\ell_i^0\) of the same direction,
then the three possible cases are:

| outliers | boundary | \(K\) | \(s=|H|-|D|/2\) | core-edge cap |
|---:|---|---:|---:|---:|
| 0 | \(D_0\) | \(r^2\) | \(p+4\) | \(3p+12\) |
| 1 | \(D_0\mathbin\triangle\ell_1^0\mathbin\triangle\ell_1\) | \(r^2-r\) | \(2r+5\) | \(6r+15\) |
| 2, ordinary | \(D_0\mathbin\triangle\ell_1^0\mathbin\triangle\ell_1\mathbin\triangle\ell_2^0\mathbin\triangle\ell_2\) | \((r-1)^2\) | \(7\) | \(21\) |
| 2, triple | same | \((r-1)^2+1\) | \(9\) | \(27\) |

The last distinction records whether \(\ell_1\cap\ell_2\) lies on a
retained pencil line.

The core bounds use the following graph lemma. If a graph has \(d\)
odd-degree vertices and \(d/2+s\) edges, delete all isolated \(K_2\)
components. In every remaining connected component with \(e\) edges and
\(b\) odd vertices, \(b\le4e/3\): this is immediate for a cyclic component,
for a tree with at least three edges it follows from \(b\le e+1\le4e/3\),
and the two-edge path is checked directly. Thus

\[
 e-b/2\ge e/3,
\]

and the remaining core has at most \(3s\) edges. In particular, the
two-outlier boundary is a huge isolated matching plus at most 21 or 27
core edges.

## Exact global-cut consequence

The global undercutter premise supplies an additional condition which is
not contained in the incidence data above. Let \(C\) be the Paley
conference matrix of order \(n=p^2+1\), let
\(\Phi=np/2\), and put

\[
 A=C\mathbin\oplus G,\qquad H=G\cup\{e\},\qquad
 B=C\mathbin\oplus H.
\]

Assume

\[
 \Phi(A)=\Phi-2,\qquad \Phi(B)\le\Phi-4,
\]

and choose \(x\in\{\pm1\}^n\) and \(\epsilon\in\{\pm1\}\) with
\(\epsilon q_A(x)=\Phi-2\). Define

\[
 \delta=\Phi-\epsilon q_C(x)=2c.
\]

The one-edge difference first forces

\[
 \epsilon f_e(x)=1,\qquad \epsilon S_H(x)=2-c.       \tag{10}
\]

Every Boolean \(\epsilon p\)-eigenvector \(y\) satisfies
\(\epsilon S_H(y)\ge3\): the \(B\)-bound gives at least two, and
\(|H|\) is odd. If \(x\) is obtained from \(y\) by switching a coordinate
set \(T\), then

\[
 2\epsilon\,c_H(T;y)
 =\epsilon S_H(y)-\epsilon S_H(x)
 =c+\epsilon S_H(y)-2\ge c+1,                       \tag{11}
\]

where \(c_H(T;y)\) is the \(C_{uv}y_uy_v\)-signed \(H\)-cut of \(T\).
Thus every sign-compatible eigenshell representation obeys

\[
 |H\cap\delta(T)|\ge {c+1\over2}.                   \tag{12}
\]

If \(H=F\mathbin{\dot\cup}R\), where \(F\) is the union of its isolated
\(K_2\) components and \(|R|=\rho\), and
\(u=\min(|T|,n-|T|)\), then

\[
 |H\cap\delta(T)|\le u+\rho.
\]

Consequently

\[
 \boxed{c\le2u+2\rho-1,\qquad \rho\le3s.}           \tag{13}
\]

Equations (10)--(13) are a genuine use of the full Boolean maximum. They
also identify the remaining issue precisely: the lattice-shell
classification gives a projected vector, but it does not automatically
give a sign-compatible small switching set \(T\). Such compatibility must
be proved rather than assumed.

## A sharper cut-surplus bound

The preceding \(u+\rho\) estimate discards the fact that vertices used by
the core are then unavailable to the isolated matching.  The following
bound keeps exactly that information.

Let a graph \(J\) have \(h\) edges, \(d\) odd-degree vertices, and surplus

\[
 s_J=h-d/2.
\]

For every cut with smaller shore of size \(u\),

\[
 \boxed{|J\cap\delta(T)|\le u+2s_J.}               \tag{14}
\]

Indeed, let \(V_+(J)\) be the vertices of positive degree and replace
\(T\) by its smaller shore.  Then

\[
\begin{aligned}
 |J\cap\delta(T)|
 &\le \sum_{v\in T\cap V_+(J)}\deg_J(v)\\
 &\le u+\sum_{v\in V_+(J)}(\deg_J(v)-1)\\
 &=u+2h-|V_+(J)|\\
 &\le u+2h-d=u+2s_J,
\end{aligned}
\]

because every odd-degree vertex belongs to \(V_+(J)\).  This proof is
independent of a chosen matching/core decomposition.  Applied to (12), it
strengthens (13) to

\[
 \boxed{c\le2u+4s_J-1.}                            \tag{15}
\]

There is also an exact way to state all of the full-Max information.  Put

\[
 K=\epsilon\,\operatorname{diag}(x)C\operatorname{diag}(x),
\]

and, for a coordinate cut \(T\), write

\[
 W_K(T)=\sum_{ab\in\delta(T)}K_{ab},\qquad
 Z_H(T)=\sum_{ab\in H\cap\delta(T)}K_{ab}.
\]

In the \(x,\epsilon\) gauge the \(B\)-cut is
\(W_B(T)=W_K(T)-2Z_H(T)\).  Since
\(\epsilon q_B(x)=\Phi-4\), the assertion that \(x\) is an absolute
maximizer of \(B\) is equivalent to the complete cut box

\[
 0\le W_K(T)-2Z_H(T)\le\Phi-4
 \quad\hbox{for every }T,                           \tag{16}
\]

or, equivalently,

\[
 {W_K(T)-\Phi+4\over2}\le Z_H(T)\le {W_K(T)\over2}.
                                                               \tag{17}
\]

These inequalities already imply the corresponding \(A\)-cut box: the
distinguished edge has \(K_e=1\), so
\(W_A(T)=W_B(T)+2\mathbf1_{e\in\delta(T)}\) and the \(A\) endpoint is
two larger.  Thus (16), rather than two independent hereditary cut
families, is the exact shared-maximizer condition.

## The zero-circle-mismatch triple is excluded

Consider the all-bad signed-triple endpoint.  In the same gauge there is a
three-set \(E\) such that

\[
 K_{ij}=1\quad(i,j\in E),\qquad
 v=\mathbf1+2\mathbf1_E,\qquad Kv=pv.              \tag{18}
\]

Let \(\Gamma(E)\) be the unique Miquelian \(\mathbf F_p\)-subline through
\(E\).  Its type is forced by the positive signed triangle in (18), so its
sparse circle vector can be oriented as

\[
 \operatorname{supp}(d)=\Gamma(E),\qquad Kd=pd,qquad
 d_i=1\quad(i\in E).                               \tag{19}
\]

To check the sign, a sparse circle eigenvector has one constant switched
edge sign \(\tau\) on its whole support and eigenvalue \(\tau p\).  The
product of the three switched edges on any triple in that circle is
\(\tau^3=\tau\).  Equation (18) makes this product positive, hence
\(\tau=1\), which proves both assertions in (19).

Define the circle-sign mismatch

\[
 \mu(x,E)=\#\{j\in\Gamma(E)\setminus E:d_j=-1\}.   \tag{20}
\]

Before passing to the \(x\)-gauge, this says to orient the sparse circle
vector to agree with \(x\) on \(E\), and count its disagreements with
\(x\) on the rest of the circle.  Hence \(\mu\) is invariant under every
signed projective automorphism of \(C\): the Boolean vector and the circle
vector acquire the same diagonal switching signs, and their coordinatewise
products are merely permuted.  Projective transitivity on \(E\) therefore
does not determine \(\mu\).

If \(\mu=0\), then

\[
 h=v-2d
\]

is a Boolean \(+p\)-eigenvector of \(K\).  It differs from \(\mathbf1\)
exactly on

\[
 T=\Gamma(E)\setminus E,qquad |T|=p-2.            \tag{21}
\]

Here \(c=3p-6\), so (12), (14), and (21) give

\[
 {3p-5\over2}
 \le |H\cap\delta(T)|
 \le p-2+2s.                                      \tag{22}
\]

In the ordinary two-outlier case \(s=7\), whose upper bound in (22) is
\(p+12\); this contradicts (22) for \(p>29\).  In the triple case
\(s=9\), the upper bound is \(p+16\), a contradiction for \(p>37\).
Consequently both zero-mismatch cases are excluded throughout the present
range \(p\ge53\).

The Kiss--Somlai construction in
`NOTE_2026-09-02_KISS_SOMLAI_TRIANGULAR_THREE_SPIKE.md` has
\(\mu=0\) explicitly.  Its circle is the augmented affine line
\(M\cup\{\infty\}\), and the circle signing agrees with its Boolean shadow
at all \(p-2\) nonspike points.  Equivalently, subtracting the hard line
from

\[
 f=\mathbf1_S+\mathbf1_L+\mathbf1_M
\]

leaves the disjoint \(0/1\) function
\(g=\mathbf1_S+\mathbf1_L\).  The line identity
\(Q\mathbf1_M=p\mathbf1_M-\mathbf1\) and
\(Qf=pf-(p+3)\mathbf1/2\) give

\[
 Qg=pg-{p+1\over2}\mathbf1,
\]

so \(w_\infty=1\), \(w_z=2g(z)-1\) is the Boolean eigenvector in
(21).  Since \(M\) meets \(S\cup L\) in exactly its two spike points, the
switch really has size \(p-2\).  Thus the explicit triangular family cannot
be the dangerous shared maximizer in either surviving near-pencil geometry.

This is an exclusion of the zero-mismatch subclass, not a classification of
the all-bad shell.  If \(\mu>0\), then \(v-2d\) has new magnitude-three
coordinates at the mismatches and is not Boolean; neither (21) nor (22)
follows.  Proposition 15.639 classifies the projected signed-triple shell,
but it does not identify all Boolean representatives in its affine
eigenspace fibre.  A separate argument is still required to exclude those
positive-mismatch completions.

There is an exact mismatch ladder, but it does not force \(\mu=0\).  Put

\[
 P_\Gamma=\{j\in\Gamma(E)\setminus E:d_j=1\},\qquad
 N_\Gamma=\{j\in\Gamma(E)\setminus E:d_j=-1\},
\]

so \(|P_\Gamma|=a=p-2-\mu\) and \(|N_\Gamma|=\mu\).  The integral
eigenvector

\[
 z=v-2d
\]

has Boolean shadow \(h=\mathbf1-2\mathbf1_{P_\Gamma}\) and exactly
\(\mu\) entries equal to three, on \(N_\Gamma\).  Since
\(z=h+2\mathbf1_{N_\Gamma}\) and the switched circle is a positive
clique on \(N_\Gamma\),

\[
 \delta_K(h):=\Phi-q_K(h)
 =2p\mu-2\mu(\mu-1)
 =2\mu(p+1-\mu).                                  \tag{23}
\]

On the other hand \(h\) is obtained from \(\mathbf1\) by the cut
\(P_\Gamma\), while \(\delta_K(\mathbf1)=6p-12\).  Therefore

\[
 W_K(P_\Gamma)
 ={\delta_K(h)-(6p-12)\over2}
 =a(\mu-3).                                       \tag{24}
\]

The lower half of the exact Max box (16) consequently says only

\[
 Z_H(P_\Gamma)\le {a(\mu-3)\over2}.               \tag{25}
\]

For \(\mu=1\), this forces at least \(p-3\) crossing edges, while (14)
allows \(p-3+2s\).  For \(\mu=2\), it forces only \((p-4)/2\), and for
\(\mu\ge3\) the right side of (25) is nonnegative.  Thus the exact
circle iteration and the odd-degree surplus bound, by themselves, do not
exclude any positive value of \(\mu\).  In particular, (23)--(25) must not
be promoted into a positive-mismatch close.

There is a boundary-sensitive sharpening of (14).  If \(D\) is the
odd-degree set of a graph \(J\), then

\[
 \boxed{|J\cap\delta(T)|
 \le \min\{|D\cap T|,|D\setminus T|\}+2s_J.}       \tag{26}
\]

Indeed,

\[
 \sum_v\bigl(\deg_J(v)-\mathbf1_D(v)\bigr)=2|J|-|D|=2s_J.
\]

All summands are nonnegative.  Bounding a cut by the degree sum on either
shore proves (26).  At the present circle cut, \(|D\setminus
P_\Gamma|>|D\cap P_\Gamma|\) throughout the live range, so (25)--(26)
give the exact necessary boundary incidences

\[
 |D\cap P_\Gamma|+2s
 \ge \left\lceil{(p-2-\mu)(3-\mu)\over2}\right\rceil,
 \qquad \mu=0,1,2.                                \tag{27}
\]

In particular, a \(\mu=1\) completion must satisfy

\[
 |P_\Gamma\setminus D|\le2s,
\tag{28}
\]

so at most fourteen matched circle points can lie outside \(D\) in the
ordinary geometry, and at most eighteen in the triple geometry.  For
\(\mu=2\), (27) becomes

\[
 |D\cap P_\Gamma|\ge {p-3\over2}-2s.              \tag{29}
\]

This is a genuine coupling of the Boolean completion to the near-pencil
boundary, stronger than the size-only surplus bound.  It is still a
necessary condition rather than an exclusion: the near-pencil theorem
does not control which circle points belong to \(P_\Gamma\).

Finally, repeating the same circle operation cannot lower a positive
mismatch.  Switch to the Boolean shadow \(h\), putting

\[
 K'=\operatorname{diag}(h)K\operatorname{diag}(h),\qquad
 v'=\operatorname{diag}(h)z=\mathbf1+2\mathbf1_{N_\Gamma}.
\]

The sparse vector

\[
 d'=-\operatorname{diag}(h)d
\]

satisfies \(K'd'=pd'\); it is positive on all of \(N_\Gamma\), positive
on \(P_\Gamma\), and negative precisely on the original three-set \(E\).
If \(\mu\ge3\), every triple of the new spike set \(N_\Gamma\) determines
the same circle \(\Gamma\), and its required orientation is exactly
\(d'\).  Its new mismatch is therefore three.  Moreover, switching back by
\(h\) sends \(v'-2d'\) to the original \(v\).  Thus the operation is the
involution

\[
 (E,3)\longleftrightarrow(N_\Gamma,\mu),           \tag{30}
\]

not a descending mismatch iteration.  For \(\mu=1,2\) there is no triple
of new spikes on which to repeat it.  Equations (26)--(30) isolate the
remaining information needed: a close must rule out the boundary
near-containment (28), strengthen (29), or use a condition beyond this
single-circle involution.

## What this proves—and what it does not

This theorem replaces the unrestricted lower-endpoint boundary problem by
three explicit two-half geometries, and the sharpened global-cut argument
excludes their zero-circle-mismatch all-bad triples, including the explicit
Kiss--Somlai triangular family.  It does **not** exclude the same geometries
with a positive-mismatch completion.  Incidence parity, the directional
parallel-count quotas, and the matching/core decomposition alone are
therefore still not a proof of residual (ii); the remaining positive-
mismatch fibre requires another use of the global Boolean maximum or an
equivalent nonlinear condition.

## Same-day positive-mismatch follow-up

`NOTE_2026-09-02_CIRCLE_BOUNDARY_GEOMETRY.md` applies (28) to every square
circle.  For \(\mu=1\), all affine circles except the two removed pencil
lines are excluded for every live \(p\).  Generic finite circles are excluded
from \(p=73\) in the ordinary geometry and from \(p=89\) in the triple
geometry; noncentral circles through the pencil point are excluded from
\(p=53\) ordinary and \(p=61\) triple.  The two uniform survivors are exactly
the removed pencil lines and the finite square norm circles centered at the
pencil point.

The survivors have been audited rather than silently promoted to a close.
`NOTE_2026-09-02_REMOVED_LINE_MU1_CUT_AUDIT.md` gives every negative
circle-supported cut and shows local feasibility.  The centered case obeys
the exact saturation identity

\[
 \sum_{v\in P}\bigl(\deg_H(v)-\mathbf1_D(v)\bigr)
   =q+g+2(U+I)\le2s,
\]

proved in `NOTE_2026-09-02_CENTERED_MU1_SATURATION.md`.
Finally, `NOTE_2026-09-02_CIRCLE_COMPLEMENT_PHASE_PAIR.md` supplies two
unconditional Boolean eigenvectors from each circle and sharpens one member's
full-Max score according to its boundary phase.  These results locate the
remaining \(\mu=1\) blocker precisely: one needs a directional bound on the
two complement-word \(H\)-scores, or an equivalent simultaneous cut using
vertices outside the spike circle.  Repeating circle-incidence counts or the
conditional affine-alias argument cannot provide that bound.

`NOTE_2026-09-02_CENTERED_COMPLEMENT_DIRECTION_CAPACITY.md` also rules out
the separate-direction version of that idea: its exact joint sign table has
enough double-negative edges in every direction to absorb the full branch-B
parallel quota.  Any bound on the complement-word scores must therefore use
the **common cross-direction simple Radon lift**, not independent directional
capacity.

The global pencil sum also fails for a structural reason.  Weighting the
paired inequalities over all pencil directions by their eigensigns does make
the hard and opposite Max bounds point the same way, but the coefficient of
an edge \(uv\) contains the projective quartic correlation

\[
 \sum_{[g]}\eta\!\left(N(g)\det(g,u-a)\det(g,v-a)\right).
\]

It depends on the full triangle \((a,u,v)\), not only the edge direction or
the branch-B parallel count.  Proposition 15.634's multiple nonzero circle-
operator eigenspaces certify that this channel is not scalar, and its Weil
bound is too large after summing over \(|H|\asymp p^2\).  Thus an all-pencil
average does not bypass the common cross-direction lift either.
