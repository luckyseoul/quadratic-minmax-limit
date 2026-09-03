# Compact rays: odd blindness and the joint even-moment gate

**Status:** proved method barriers, separate degree-six/eight compatibility,
an all-prime bounded-compact odd-Radon theorem, and an exhaustive
arbitrary-compact local zero-form no-go for one row profile. Residual (ii)
remains open.

This note attacks the higher-moment gate left by Propositions 15.759--15.761
for the two compact survivor rays of Proposition 15.758. Its all-prime parts
are symbolic; Section 4 also records one explicitly classified finite
certificate for the one-compact/six-AE row first encountered at
$p=31,t=69$.

The exact outcome is:

- every odd moment is blind on both full rays;
- degree six passes on the full branch-B ray and lower branch-C endpoint;
- degree eight passes on those fronts as a **separate** projection;
- no one $\mathbf F_p$ labelling is yet known to pass degrees six and eight
  together;
- for $p=4r+3$, every row with $b$ arbitrary compact atoms and $r-1$
  all-equal atoms is central under zero odd forms whenever $3b\le r+2$;
- at $p=31$, no arbitrary compact atom plus six all-equal atoms can make
  that row zero in every odd degree and in degrees six and eight. The same
  profile appears in every balanced local allocation for $69\le t\le99$.

Thus the live congruence gate is the joint degree-six/eight system, not degree
ten. The new $p=31$ result is a balanced-profile **zero-global-form** method
no-go, not a global obstruction. The dominance construction below exists only
over an algebraic closure or a finite extension: it does not supply admissible
$\mathbf F_p$ labels or form coefficients. Nonzero forms coupled across
directions therefore remain a live possibility, not a constructed finite-field
lift. The signed Boolean affine-box intersection is still undecided.

## 1. Every odd moment is rowwise blind

For $2\le d\le p-1$ and $0\le k<d/2$, put

\[
 Q_{d,k}(s,t)=(s-t)^2(st)^k(s+t)^{d-2-2k}.                 \tag{1}
\]

For odd $d$,

\[
 Q_{d,k}(a,-a)=0,\qquad
 Q_{d,k}(a,0)+Q_{d,k}(-a,0)=0.                            \tag{2}
\]

Hence a compact atom on $(a,-a,0)$, with $0$ distinguished, and an
all-equal atom on the same labels both have zero odd moment. For
$d\le p-2$, a unit star is a field sum of a polynomial of degree at most
$p-2$ and also vanishes. Therefore the omitted-pair graph
$-S_a-S_{-a}+\{a,-a\}$ vanishes too.

These identities are atomwise and independent of every count and direction
sign. Both full rays satisfy all odd rows

\[
 d=3,5,\ldots,p-2.                                        \tag{3}
\]

In particular degree five cannot exclude either ray.

## 2. Degree six

Writing $x=a^6$, the three channels $k=0,1,2$ are

\[
\begin{array}{c|ccc}
\text{atom on }(a,-a,0)&0&1&2\\ \hline
\text{compact, }0\text{ distinguished}&-2x&0&4x\\
\text{all-equal}&2x&0&4x\\
\text{omitted pair}&0&0&4x.
\end{array}                                                \tag{4}
\]

For $H_6=\{a^6:a\in\mathbf F_p^*\}$, Cauchy--Davenport gives

\[
 |N H_6|\ge\min\{p,N|H_6|-(N-1)\}.                        \tag{5}
\]

### Full branch B, $p=4r+1$

Choose a linear form $s_L=L(v)$ with its one zero in a hard direction.
Give a hard compact group sixth-power sum $-s_L^6/2$. On an opposite row
let $b_L=Q_L-r$. Since $-1$ is a sixth power, give its compact atoms sum
$0$ for even $b_L$ and $s_L^6$ for odd $b_L$. Give the omitted pair
$s_L^6$ or $-s_L^6$, respectively, and the $r-2$ all-equal atoms minus
one half the omitted value. Including the outer opposite sign gives

\[
 F_{6,0}=L^6,\qquad F_{6,1}=0,\qquad F_{6,2}=-2L^6.        \tag{6}
\]

Equation (5) supplies the exact sums. Its only live short exception is
$p=37$, where

\[
 H_6=\{1,10,11,26,27,36\},                                \tag{7}
\]

and the needed targets $0,\pm1/2$ at counts six and seven are:

\[
\begin{array}{c|l}
0\ (6)&1+1+1+36+36+36\\
0\ (7)&1+1+1+1+11+11+11\\
18=-1/2\ (6)&1+1+1+26+27+36\\
18\ (7)&1+1+1+1+26+26+36\\
19=1/2\ (6)&1+1+1+1+26+26\\
19\ (7)&1+1+1+1+26+27+36.
\end{array}                                                \tag{8}
\]

Appending $1+36=0$ handles every larger count. The full branch-B ray
therefore passes every odd moment and degree six with one labelling.

### Lower branch-C endpoint, $p=4r+3$

Here the hard compact counts are $r-3,r-2$, the opposite all-equal count is
$r-1$, and the opposite compact count is zero. Make all sixth-power sums
zero. Equation (5) covers every live order except $p=31,r=7$, where use

\[
(0,1,8),(0,-1,-8),(1,5,14),(-1,-5,-14)                    \tag{9}
\]

for the four compact atoms, antipodal scales

\[
1,2,4,8,3                                                     \tag{10}
\]

for five compact atoms, and antipodal scales

\[
1,1,1,4,8,3                                                   \tag{11}
\]

for six all-equal atoms. Their degree-six sums are zero; (9) is two
negation pairs and (10)--(11) are antipodal, so all odd moments vanish too.

## 3. Degree eight separately

Writing $x=a^8$, the four channels are

\[
\begin{array}{c|rrrr}
\text{atom on }(a,-a,0)&0&1&2&3\\ \hline
\text{compact, }0\text{ distinguished}&-2x&0&0&-4x\\
\text{all-equal}&2x&0&0&-4x\\
\text{omitted pair}&0&0&0&-4x.
\end{array}                                                \tag{12}
\]

Let $H_8=(\mathbf F_p^*)^8$.

### Full branch B

For every live prime $p=4r+1$, there exist $a,c\in H_8$ with

\[
 c+2a=1.                                                   \tag{13}
\]

This is symbolic. If $p\equiv1\pmod {16}$, take $(a,c)=(1,-1)$. If
$p\equiv5\pmod8$, write $p=s^2+t^2$, $s\equiv1\pmod4$; the number of
solutions is the positive order-four cyclotomic number

\[
 (1,0)_4=(3,0)_4=\frac{p-3-2s}{16}.                       \tag{14}
\]

If $p\equiv9\pmod {16}$, write $p=x^2+4y^2$, $x\equiv1\pmod4$. According
as 2 is or is not quartic, the relevant order-eight number is

\[
 \frac{p-15-2x}{64}
 \quad\text{or}\quad
 \frac{p-7+6x}{64}.                                       \tag{15}
\]

The bounds $|s|,|x|\le\sqrt p$ prove positivity in the live range; the only
boundary check is $p=41,x=5$, where the second count is one. Formulas
(14)--(15) are from the cyclotomic appendix of
[Huczynska--Johnson, arXiv:2201.07553](https://arxiv.org/abs/2201.07553).

One pair from (13) handles every exact compact count $b$:

\[
 1=c^b+2a(1+c+\cdots+c^{b-1}).                            \tag{16}
\]

After multiplication by $L(v)^8$, use $c^bL(v)^8$ for the omitted pair and
the $b$ terms $ac^jL(v)^8$ for compact atoms. The all-equal sum is
$-c^bL(v)^8/2$ and the hard sum is $-L(v)^8/2$. This yields

\[
 F_{8,0}=L^8,\qquad F_{8,1}=F_{8,2}=0,\qquad F_{8,3}=2L^8. \tag{17}
\]

Cauchy--Davenport covers the exact sums except for these displayed
eighth-power-residue witnesses:

\[
\begin{array}{c|c|l}
p&N&\text{target}\\ \hline
29&4&-1/2:1+1+16+25;\quad0:1+1+7+20\\
41&7&-1/2:1+1+1+1+10+10+37;\quad0:1+1+1+1+1+18+18\\
41&8&-1/2:1+1+1+1+1+1+18+37;\quad0:1+1+1+1+1+10+10+16\\
41&9&-1/2:1+1+1+1+1+10+10+18+18;\quad0:1+1+1+1+1+1+1+16+18.
\end{array}                                                \tag{18}
\]

The zero witnesses cover the hard direction where $L(v)=0$. This labelling
also kills every odd moment, but its sixth-power sums are not known to be
those in Section 2.

### Lower branch-C endpoint

For $p\equiv3\pmod4$, $\gcd(8,p-1)=2$, so $H_8$ is the nonzero-square
subgroup and Cauchy--Davenport gives $3H_8=\mathbf F_p$. All exact counts
$r-3,r-2,r-1$ are at least three, so they can separately be assigned
eighth-power sum zero. This does not show the same scales satisfy degree six.

## 4. Interior branch C and the first joint gate

For centered antipodal degree-six atoms define

\[
 H=F_{6,2}+2F_{6,0}.                                      \tag{19}
\]

It vanishes on all $m=2r+2>6$ hard directions, hence identically. On an
opposite row, if $A$ and $B$ are the all-equal and compact sixth-power sums,

\[
 F_{6,0}=-2A+2B,\quad F_{6,2}=-4A-4B,\quad H=-8A.          \tag{20}
\]

Thus $A=0$. A compact count $b=0$ forces $F_{6,0}=0$, while $b=1$ forces
$F_{6,0}=2a^6\ne0$. With

\[
 \delta=t-(2r^2-4r-2)=\sum_L b_L,                         \tag{21}
\]

the balanced allocation has $m-\delta$ zero rows when $1\le\delta<m$.
Consequently $1\le\delta\le m-7$ defeats this **centered** construction.
The first case is $p=31,t=69$, with $b$-multiset $\{1,0^{15}\}$.

This is not an arbitrary-label obstruction. For a compact atom on $(z,1,0)$,

\[
 Q_{6,2}^{\rm comp}+2Q_{6,0}^{\rm comp}
 =z(z+1)^2(4z^2-9z+4),                                   \tag{22}
\]

so $z=-1$ is a special root.

Indeed $p=31,t=69$ has an escape passing every odd degree and degree six.
Use compact atom $(1,-1,0)$ and six all-equal atoms in three negation pairs:

\[
\begin{split}
 &(0,1,2),(0,-1,-2),\\
 &(0,1,24),(0,-1,7),\\
 &(0,4,23),(0,-4,8).
\end{split}                                                \tag{23}
\]

The compact degree-six vector is $(29,0,4)$. The three representatives have
vectors $(22,18,4),(24,23,5),(17,21,20)$, summing to $(1,0,29)$.
Doubling and adding the compact gives zero; negation kills all odd moments.
The other rows reuse (9)--(11), so the whole balanced allocation passes all
odd moments and degree six.

Block (23) has degree-eight vector $(22,17,18,26)$ and does not extend
directly. We can now exclude **every** extension that retains its centered
compact atom and requires the exceptional row itself to remain zero in all
odd degrees and in degrees six and eight.

### Odd Radon compression

Let $M(E)$ be the aggregate multiplicity of an edge $E={s,t}$ in the six
all-equal triangles. For a nonantipodal edge put

\[
 U=(s+t)^2,\qquad D=(s-t)^2,\qquad V=st=\frac{U-D}{4},
 \qquad \sigma=s+t.                                      \tag{24}
\]

The pair $E,-E$ has the same $(U,D,V)$ and opposite $\sigma$. If
$n_E=M(E)-M(-E)$, the odd moments $d=3,5,\ldots,29$ say exactly

\[
 \sum_{(U,D)\in\Omega}W(U,D)P(U,V)=0
 \quad(\deg P\le13),\qquad
 W(U,D)=n_E D\sigma,                                     \tag{25}
\]

where $\Omega=H\times H$ and $H$ is the set of 15 nonzero squares in
$\mathbf F_{31}$. The support of $W$ has size at most 18.

Here is an elementary dual Reed--Muller support argument. If a support
$S$, $|S|\le18$, were not concentrated on a line, fix $P\in S$. Unless a
line through $P$ contains at least 15 support points, the other at most 17
points can be covered by at most

\[
 \max\left\{13,\left\lceil\frac{17}{2}\right\rceil\right\}=13
 \tag{26}
\]

lines avoiding $P$: pair points from different radial classes and leave at
most the largest radial class as singletons. The product of their line
equations has degree at most 13 and isolates $P$, contradicting (25). Thus
every support point lies on a 15-point support line, which forces all of
$S$ onto one line.

Direct character summation, replayed exactly in the certificate, gives the
intersection distribution for all 992 affine lines:

\[
 |\ell\cap\Omega|=0,7,8,15
 \quad\text{with multiplicities}\quad 47,675,225,45.     \tag{27}
\]

The 45 maximal lines are precisely

\[
 U=u_0,\qquad D=d_0,\qquad D=aU
 \quad(u_0,d_0,a\in H).                                  \tag{28}
\]

On their 15 square parameters, Vandermonde interpolation makes the unique
degree-$13$-orthogonal weights proportional to the parameter itself.
For $D=d_0$ or $D=aU$, this forces $n_E$ to traverse all 15 nonzero sign
classes, impossible because $|n_E|\le6$. For $U=u_0$, it forces $n_E$ to
be a nonzero constant. Absolute value at least two needs at least 30 edge
occurrences. Absolute value one requires the 15 fixed-sum edges, a perfect
matching on 30 vertices; the three remaining edge occurrences cannot turn
all 30 odd degrees even. Therefore $W=0$ and

\[
 M(E)=M(-E)\quad\text{for every edge }E.                 \tag{29}
\]

This is aggregate edge symmetry, not trianglewise pairing.

### Classification of the six triangles

Compare the six-triangle multiset $X$ with $-X$. Equation (29) makes
$(X,-X)$ a $2$-$(v,3)$ trade. Cancel common blocks. A nonempty reduced
trade has volume at least four; volume four is the unique Pasch trade,
volume five does not exist, and volume six has four simple types: the
6-cycle, semihead, trade-X, and trade-Y. The small-trade templates are
listed in [Grannell--Griggs--Knor--Pisanski,
*Small surface trades in triangular embeddings*](https://grannell.net/Papers/surftrade.pdf);
the minimum-volume and no-volume-five facts go back to Hwang and are also
summarized in [Khosrovshahi--Mahmoodian,
*Classification of simple 2-(6,3) and 2-(7,3) trades*](https://ajc.maths.uq.edu.au/pdf/19/ocr-ajc-v19-p55.pdf).

There is no repeated-block exception at these volumes. If a positive block
$T$ has multiplicity $m\ge2$, the negative leg needs at least $m$ blocks
through each of the three pairs of $T$, hence volume at least $3m$. The only
remaining case is $m=2$ and volume six. Its six negative blocks then give
the three vertices of $T$ total incidence 12. On the positive side, after
two copies of $T$, none of the four remaining blocks may contain a pair of
$T$, so their total incidence on those vertices is at most four: the total
is at most $6+4=10$, a contradiction.

Consequently there are only three cases.

1. The reduced trade is empty. The six blocks consist of three negation
   pairs; two pairs and two invariant centered triangles; one pair and four
   invariant triangles; or six invariant triangles.
2. The reduced trade is a Pasch four-trade, plus either one negation pair or
   two invariant triangles.
3. The reduced trade is one of the four volume-six types.

All three cases have exact finite certificates. For case 1, three pair
representatives would need sum

\[
 (1,0,29\mid1,0,0,2).                                     \tag{30}
\]

The calculation checks all $\binom{31}{3}=4495$ distinct-label triangles,
10,104,760 unordered pair sums, 2,543,460 distinct pair-sum vectors, and the
120, 925, and 961 invariant-triangle sum states of sizes two, four, and six.
No category completes the target. Its deterministic evidence hash is

~~~text
26bea31c9906b005ff4fc1dc0121d43eb07ef7f62369b90b902026ae0d293c95
~~~

For case 2, one Pasch parametrization is

\[
 (a,-a,b),\ (a,-b,c),\ (-a,-b,-c),\ (b,c,-c).            \tag{31}
\]

Its negative is the mate leg. All seven fixed-point-free involutions from a
standard Pasch leg to its mate give 208,537 parameter assignments, 165,660
internally valid assignments, and 6,910 distinct concrete block multisets,
including quotient collisions and repeated blocks. None can be completed by
one negation pair or two invariant triangles. The residual-transcript hash is

~~~text
40889ecbc7e92660d045e547a7f532b1aaa1dcf5519c9185ef02f0f3eea910ce
~~~

For case 3, every isomorphism from the positive to negative leg of each of
the four templates is imposed through
$\lambda(\pi v)=-\lambda(v)$. Allowing all cross-vertex collisions gives
3,042,008 assignments, 2,164,860 internally valid assignments, and 169,940
distinct concrete six-block multisets. None hits the joint degree-six/eight
target. The concrete-multiset hash is

~~~text
78ee0fc05757a9d332a8d2da3605a921b28207aafca746c806bf17e043f26dd0
~~~

Thus, with compact $(1,-1,0)$ fixed, no six all-equal triangles can make this
one row zero in every odd degree through 29 and in degrees six and eight.
This is exhaustive for that **local zero-global-form mechanism** with the
centered compact fixed. It does not by itself exclude a noncentered compact
atom, nonzero global forms coordinated across directions, or another global
moment realization. The next subsection extends the odd reduction to
arbitrary compact labels, and the finite certificate following it closes the
first of those escapes.

### Arbitrary compact labels in the odd Radon step

The odd-moment reduction itself does extend to a noncentered compact atom.
Write its signed edge chain as

\[
 K(a,b;c)=\{a,b\}-\{a,c\}-\{b,c\}.
\]

Together with six positive triangle boundaries this has at most 21 signed
edge occurrences, so the support of the same function $W$ has size at most
21. The line-isolation argument above works through support 27. Indeed, for
a chosen support point, points in different radial classes can be paired on
lines avoiding it. If no radial class has more than 13 points, at most 13
lines suffice. If one class has at least 14 points, choose the point to
isolate off that long line; that line and the remaining off-line singletons
again use at most 13 factors. Thus a nonzero $W$ is still supported on one
of the 45 maximal lines in (28).

One all-equal atom contributes at most one to any edge-orbit difference,
whereas $K(a,b;c)$ contributes at most two. Hence now $|n_E|\le8$.
Horizontal and diagonal lines remain impossible: their weights require all
15 projective nonzero classes, while integer representatives in $[-8,8]$
give only eight.

On a vertical line the invariant condition is
$n_E(s+t)=\lambda$. Scaling the labels reduces its fixed edge sum to one.
The 21-occurrence bound makes $|n_E|=1$ at its 15 points. Start with the
corresponding fixed-sum matching. The only negative
physical occurrences are the two negative compact edges, so at most two
matching edges can be reversed to their negatives. An exact parity check
(or the four-endpoint toggle bound) shows that at least 22 of the original
30 matching vertices remain odd. After subtracting the signed 15-edge
baseline, the residual chain has coefficient $\ell_1$ mass at most six, so
its parity support can toggle at most 12 vertices. At least ten odd
vertices survive, contradicting the fact that the total chain is the
mod-two sum of seven triangle boundaries. Therefore the total signed edge
chain is centrally symmetric: every nonantipodal edge-orbit difference is
zero.

This does **not** force the compact atom to be centered. For example, the
noncentered compact $(0,1;2)$ and the six all-equal triangles

\[
 (0,2,30),(0,2,30),(0,1,2),(0,1,29),(0,29,30),(1,2,29)
\]

have central total edge chain and annihilate every odd row through degree
29. Their degree-six and degree-eight vectors are respectively
$(26,26,5)$ and $(13,4,30,6)$, so this example does not solve the joint even
gate. The remaining arbitrary-compact local problem is now the exact finite
central-chain completion problem, not an unrestricted odd-moment search.

### All-prime bounded-compact odd-Radon centrality

The preceding centrality reduction is neither special to 31 nor limited to
one compact atom. Let $p=4r+3$, $r\ge7$, put
$h=2r+1=(p-1)/2$, and consider a row made from $b$ arbitrary compact atoms
and $r-1$ positive all-equal triangles. Its signed edge chain $C$ has at
most

\[
 N=3(r+b-1)                                                \tag{29a}
\]

signed edge occurrences. For a nonantipodal edge, retain the coordinates in
(24), orient its antipodal orbit, and put $n_E=C(E)-C(-E)$. The odd channels
$d=3,5,\ldots,p-2$ say exactly that

\[
 W(U,D)=n_ED(s+t)
\]

is orthogonal on $\Omega=H\times H$ to every polynomial in $U,D$ of total
degree at most $2r-1$.

A nonzero such word with support at most $4r-1$ must lie on one affine line.
Indeed, if no line contains $h$ support points, the other points relative to
a chosen support point split into radial classes of size at most $2r-1$.
Pairing different classes covers them with at most $2r-1$ lines avoiding the
chosen point; the product of their equations isolates it. Once a line
contains $h$ support points, isolating any off-line point forces the entire
support onto that line. Thus the support theorem applies exactly when

\[
 N\le4r-1\quad\Longleftrightarrow\quad 3b\le r+2.          \tag{29b}
\]

For a nonvertical line $D=AU+B$ with $AB\ne0$, direct character expansion
gives

\[
 \#\{U\in H:AU+B\in H\}
 =\frac{p-2-\chi(A)-\chi(B)-\chi(-B/A)}4.
\]

Because $\chi(-1)=-1$, this is $r+1$ when $A,B$ are both nonsquares and
$r$ otherwise. Hence the only $h$-point lines are

\[
 U=u_0,\qquad D=d_0,\qquad D=aU
 \quad (u_0,d_0,a\in H).
\]

On the $h$ square parameters the unique degree-$(2r-1)$ orthogonal weight
is proportional to the parameter itself. One positive all-equal triangle
changes a fixed orbit difference by at most one and one compact atom by at
most two, so

\[
 |n_E|\le B=(r-1)+2b.                                    \tag{29c}
\]

On $D=d_0$ and $D=aU$, the values $n_E$ must traverse all $h$ projective
nonzero classes. This is impossible when

\[
 B<h\quad\Longleftrightarrow\quad 2b<r+2,                \tag{29d}
\]

which follows from (29b). On $U=u_0$, the congruence instead makes all
$n_E$ equal to one nonzero integer $k$ on the fixed-sum matching of $h$
edges. Here $B<h$ gives $2B<p$, so the bounded congruent values really are
one integer. Also

\[
 N<2h\quad\Longleftrightarrow\quad3b<r+5,
\]

and hence $|k|=1$.

Fix the sign of $k$. Every one of the $h$ matching orbits then needs an edge
occurrence aligned with that sign. An all-equal triangle supplies at most
one: two of its pair sums cannot equal the same fixed nonzero sum unless two
labels coincide. A compact atom supplies at most two: for

\[
 K(a,b;c)=\{a,b\}-\{a,c\}-\{b,c\},
\]

alignment of all three signed edges would force $a=b$. Thus the entire row
has aligned capacity at most

\[
 (r-1)+2b=B<h,                                           \tag{29e}
\]

a contradiction. This aligned-incidence step covers the boundary case
$3b=r+2$ as well.

There is an independent parity check for the strict interior. At most one
matching edge can reverse per compact atom: if both negative compact edges
reversed, their equal fixed sums would again force $a=b$. The selected
baseline therefore has at least $4r+2-4b$ odd vertices. Group its selected
occurrences by the $r+b-1$ underlying triangles. A triangle with zero or
three selected edges is Eulerian, while one or two selected edges have only
two odd endpoints. Hence its parity support is at most $2(r+b-1)$. The
strict inequality

\[
 4r+2-4b>2(r+b-1)\quad\Longleftrightarrow\quad3b<r+2
\]

is already contradictory. At equality the sharper aligned count (29e)
finishes the proof. Therefore, under (29b), $W=0$. Since
$|n_E|\le B<p$, this lifts from characteristic $p$ to the integer equality

\[
 C(E)=C(-E)\qquad\hbox{for every edge orbit}.             \tag{29f}
\]

This gives an explicit balanced branch-C band. Put

\[
 t_0=2r^2-4r-2,\qquad m=2r+2,\qquad
 \delta=t-t_0=\sum_L b_L.
\]

The deterministic balanced allocation has
$b_L\in\{\lfloor\delta/m\rfloor,\lceil\delta/m\rceil\}$.
Writing

\[
 b_{\max}=\left\lfloor{r+2\over3}\right\rfloor,
\]

every opposite row is centrally symmetric whenever all odd global forms
vanish throughout

\[
 0\le\delta\le m b_{\max},\qquad
 t_0\le t\le t_0+m b_{\max}.                            \tag{29g}
\]

The next balanced profile has one row with $b=b_{\max}+1$, so (29g) is the
exact initial band certified by this support argument. It is a theorem about
the displayed balanced allocation, not every distribution with the same
total compact count.

For $p=31$, $r=7$, $m=16$, and $t_0=68$. Every balanced profile with

\[
 69\le t\le99                                             \tag{29h}
\]

contains a row with $b=1$: the counts are zero/one through $t=83$, all one
at $t=84$, and one/two through $t=99$; at $t=100$ they are all two. The
exhaustive certificate below depends only on the one-compact/six-AE atom
profile, not on $t$. It therefore excludes simultaneous zero odd,
degree-six, and degree-eight global forms for every balanced local profile
in (29h), not only the first occurrence $t=69$.

These are structural reductions, not residual-(ii) closure. They do not
treat nonzero odd forms, unbalanced allocations, common nonzero
degree-six/eight forms, or the signed Boolean lift.

### Full-balanced-ray maximal-line exclusion

The support-isolation step in (29b) stops when $3b>r+2$, but the exclusion
of a maximal line, *conditional on the support occupying one*, extends much
farther.  In fact it covers the full balanced range $0\le b\le r$.  Keep
$h=2r+1$ and write

\[
 N=3(r+b-1)\le 6r-3=3h-6<3h.             \tag{29i}
\]

On a line $D=d_0$ or $D=aU$, the unique degree-$(h-2)$ dual relation makes
the $h$ nonzero integers $n_E$ represent all classes in
$\mathbf F_p^*/\{\pm1\}$.  Choosing the least absolute integer in each
class gives the sharp lower bound

\[
 \sum_E|n_E|\ge1+2+\cdots+h={h(h+1)\over2}>3h-6.         \tag{29j}
\]

Thus neither of these two maximal-line types can occur, without requiring
the old pointwise bound $B<h$.

It remains to strengthen the vertical-line argument.  Choose a square root
$\sigma$ of $u_0$.  After orienting every orbit by the edge whose sum is
$\sigma$, the line relation says

\[
 n_E\equiv\kappa\pmod p                                  \tag{29k}
\]

for one nonzero residue $\kappa$.  Put
$a=\min(|\kappa|,p-|\kappa|)$, so $1\le a\le h$.  From (29i),
$ha\le N<3h$, hence $a\le2$.  If even one coordinate uses a different
integer lift of (29k), its absolute value is at least $p-a$, and the total
mass is at least

\[
 (h-1)a+(p-a)\ge (h-1)+(p-1)=3h-1>N.                    \tag{29l}
\]

Consequently all $h$ orbit differences are the same integer $k$, with
$|k|=1$ or $2$.

The unit case has a parity obstruction that survives arbitrary $b$.
Reduce the signed chain modulo two and project its vertex set by
$x\sim-x$.  The signs in a compact atom disappear modulo two, so both an
all-equal atom and a compact atom project to a triangle boundary.  Therefore
the projected chain is Eulerian.  For a nonantipodal orbit the projected
coefficient is

\[
 C(E)+C(-E)\equiv C(E)-C(-E)=n_E\pmod2;
\]

self-antipodal edges become loops and do not affect degrees.  The $h$
fixed-sum edges $\{x,\sigma-x\}$ have degree two at every projected vertex
except $[0]$ and $[\sigma/2]$, where they have degree one.  Thus their
projected graph has exactly two odd vertices, contradicting Eulerianity when
$|k|=1$.

If $|k|=2$, every one of the $h$ orbits needs at least two edge occurrences
aligned with the sign of $k$.  As in (29e), an all-equal atom supplies at
most one aligned occurrence and a compact atom at most two.  Hence

\[
 2h\le\hbox{aligned occurrences}\le(r-1)+2b
 \le3r-1<4r+2=2h,                                      \tag{29m}
\]

again impossible.  All three maximal-line supports are therefore excluded
for every $0\le b\le r$.

This is deliberately only a **conditional maximal-line exclusion**.  Past
$3b=r+2$, (29b) no longer forces a nonzero word onto one line; conic and
cubic low-weight configurations remain possible.  Thus (29i)--(29m) do not
extend the centrality interval (29g) by themselves.

#### Couvreur peeling for a line contained in a larger support

There is nevertheless an exact way to control the distinction between
"supported on a line" and "contains a line."  We use A. Couvreur,
*The dual minimum distance of arbitrary-dimensional algebraic-geometric
codes*, J. Algebra **350** (2012), 84--107, Theorem 3.8 (arXiv
`0905.2345`, DOI `10.1016/j.jalgebra.2011.09.030`).  For $d\ge2$ its
first two linked-configuration thresholds say that the smallest
$d$-linked set has size $d+2$ exactly when $d+2$ points are collinear, and,
if there is no such line, has size at least $2d+2$.  Lemma 2.13 of that
paper says linkedness and minimal linkedness are preserved by algebraic
field extension.  We may therefore homogenize the distinct points of
$\Omega\subset\mathbf A^2(\mathbf F_p)$ and apply the theorem over the
algebraic closure.  Repeated edge occurrences cause no problem: they have
already been aggregated into one nonzero coefficient at each support
point.

Let $m=h-2$, let $S$ support a degree-$m$ dual word, and assume
$|S|\le3m=3h-6$.  Suppose $S$ contains $h=m+2$ points on a line $L$.  If
$R=S\setminus L$ is nonempty and $\ell_L$ defines $L$, then

\[
 (W(P)\ell_L(P))_{P\in R}                               \tag{29n}
\]

is a full-support degree-$(m-1)$ dual word: test the original relation
against $\ell_LG$ for every $\deg G\le m-1$.  But

\[
 |R|\le2h-6<2(m-1)+2=2h-4.                              \tag{29o}
\]

Couvreur's thresholds force $h-1$ points of $R$ onto a second line $L_2$.
If $Q=R\setminus L_2$ were nonempty, multiplying (29n) by an equation of
$L_2$ would make $Q$ degree-$(m-2)$ linked.  Hence

\[
 |Q|\ge(m-2)+2=h-2,\qquad
 |S|\ge h+(h-1)+(h-2)=3h-3,                             \tag{29p}
\]

contrary to the support bound.  Thus $R\subset L_2$ and
$h-1\le|R|\le h$.  A nonmaximal line meets $H\times H$ in at most
$r+1<h-1$, so $L_2$ is maximal.  Every line-containing support in the
Couvreur range has now been reduced to one maximal line or to two maximal
lines, with no census and no minimal-support assumption on $W$.

#### Coefficients on the two maximal lines

The two-line alternatives can also be excluded throughout $0\le b\le r$.
The Cartesian duality formula for $H$, the roots of $X^h-1$, is

\[
 C_{H^2}(h-2)^\perp
   =\operatorname{diag}\!\left({UD\over h^2}\right)C_{H^2}(h-1).
                                                               \tag{29q}
\]

Thus write $W(U,D)=UDf(U,D)/h^2$ with $\deg f\le h-1$; after orienting by
$\sigma^2=U$ one has $n_E=f(U,D)\sigma/h^2$.

First take maximal lines from different families.  They meet at a point of
$H^2$, so their union has $2h-1$ points.  A degree-$(h-2)$ form vanishing
on the union vanishes identically on each component and is divisible by
their product.  The evaluation rank is therefore

\[
 { (m+1)(m+2)-(m-1)m\over2}=2m+1=2h-3,                 \tag{29r}
\]

and the dual space on the union has dimension two.  It is spanned by the
two individual line relations.  At least one component is nonvertical; on
its $h-1$ points away from the intersection, $n_E$ occupies $h-1$ distinct
projective residue classes.  Consequently

\[
 \sum_E|n_E|\ge1+\cdots+(h-1)={h(h-1)\over2}>3h-6.       \tag{29s}
\]

Now take two lines in the same family.  Since $f$ vanishes on the other
$h-2$ lines of that family, it factors as their product times one affine
linear factor.  For two vertical lines this reads

\[
 f=P(U)(AU+BD+C),qquad \deg P=h-2.                     \tag{29t}
\]

If $B\ne0$, the values of $n_E$ on the full first component are affine and
injective in $D\in H$.  The least absolute mass of $h=2r+1$ distinct
nonzero residues is

\[
 2(1+\cdots+r)+(r+1)=(r+1)^2>6r-3=3h-6.                \tag{29u}
\]

If $B=0$, both nonempty components are full and their coefficients are
constant.  Their canonical absolute values have sum at most
$\lfloor N/h\rfloor\le2$, so both are units.  Replacing even one unit by
its other integer lift costs at least $4h-1>N$, hence the actual integers
are constant units.  Modulo two, the two projected fixed-sum graphs cancel
their common odd vertex $[0]$ and retain the distinct odd vertices
$[\sigma_1/2]$ and $[\sigma_2/2]$, contradicting Eulerianity.

For two horizontal lines, (29t) with $U,D$ exchanged shows on the full
first component that equality of projective coefficient classes is governed
by

\[
 z=U(A+BU)^2.                                            \tag{29v}
\]

For two diagonals, factoring
$\prod_{a\in H\setminus\{a_1,a_2\}}(D-aU)$ and using $U^h=1$ gives

\[
 zU^3=(A+BU)^2.                                         \tag{29w}
\]

Each equation has at most three solutions in $U$ for a fixed nonzero
projective class.  If $h=3q+s$, $0\le s<3$, the least possible integer mass
when every class is repeated at most three times is

\[
 L_3(h)=3{q(q+1)\over2}+s(q+1)>3h-6\qquad(h\ge15).       \tag{29x}
\]

This excludes all same-family pairs as well.  Combining (29i)--(29x), no
realizable word of support at most $3h-6$ can contain $h$ collinear points
anywhere on the full balanced branch.  The remaining Couvreur alternatives
are a support containing $2h-2$ points on a conic but no $h$ collinear
points, or, at the exact boundary $3h-6$, a minimally linked cubic/degree-
$(h-2)$ complete intersection.  The conic coefficient exclusion remains
open.  The next subsection excludes the cubic for every $p\ge31$ but still does
not enlarge the centrality interval (29g).

#### Boundary cubic units exclude every cubic

The exact boundary cubic alternative admits an all-prime exclusion.
Let $m=h-2$ and suppose $S$ is the complete intersection of a cubic $F$ and
a degree-$m$ curve $G$ with no common component.  Its size is

\[
 |S|=3m=3h-6.                                            \tag{29y}
\]

The balanced occurrence budget reaches (29y) only for $b=r$, where
$N=3(r+r-1)=3h-6$.  Every support coordinate is a nonzero integer orbit
difference, so equality of support and $\ell_1$ budgets forces

\[
 |n_E|=1\qquad(E\in S).                                  \tag{29z}
\]

The $3m$ distinct points exhaust Bezout's intersection number.  Thus the
intersection is reduced and transverse, and the unique Cayley--Bacharach
relation has the Euler--Jacobi residue form

\[
 W(P)={\lambda\over J(F,G)(P)},qquad
 J(F,G)=F_U G_D-F_D G_U.                                 \tag{29aa}
\]

Assume first that $F$ is reducible, say $F=LQ$ for a line $L$ and a conic
$Q$.  A support point cannot lie at $L\cap Q$, where $F$ is singular.
Componentwise Bezout therefore puts exactly $m=h-2$ support points on $L$
and $2m=2h-4$ on $Q$.  Since

\[
 r+1<h-2,                                                 \tag{29ab}
\]

the line-intersection table following (29b) makes $L$ maximal.  Parameterize
it by $z\in H$.  Its support is $H\setminus\{a,b\}$, and

\[
 G|_L=c\,{z^h-1\over(z-a)(z-b)},qquad
 {1\over(G|_L)'(z)}=c' z(z-a)(z-b)                       \tag{29ac}
\]

at every support parameter.  Along $L$, (29aa) has the additional normal
factor $Q|_L(z)$.

If $L$ is vertical, $U=u_0$ and $D=z$, equations (29z)--(29ac) give

\[
 n_E=c\,{(z-a)(z-b)\over Q|_L(z)},qquad
 Q|_L(z)^2=c^2(z-a)^2(z-b)^2.                            \tag{29ad}
\]

The second identity holds at $h-2>5$ points while both sides have degree at
most four, hence it is a polynomial identity.  Therefore

\[
 Q(u_0,z)=c_0(z-a)(z-b).                                 \tag{29ae}
\]

The conic meets $L$ exactly at the two omitted grid points.

For a horizontal line, the corresponding squared identity is

\[
 Q|_L(z)^2=c^2z(z-a)^2(z-b)^2,                           \tag{29af}
\]

and for a diagonal it is

\[
 (z-a)^2(z-b)^2=c^2zQ|_L(z)^2.                           \tag{29ag}
\]

Again $h-2>5$ makes these polynomial identities.  Each is impossible at
$z=0$ by the odd multiplicity of the displayed factor $z$ (recall
$a,b\in H$ are nonzero).  Thus every line component of a reducible boundary
cubic must be vertical.

If $Q$ split into two further lines, the same argument would make all three
lines vertical.  But the product of the other two vertical equations is
constant on the first line, contradicting the nonconstant quadratic in
(29ae).  Hence the only remaining reducible candidate has the form

\[
 F=(U-u_0)Q(U,D),\qquad Q\text{ absolutely irreducible}. \tag{29ah}
\]

The conic contains its $2h-4$ support points plus the two grid points in
(29ae), so

\[
 |Q\cap(H\times H)|\ge2h-2=2m+2.                        \tag{29ai}
\]

Every line component carries $m>2$ rational support points, hence descends
to $\mathbf F_p$; the remaining conic does as well.  Section 2 of
`NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md` applies to (29ai) and gives
the tangent normal form

\[
 U=u z^2,\qquad D=d(z-1)^2,\qquad
 Q\cap(H\times H)=\{z\in\mathbf F_p:z\ne0,1\}.             \tag{29aj}
\]

Write the vertical line as $U=u\gamma^2$.  Its two intersections with $Q$
are $z=\gamma,-\gamma$.  The conic component has $2m=p-5$ support
points, whereas (29aj) has $p-2$ grid parameters.  Thus there is exactly one
further omitted parameter $e\ne0,1,\pm\gamma$, and, with

\[
 T=\mathbf F_p\setminus\{0,1\},\qquad
 P_T(x)=\prod_{t\in T}(x-t)={x^p-x\over x(x-1)},
\]

the restriction $g(z)=G(uz^2,d(z-1)^2)$ is

\[
 g(z)=c\,{P_T(z)\over(z^2-\gamma^2)(z-e)}.                 \tag{29ak}
\]

At every root of $g$ in the support,

\[
 {1\over g'(z)}=c' z(z-1)(z^2-\gamma^2)(z-e),             \tag{29al}
\]

because $P_T'(z)=-1/(z(z-1))$.  For the implicit equation
$Q=(U/u+1-D/d)^2-4U/u$, its Hamiltonian tangent
$(Q_D,-Q_U)$ is a constant multiple of
$(U',D')=(2uz,2d(z-1))$.  Hence $J(Q,G)|_Q$ is a constant multiple of
$g'(z)$.  Since $F=(U-u\gamma^2)Q$, on $Q$ one has

\[
 J(F,G)=u(z^2-\gamma^2)J(Q,G).
\]

The line factor cancels the matching factor in (29al).  Equations
(29aa) and (29z), with the coherent square root $\sigma=\sqrt u\,z$, now
give

\[
 W(z)=c_0z(z-1)(z-e),\qquad
 n(z)=c_1{z-e\over z-1}.                                  \tag{29am}
\]

Squaring $|n(z)|=1$ gives

\[
 (z-e)^2=c_2(z-1)^2                                      \tag{29an}
\]

at all $p-5>2$ conic-support parameters.  Both sides have degree two, so
(29an) is a polynomial identity.  Its leading and linear coefficients force
$c_2=1$ and $e=1$, contradicting $e\in T$.  Therefore **every reducible
boundary cubic is excluded**.  This is special to the boundary complete-
intersection weights and does not exclude the separate conic-supported
alternative.

The cubic may initially be supplied over the algebraic closure, but descent
causes no gap here.  It and its Frobenius conjugate both contain the more
than nine rational points of $S$.  If the cubic is geometrically integral,
Bezout forces the two cubics to coincide, so it is defined over
$\mathbf F_p$ up to scale.

For completeness, a geometrically integral singular cubic has at most
$p+2$ rational points, less than $3h-6$ for every $p\ge31$.  It remains to
treat a smooth cubic $C$.  The coordinate function

\[
                         U={X\over Z}\in\mathbf F_p(C)       \tag{29ao}
\]

is nonconstant and has degree at most three, because it comes from a pencil
of line sections of the plane cubic.  Its degree is not one on a genus-one
curve, so $\deg U\in\{2,3\}$.  Moreover $U$ is geometrically nonsquare.  If
$U=v^2$ over the algebraic closure, then
$\deg U=2\deg v$; the only numerical possibility would be
$\deg U=2$, $\deg v=1$, again impossible on a genus-one curve.  This argument
also covers the case where the coordinate vertex lies on $C$ and cancels one
zero and pole; no unsupported degree-three assumption is needed.

Normalize the connected double cover

\[
                   \pi:\widetilde C\longrightarrow C,
                   \qquad y^2=U.                           \tag{29ap}
\]

The zero and pole divisors of $U$ each have degree at most three, so at most
six geometric points have odd valuation.  Riemann--Hurwitz therefore gives

\[
             2g(\widetilde C)-2\le6,
             \qquad g(\widetilde C)\le4.                   \tag{29aq}
\]

Every point of $S\subset H\times H$ has $U$ a nonzero square and hence has
two distinct $\mathbf F_p$-points above it.  Consequently

\[
 \#\widetilde C(\mathbf F_p)\ge2|S|=3p-15.
\]

On the other hand, Weil and (29aq) give

\[
 \#\widetilde C(\mathbf F_p)\le p+1+8\sqrt p.
\]

These bounds contradict each other because

\[
 3p-15>p+1+8\sqrt p
 \quad\Longleftrightarrow\quad
 p-8>4\sqrt p,                                           \tag{29ar}
\]

which holds at $p=31$ and increases thereafter.  Thus the smooth case is
impossible as well.  Combining the reducible, singular, and double-cover
arguments, **every boundary cubic support is excluded for every
$p=4r+3\ge31$**.  The independent high-intersection conic alternative remains
open, so this result still does not assert centrality or residual-(ii)
closure.

### Exact arbitrary-compact central-chain certificate

There are

\[
 \binom{31}{2}(31-2)=13{,}485
\]

labelled compact atoms $(a,b;c)$, with the positive pair unordered.
Multiplication of all three labels by $\mathbf F_{31}^*$ preserves centrality
and the vanishing of every degree-six/eight channel. The labelled atoms split
into 450 scaling orbits: 449 free orbits of size 30 and the one centered orbit
of size 15. In lexicographic canonical order the centered representative is
index 435, $(1,30;0)$.

For each noncentered representative the exact search stores the sparse
edge-orbit imbalance and all seven even channels. At a nonzero imbalance
coordinate, some remaining triangle must contribute with the opposite sign,
so the search branches over every such triangle (28 or 29 candidates,
according to the orbit). With $r$ triangles left, the necessary prunes are

\[
 \|\partial\|_1\le3r,\qquad
 \|\partial\|_\infty\le r.
\]

States are memoized by depth, sparse boundary, and the seven moment values.
When a partial boundary first becomes zero, at most five blocks remain. The
complete trade classification through volume five leaves only an
invariant/negation-pair core or a Pasch four-trade, possibly plus one
invariant block.
The exact catalog sizes are

~~~text
invariant blocks:                  15
fixed-block sumsets, sizes 0..5:   1, 15, 120, 535, 925, 961
negation-pair vectors:             2,255
two-pair vectors:                  2,543,460
Pasch even vectors:                3,725
~~~

There is one final necessary-only acceleration. If one or two triangles
remain with required even vector $v$, then $2v$ must lie in the one-pair or
two-pair catalog, respectively. This cannot discard a completion because two
is invertible modulo 31 and even moments are invariant under global negation.

The exhaustive replay visits 317,916,856 DFS states. All 449 noncentered
orbits are infeasible. The sole centered orbit is delegated to the preceding
independently exhaustive centered theorem, which is also infeasible. Thus no
arbitrarily labelled compact atom plus six all-equal atoms makes this row zero
in all odd degrees $3,5,\ldots,29$ and in degrees six and eight.

The audited evidence chain is:

~~~text
source SHA-256:
1dcfce7b5765630655d049413c4d9138c544a6d05fe19e3308a9a20a2880d1f2
binary SHA-256:
4622dbcb2afbfdc4c0da3588e42cba86542a98763fb89429bed7ea2185915955
raw log SHA-256:
f3f77607181287095aa69644649d14d7b9b5e3a8f24044477b667549ef0512e3
normalized status SHA-256:
ad3bf3c97b378c9cdebb0b77d486cced544199750ad689060bd2a24f6a2210cb
merge JSON file SHA-256:
c7f5dea5811a8d2aa25d7bd3224b1fceae3fce73bb49fd4c8fe3f335e2e71c2f
merge payload SHA-256:
efcab50a9f0c67bb00aa6e11a53959205f4213f266072837f1f50fe87ef86459
~~~

As an independent regression rather than a proof premise, the earlier solver
without the final doubled-moment prune was run through indices 0--434. Its 16
archived logs cover the ordered ranges

~~~text
0:55, 55:71, 71:97, 97:123, 123:149, 149:175, 175:201, 201:227,
227:253, 253:279, 279:305, 305:331, 331:357, 357:383, 383:409,
409:435.
~~~

The first range is the completed prefix of an interrupted `0:113` owner.
All 435 verdicts and compact representatives agree byte-for-byte with the
full v2 prefix. The normalized payload (without a final newline) has SHA-256
`8b6b6277cb63561f744865ecc6aa7012dacc20be7d062c6b53d8670cfd7d75fd`;
the archived text file including its final newline has SHA-256
`86885254554eef4d616fc20aa937b1bde3fba35fb3cb2b8203d9454d5aad8d73`.
The comparison JSON file and its canonical payload have hashes
`32b1d64679b239fbc001eaf0182cb1978183880547a84c81541cee63a8483ce7`
and `b0fbbba6f26ef2f2579b6e72f8656c963091a46b99e21e5c5e30f55276f01890`.

This closes arbitrary compact labels only for the local
one-compact/six-AE row profile under the requirement that its odd,
degree-six, and degree-eight global forms are all zero. By (29h) that same
profile excludes the balanced local allocations for every $69\le t\le99$.
It does not exclude nonzero forms coordinated across directions or
unbalanced allocations, does not construct or exclude the Boolean
common-edge lift, and does not close residual (ii).

### Why no universal degree-six/eight identity can finish this gate

There is also a precise obstruction to extending the centered identity
(19) to arbitrary labels. For four triples $(a_i,b_i,0)$, map their sum to
the seven channels

\[
 (6,0),(6,1),(6,2),(8,0),(8,1),(8,2),(8,3).
\]

At

\[
 ((2,1),(3,2),(4,3),(5,4)),
\]

the Jacobian minor in the variables
$(b_1,a_2,b_2,a_3,b_3,a_4,b_4)$ is

~~~text
four compact atoms:
226534996574208000 = 2^28 * 3^9 * 5^3 * 7^3

four all-equal atoms:
220242357780480000 = 2^26 * 3^7 * 5^4 * 7^4
~~~

Both maps are therefore dominant on the distinct-label locus in every
characteristic at least 11. Every balanced branch-C hard row for $r\ge7$
has at least four compact atoms, and every opposite row has at least six
all-equal atoms; freezing all additional atoms preserves dominance.

This also answers the purely algebraic projective-interpolation question.
Over $\overline{\mathbf F}_p$, every signed row-type image contains a dense
open subset of the seven-dimensional value space. The finite intersection
of those opens is nonempty. Choose a common value $v$ and a linear form
$\ell$ over a finite extension with no zero on
$\mathbf P^1(\mathbf F_p)$. Scaling the base labels in direction $L$ by
$\ell(L)$ gives

\[
 C_{6,k}(L)=v_{6,k}\ell(L)^6,\qquad
 C_{8,k}(L)=v_{8,k}\ell(L)^8,
\]

which are genuine nonzero common binary forms. Thus there can be no
universal polynomial relation among these seven unrestricted atom channels,
and projective interpolation alone creates no algebraic contradiction.

The field extension is essential to the scope: this does **not** construct
$\mathbf F_p$ labels or $\mathbf F_p$ form coefficients. Finite-field
rational image holes, nonzero odd forms, the remaining moment hierarchy,
and the Boolean common-edge lift all remain live.

## 5. Remaining condition inside the antipodal construction

For one antipodal scale $u$, the joint even datum is

\[
 V(u)=(u^6,u^8).                                          \tag{32}
\]

The next task is a two-coordinate exact sumset. On branch B the same labels
must, for $q=6,8$, realize

\[
 \sum_{\rm hard}u^q=-F_q/2,\qquad
 A_q=-C_q/2,\qquad F_q=C_q+2B_q,                           \tag{33}
\]

where $(C_6,C_8)=V(c)$ comes from the omitted pair, $B$ is an exact
$b$-term sum of $V(u)$, and $A$ is an exact $(r-2)$-term sum. Branch C has
the analogous system without the omitted pair. This is the exact remaining
condition for this antipodal scale construction, not for arbitrary triangle
labels. The global moment-and-Boolean gate is broader.

Even solving (33) would not construct an integral preimage, much less

\[
 (z_0+\ker_{\mathbf Z}R)\cap\prod_e\{0,\tau_e\}.           \tag{34}
\]

The signed Boolean affine box remains open.

Replay:

~~~bash
cp evidence/p31_arbitrary_compact_fiber.cpp \
  /tmp/p31_arbitrary_compact_fiber_v2.cpp
g++ -O3 -std=c++20 /tmp/p31_arbitrary_compact_fiber_v2.cpp \
  -o /tmp/p31_arbitrary_compact_fiber_v2
/tmp/p31_arbitrary_compact_fiber_v2 0 450

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
  tests/test_compact_ray_moment_gate.py
~~~
