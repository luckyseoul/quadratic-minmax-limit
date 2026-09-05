# Affine-parity cube floor and a closed middle-boundary branch

Date: 2026-09-04.

Classification: **proved infinite-family local branch exclusion**. For every
odd $p\geq29$, a nonnegative integral quadratic $A$ on
$J(p,(p+1)/2)$ of mass

\[
2p\,\mathbb E A\in\{2p+4,\,2p+6\}
\]

cannot have an affine parity whose **even support representative** has
size $6\leq b\leq p-5$. The result is independent of the parity phase and
does not require primality. It applies in particular to both low masses at
the first uncovered residual layer $k=5p+5$, but does not alone close that
whole layer or a global acceptance predicate.

This is a new use of integrality and affine parity on every paired cube.
It does not extend the strict mass interval in Proposition 15.774, classify
the endpoint $2p-10$, or run a cube, graph, prime, or equality-family census.
No value of the actual parallel count or transported signed total is assumed.

## 1. Dimension-free cube parity mean

Let $f:\{0,1\}^N\to\mathbb Z_{\geq0}$ have degree at most two and parity

\[
f(x)\equiv\eta+\sum_{i\in D}x_i\pmod2,\qquad d=|D|.
\]

Then

\[
\boxed{\mathbb E f\geq\frac{\lceil d/2\rceil}{2}
                    =\frac{d+1_{\{d\text{ odd}\}}}{4}.}
\tag{1}
\]

Integer values make all multilinear coefficients integral. The affine
parity makes every quadratic coefficient even. Write

\[
f(x)=a+\sum_i l_i x_i+2\sum_{i<j}b_{ij}x_ix_j,
\qquad a,l_i,b_{ij}\in\mathbb Z.
\]

With $z_i=2x_i-1$, homogenize the linear terms using an extra sign $z_0$:

\[
\widetilde f(z_0,z)
=c+\frac12\sum_{0\leq i<j\leq N}w_{ij}z_i z_j,
\quad
w_{ij}=b_{ij}\ (i,j>0),\quad
w_{0i}=l_i+\sum_{j\ne i}b_{ij},
\tag{2}
\]

where $c=\mathbb E f$. The values of $\widetilde f$ are those of
$f((1+z_0z_i)/2)_i$, so they are nonnegative. At an ordinary vertex $i$,
the weighted degree is congruent to $l_i$ modulo two. At vertex zero,
the weighted degree is congruent to $\sum_i l_i=d$ modulo two.
There are therefore exactly $d+(d\bmod2)$ odd weighted degrees.

Choose a global minimum of $q(z)=\sum_{i<j}w_{ij}z_i z_j$.
Every local field $z_i\sum_j w_{ij}z_j$ is nonpositive; otherwise
flipping that sign would decrease $q$.
Each odd-degree local field is an odd integer and hence at most $-1$.
Summing local fields gives

\[
2q_{\min}\leq-d-(d\bmod2).
\]

Nonnegativity of (2) now gives $c\geq-q_{\min}/2$, proving (1).
The homogenizing vertex is essential: dropping its odd degree would
lose the improvement used below.

## 2. Paired cubes on the middle slice

Let $p\geq5$ be odd, $m=(p+1)/2$, $q=m-1$, and suppose

\[
A:J(p,m)\to\mathbb Z_{\geq0},\qquad \deg A\leq2,\qquad
A(X)\equiv\eta+|X\cap B|\pmod2,
\]

where $b=|B|$ is **even**. This convention entails no loss for an arbitrary
affine parity: if its original support has odd size, replace it by its
complement and replace $\eta$ by $\eta+m$ modulo two. Since $p$ is odd,
the new support is even and represents the same parity on this slice.

Fix a slice point $X$, and put $a=|X\cap B|$. Select uniformly one
unmatched vertex of $X$, and match the other $q$ vertices bijectively
to $X^c$, uniformly over all bijections. The choices of one endpoint
from each pair form a $q$-cube through $X$. A cube coordinate is
parity-active exactly when its pair meets $B$ once. Let $d$ be the
number of these coordinates.

Counting the two types of crossing pairs gives

\[
\mathbb E d
=\frac{a(q-b+a)+(m-a)(b-a)}m
=b-\frac{a(2b+1-2a)}m.
\tag{3}
\]

Since $b$ is even, $d$ is odd exactly when the unmatched vertex belongs
to $B$. Thus

\[
\mathbb P(d\text{ odd})=\frac am.
\tag{4}
\]

Combining (1), (3), and (4), the average mean of these restricted cubes
is at least

\[
\begin{aligned}
\mathbb E_{\rm cubes}\mathbb E_{\rm cube}A
&\geq\frac14\left(b-\frac{2a(b-a)}m\right)\\
&\geq\frac{b(p+1-b)}{4(p+1)}.
\end{aligned}
\tag{5}
\]

The last inequality is $a(b-a)\leq b^2/4$. Its exact excess is
$(2a-b)^2/[4(p+1)]$. In particular there is no assumption that a
minimum point intersects $B$ in exactly half its vertices.

The existing averaged-cube operator from Proposition 15.688 is

\[
TA(X)=\frac{A(X)+p\,\mathbb E A}{p+1}.
\tag{6}
\]

For completeness, (6) can be checked on $1,x_i,x_ix_j$. The averaged
cube means of $x_i$ are $(m+1)/(2m)$ inside $X$ and $1/2$ outside.
For $x_ix_j$, the mean is $(m+2)/(4m)$ when both indices are in $X$,
and $1/4$ otherwise. These equal (6), using the corresponding
middle-slice means $m/p$ and $m(m-1)/[p(p-1)]$.

Apply (5)--(6) at a point minimizing $A$. This proves the general bound

\[
\boxed{b(p+1-b)\leq4\bigl(\min A+p\,\mathbb E A\bigr).}
\tag{7}
\]

## 3. The two first-uncovered masses

Suppose $2p\,\mathbb E A=2p+s$ with $s=4$ or $6$. In the ranges below,
$\mathbb E A=1+s/(2p)<2$. Integrality implies $\min A\leq1$, so (7) gives

\[
b(p+1-b)\leq4p+2s+4.
\tag{8}
\]

For every even $6\leq b\leq p-5$, concavity and the two endpoints give

\[
b(p+1-b)\geq6(p-5).
\tag{9}
\]

The difference of the lower and upper bounds is exactly

\[
6(p-5)-(4p+2s+4)=2(p-17-s).
\tag{10}
\]

Consequently this whole boundary interval is excluded at mass $2p+4$
for $p>21$ (in particular every prime $p\geq23$), and at mass $2p+6$
for $p>23$ (in particular every prime $p\geq29$). Both masses are therefore
excluded on the whole middle-boundary interval for every odd $p\geq29$,
in both parity phases.

Equality at the numerical thresholds is not excluded by this argument.
The remaining even boundary candidates from this theorem alone are
$0,2,4,p-3,p-1$. Other established local arguments may remove some of
them, but are not silently included in this claim.

## 4. Scope and exact regression

The new source is
[e1_gmin_m4_affine_parity_cube_floor.py](../src/e1_gmin_m4_affine_parity_cube_floor.py),
with tests in
[test_affine_parity_cube_floor.py](../tests/test_affine_parity_cube_floor.py).

The source uses exact integers and Fraction values. Tests cover the
homogenizing odd-degree term, the monomial operator identity, the exact
intersection excess, odd-support phase normalization, strict endpoint
preservation, and failure injection for the cube/operator dependencies.
They verify algebra and scope, not a finite classification.

No first-layer closure or all-size closure flag is set. In particular,
the scalar survivor at $H=5p+6$ does not by itself force $P=5$ or
$T=\pm1$, and neither assumption is used here.
