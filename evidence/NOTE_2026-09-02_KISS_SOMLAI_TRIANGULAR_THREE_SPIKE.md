# Kiss--Somlai triangular set gives an all-bad three-spike family

**Date:** 2026-09-02
**Status:** proved all-prime construction/barrier for odd primes $p\ge5$;
residual (ii) remains open

This note records an explicit realization of the all-bad three-coordinate
branch at defect $6p-12$. It is not a finite-prime inference and it is not an
exclusion: it proves that this branch is populated for every relevant prime.
The exact test module only guards the formulas at representative primes.

The geometric input is the triangular set introduced by Gergely Kiss and
Gábor Somlai in *Special directions on the finite affine plane*,
[arXiv:2109.13992](https://arxiv.org/abs/2109.13992), published in
*Designs, Codes and Cryptography* 92 (2024), 2587--2597,
[doi:10.1007/s10623-024-01404-y](https://doi.org/10.1007/s10623-024-01404-y).
Their Section 4 proves that the set below is equidistributed in all but three
directions; Section 5 gives the equivalent weighted-line/Fourier formulation.

## Paley and sign conventions

Put $V=\mathbb F_{p^2}$, write $\chi$ for its quadratic character extended by
$\chi(0)=0$, and define the finite Paley convolution operator

\[
 (Qf)(z)=\sum_{w\in V}\chi(z-w)f(w).
\]

The normalized symmetric Paley conference matrix is

\[
 C=\begin{pmatrix}0&\mathbf1^T\\ \mathbf1&Q\end{pmatrix}.
\tag{1}
\]

There are three different signs here, and conflating them reverses the
construction.

1. A geometric affine-line direction is $d\mathbb F_p$. Since every element
   of $\mathbb F_p^*$ is a square in $\mathbb F_{p^2}$, its spatial direction
   character $\chi(d)$ is well-defined. For every affine line
   $L=a+d\mathbb F_p$, the quadratic-character line sum gives

   \[
   Q\mathbf1_L=\chi(d)(p\mathbf1_L-\mathbf1).
   \tag{2}
   \]

   Thus a **square spatial direction** has $+p$ line-contrast eigensign.

2. With Fourier convention

   \[
   \widehat f(\xi)=\sum_z f(z)
      \exp(-2\pi i\operatorname{Tr}(\xi z)/p),
   \]

   the support of a line contrast lies on the trace annihilator of
   $d\mathbb F_p$, not on the spatial direction itself. If $\eta$ is the
   Legendre character of $\mathbb F_p$, a trace-zero direction has character
   $-\eta(-1)$. Hence the annihilator character is

   \[
   \chi(d^\perp)=-\eta(-1)\chi(d).
   \tag{3}
   \]

3. The quadratic Gauss sum over $\mathbb F_{p^2}$ for this convention is
   $-\eta(-1)p$. Multiplying it by (3) gives the Fourier multiplier

   \[
   (-\eta(-1)p)\chi(d^\perp)=p\chi(d),
   \tag{4}
   \]

   exactly as in the direct integer identity (2). For $p=1\pmod4$, a square
   spatial direction has a nonsquare Fourier annihilator and Gauss sign
   $-p$; its eigensign is nevertheless $+p$.

## The augmented triangular function

Use integer representatives $0,1,\ldots,p-1$ and set

\[
 S=\{(a,b)\in\mathbb F_p^2:0\le b<a\le p-1\}.
\]

Its three special geometric direction classes are horizontal, vertical, and
slope one:

\[
 D_0=\{\langle(1,0)\rangle,\langle(0,1)\rangle,
          \langle(1,1)\rangle\}.
\tag{5}
\]

Let

\[
 L=\{(0,b):b\in\mathbb F_p\},\qquad
 M=\{(a,p-2):a\in\mathbb F_p\},
\]

and define the integral function

\[
 f_0=\mathbf1_S+\mathbf1_L+\mathbf1_M.
\tag{6}
\]

There are exactly two overlaps:

\[
 u=(0,p-2)=L\cap M,\qquad
 v=(p-1,p-2)=S\cap M.
\tag{7}
\]

Indeed $S\cap L=\varnothing$, $S\cap M=\{v\}$, and
$L\cap M=\{u\}$. Consequently $f_0(u)=f_0(v)=2$, while every other
value of $f_0$ is zero or one. Also

\[
 \sum_z f_0(z)=\frac{p(p-1)}2+2p=\frac{p(p+3)}2.
\tag{8}
\]

For $p\ge5$, there are $(p+1)/2\ge3$ square projective directions. Choose
three distinct ones. Given representatives of the first two, their scalings
in $\mathbb F_p^*$ can be chosen so that their sum lies in the third
direction. Therefore there is a linear $T\in\mathrm{GL}(2,p)$ such that

\[
 T\langle(1,0)\rangle,\quad T\langle(0,1)\rangle,\quad
 T\langle(1,1)\rangle
\]

are all square spatial directions. This is an all-prime existence argument;
the audit module independently finds the lexicographically first such $T$
instead of assuming a closed formula.

Set $f=f_0\circ T^{-1}$. Kiss--Somlai equidistribution, or equivalently
their weighted-line description, puts the nonconstant part of $f$ in the sum
of the three line-contrast spaces in (5). After $T$, all three have
eigenvalue $+p$ by (2). Thus $Qf=pf-c\mathbf1$. Summing over $V$, using
$Q\mathbf1=0$, and applying (8) determines $c$ exactly:

\[
 \boxed{Qf=pf-\frac{p+3}{2}\mathbf1.}
\tag{9}
\]

## The integral eigenvector and its Boolean shadow

Define $y\in\mathbb Z^{p^2+1}$ by

\[
 y_\infty=3,\qquad y_z=2f(z)-1\quad(z\in V).
\tag{10}
\]

Equation (8) gives

\[
 \sum_{z\in V}y_z=2\sum_z f(z)-p^2=3p=p y_\infty.
\]

For a finite coordinate, (9) gives

\[
 (Cy)_z=3+Q(2f-\mathbf1)(z)
       =3+2pf(z)-(p+3)
       =p(2f(z)-1)=py_z.
\]

Hence

\[
 \boxed{Cy=py.}
\tag{11}
\]

By (7), the only coordinates of magnitude three are
$E=\{\infty,T(u),T(v)\}$, and all three values are $+3$. The three
conference edges on $E$ are positive. The two infinity edges are $+1$ by
normalization (1). Moreover,

\[
 T(v)-T(u)=T(-1,0),
\]

which lies in the square image of the horizontal direction; therefore
$\chi(T(v)-T(u))=+1$. Thus

\[
 C[E]=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix}.
\tag{12}
\]

Let $x\in\{\pm1\}^{p^2+1}$ be obtained from $y$ by replacing its three
$+3$ entries with $+1$, and let $z=-\mathbf1_E$. Then

\[
 y=x-2z,\qquad z_i=-x_i\quad(i\in E).
\tag{13}
\]

So this is the **all-bad** sign choice at all three spike coordinates. It is
also important that the endpoint triple is positive here. Proposition
15.639 uses the opposite projector convention when it calls its signed
triples negative: for the present $+p$ completion, the defect projector is
$(I-C/p)/2$, and (12) is the required positive signed triangle.

Finally, $\|y\|^2=p^2+1+24$, $\mathbf1_E^Ty=9$, and, by (12),
$\mathbf1_E^TC\mathbf1_E=6$. Since $x=y-2\mathbf1_E$ and $Cy=py$,

\[
 x^TCx
 =p\|y\|^2-4p\mathbf1_E^Ty
    +4\mathbf1_E^TC\mathbf1_E
 =p(p^2+1)-12p+24.
\]

Writing $\Phi=p(p^2+1)/2$ and $q_C(x)=x^TCx/2$, this proves

\[
 \boxed{\Phi-q_C(x)=6p-12.}
\tag{14}
\]

## Signed-projective normal form for the triple datum

The construction is not confined to the displayed positive triangle. Use
the canonical homogeneous representatives

\[
 r_\infty=(1,0),\qquad r_u=(u,1),
\]

so that $C_{PQ}=\chi(\det(r_P,r_Q))$. For
$g\in\mathrm{GL}(2,p^2)$, write

\[
 g r_P=\lambda_P r_{gP}.
\]

Taking quadratic characters of determinants gives the exact switching
factor

\[
 \boxed{C_{gP,gQ}=\chi(\det g)\chi(\lambda_P)
        \chi(\lambda_Q)C_{P,Q}.}
\tag{15}
\]

If $\det g$ is square, the signed permutation

\[
 (U_gw)_{gP}=\chi(\lambda_P)w_P
\tag{16}
\]

commutes with $C$. Multiplying (15) around a triangle cancels the three
vertex-switching factors and shows that its conference sign is multiplied by
$\chi(\det g)$.

Now $\mathrm{PGL}(2,p^2)$ is sharply three-transitive. Therefore the unique
projectivity carrying one ordered positive conference triangle to another
has square determinant and lies in $\mathrm{PSL}(2,p^2)$. The count makes
the same conclusion transparent. Since $C^3=p^2C$, the signed sum of all
unordered triangle products is $\operatorname{tr}(C^3)/6=0$, so exactly
half of the triangles are positive. Writing $q=p^2$,

\[
 6\cdot\frac12{q+1\choose3}
 =\frac{q(q^2-1)}2
 =|\mathrm{PSL}(2,q)|.
\tag{17}
\]

No nonidentity projectivity fixes three points, so the action on ordered
positive triangles is regular. Finally, replacing a matrix lift $g$ by
$\mu g$ for a nonsquare $\mu$ leaves the projectivity and determinant
character unchanged but reverses every $\chi(\lambda_P)$. Thus (16) also
matches either global sign of the signed triple.

Consequently all shell data

\[
 z=\pm e_i\pm e_j\pm e_k,\qquad
 C_{ab}z_az_b=+1\quad(a,b\in\{i,j,k\},\ a\ne b),
\]

form one signed-$\mathrm{PSL}(2,p^2)$ orbit. Transporting (10)--(13) gives
at least one triangular integral/Boolean completion for every such datum.
This is the precise normal-form statement. It does **not** say that every
possible Boolean completion of a fixed signed triple is unique or lies in
the same stabilizer orbit.

## Exact audit and consequence

The exact integer implementation is
**src/residual_kiss_somlai_three_spike.py**; its fail-when-wrong tests are in
**tests/test_residual_kiss_somlai_three_spike.py**. The tests use
$p=5,7,11,13$, including both congruence classes modulo four so that the
Fourier-annihilator/Gauss-sign cancellation in (3)--(4) is checked. They
verify the Kiss--Somlai line sums, choose $T$ deterministically, check (2),
(9), (11), every coordinate, all three spike edges, the all-bad signs,
(14), the exact switching factor (15), and the signed-projective orbit count.

This construction eliminates a possible shortcut: the all-bad three-spike
branch cannot be ruled out by conference eigenspace integrality, by the
signed-triangle shell conditions, or by an assertion that the endpoint has
no global all-prime realization. It does **not** construct a common residual
graph $H$, does **not** verify the global $A/B$ maximum inequalities, and
does **not** close residual (ii).
