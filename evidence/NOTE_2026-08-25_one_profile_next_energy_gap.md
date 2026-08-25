# The zero-common-sum gap after the complete third shell

Date: 2026-08-25. Proposition 15.637. Put

\[
k={p-1\over2},\qquad R=k+1={p+1\over2}.
\]

The first possible even profile energy after the complete third shell is
$E=p+3=2k+4$, corresponding to scaled dual norm $2(p+3)$. This note
excludes the entire zero-common-sum channel at that energy for every
$p\ge11$. The profile-balancing bound leaves only
$|t|\in\{2,p-1,p+1\}$ among nonzero common sums. Those cases remain, so
this does not classify the next shell or prove R1.

If $h<R$ profiles are active, the MDS/Newton mass bound and integer energy
give

\[
M\ge h(R-h),\qquad E\ge2M.                         \tag{1}
\]

For $k\ge5$, (1) leaves only $h=1$ and $h=R-1=k$; $h=R$ is the
separate full-direction case. We now exclude $h=1$.

For one profile, (1) gives $M\ge k$, while $E=2k+4$. If $a_s$ are
its integer entries, then

\[
E-2M=\sum_s |a_s|(|a_s|-1).                       \tag{2}
\]

The choice $M=k+2$ would require $p+3$ distinct nonzero entries and is
impossible. Thus (2) leaves exactly two multiplicity patterns.

## One doubled entry

If $M=k+1$, exactly one entry has magnitude two and all other nonzero
entries have magnitude one. The positive and negative root multisets each
have size $k+1$, are disjoint, and together cover every field element,
with one root $\alpha$ repeated. Let their monic root polynomials be
$A,B$. Equal power sums through degree $k-1$ imply that $A-B$ is
linear, and

\[
AB=(X^p-X)(X-\alpha).                              \tag{3}
\]

Put $T=(A+B)/2$, and reverse it as
$U(y)=y^{k+1}T(1/y)$. The reverse of (3), together with the square of the
linear difference, shows that $U^2$ agrees through degree $p-2$ with
$1-\alpha y$. The unique formal square root with constant term one is,
in characteristic $p$,

\[
\sqrt{1-\alpha y}=(1-\alpha y)^{k+1}\quad(\bmod y^p).
\]

Since $U$ has degree $k+1$, this forces

\[
T(X)=(X-\alpha)^{k+1}.
\]

Using Frobenius in (3) now gives

\[
{(A-B)^2\over4}=T^2-AB=(X-\alpha)^2.
\]

Hence $A-B=\pm2(X-\alpha)$, and both $A$ and $B$ are divisible by
$X-\alpha$, contrary to their disjoint root supports.

## Two doubled entries

If $M=k$, exactly two entries have magnitude two. The combined support
has size $p-3$. Normalize the three omitted roots to $0,1,\rho$, and
write the two repeated roots as $\alpha,\beta$. Now $A,B$ have degree
$k$, their power sums agree through degree $k-1$, and $A-B$ is a
nonzero constant. Define

\[
N(y)=(1-\alpha y)(1-\beta y),\qquad
D(y)=(1-y)(1-\rho y).                              \tag{4}
\]

Reversing the square identity shows that the degree-$k$ reverse root
agrees through degree $2k-1$ with

\[
S(y)=\sqrt{N(y)/D(y)}=\sum_{j\ge0}u_jy^j.
\]

Therefore

\[
u_{k+1}=u_{k+2}=\cdots=u_{2k-1}=0.                 \tag{5}
\]

Differentiate $DS^2=N$:

\[
2ND S'=(N'D-ND')S.                                 \tag{6}
\]

The polynomial on the right multiplying $S$ has degree at most two,
while the leading coefficient of $ND$ is
$\alpha\beta\rho\ne0$. In the coefficient of degree $j+3$, once
$u_{j+1},\ldots,u_{j+4}$ vanish, (6) reduces to

\[
2\alpha\beta\rho\,j\,u_j=0.                       \tag{7}
\]

For $k\ge5$, (5) reaches through $u_{k+4}$. Start (7) at $j=k$ and
descend to obtain $u_1=\cdots=u_k=0$. Comparing the coefficients of
degrees one and two in $DS^2=N$ then gives $N=D$. Thus
$\{\alpha,\beta\}=\{1,\rho\}$: the repeated roots are omitted roots,
again a contradiction.

Consequently no one-active zero profile has energy $p+3$. Combining this
with (1), a zero-sum candidate could now use only

\[
\boxed{h\in\{R-1,R\}.}                             \tag{8}
\]

## The full-direction branch

Let $q_d$ be the degree-$d$ binary form supplied by the profile glue. The
moments of every signed pair $\delta_a-\delta_b$ satisfy

\[
4q_1q_3-3q_2^2-q_1^4=0.                           \tag{9}
\]

For $h=R$, the energy budget makes exactly one profile have energy four;
the other $R-1=k\ge5$ profiles are signed pairs. The left side of (9) is
a binary quartic, so its five distinct zeros force it to vanish identically.
But on the exceptional profile
$\delta_a+\delta_b-\delta_c-\delta_d$ it equals

\[
-12(a-c)(a-d)(b-c)(b-d),                          \tag{10}
\]

which is nonzero because the four roots are distinct and $p\ge11$. Thus
$h=R$ is impossible.

## The one-inactive-direction branch

Now let $h=R-1=k$. The unique inactive direction is the root of the linear
form $q_1=L$. Since $q_2$ also vanishes there, write $q_2=LS$ with $S$
linear. A signed pair satisfies the lower-degree recurrence

\[
4q_3=L(3S^2+L^2).                                  \tag{11}
\]

The energy excess four is distributed either among two energy-four profiles
or one energy-six profile.

With two exceptions, the cubic difference in (11) vanishes at the $k-2$
ordinary profiles and the inactive direction: at least four points. It is
therefore identically zero, contradicting (10) at each exception.

With one energy-six exception, (11) is again an identity. The multiplicity
pattern $2\delta_a-\delta_b-\delta_c$ is impossible because its defect in
(9) is

\[
-12(a-b)^2(a-c)^2\ne0.                            \tag{12}
\]

The only remaining pattern has three distinct positive and three distinct
negative roots. Use the degree-four glue and write $q_4=LT$, where $T$ is
cubic. Every ordinary pair satisfies

\[
2T=S(S^2+L^2).                                     \tag{13}
\]

The cubic difference in (13) vanishes at the $k-1\ge4$ ordinary directions,
so it is identically zero. Equations (11) and (13) say that the exceptional
profile's first four power-sum differences equal those of one signed pair,
say $\delta_u-\delta_v$. After moving $u,v$ to the opposite sides, two
four-element multisets have equal first four power sums. Newton identities
make those multisets equal, contradicting the six disjoint roots.

This excludes $h=R-1$ as well. Therefore

\[
\boxed{\text{No zero-common-sum profile has energy }p+3
       \text{ for }p\ge11.}                       \tag{14}
\]

Since a common sum must be even at this scaled norm, the balancing
polynomial leaves exactly

\[
|t|\in\{0,2,p-1,p+1\}.
\]

The theorem above removes $t=0$; the other three magnitudes are the precise
remaining cases at this energy.

The exact CP-SAT audits at $p=11,13,17,19$ and the normalized square-shift
scan of 5,812,617 cases through $p=127$ independently found no one-profile
exception. They check the first branch; equations (3)--(14) are the proof.
