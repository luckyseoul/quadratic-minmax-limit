# The complete third Paley-dual shell

Date: 2026-08-25. Proposition 15.636. This closes the one equality case
left open by Proposition 15.635 and classifies the third dual shell for
every odd prime `p>=11`. It does not control the later theta tail or prove
R1.

Put `m=(p-1)/2`. Proposition 15.635 reduced every possible extra vector at
scaled norm `2(p+1)` to one zero-common-sum profile. Up to sign, its
positive multiset has `m` elements, one repeated, while its negative
multiset has `m` distinct elements disjoint from the positive support. The
two multisets have equal power sums through degree `m-1`.

Let `A,B` be their monic root polynomials. Newton identities (all relevant
degrees are below `p`) show that `A-B` is constant. The root supports cover
`p-2` distinct field elements; write the two omitted elements as `u,v` and
the repeated root as `alpha`. Thus

\[
 A(X)B(X)=\frac{(X^p-X)(X-\alpha)}{(X-u)(X-v)}.       \tag{1}
\]

After the affine change `X=u+(v-u)x`, put
`lambda=(alpha-u)/(v-u)`. Since the repeated root is not omitted,
`lambda` is neither zero nor one. Equation (1) becomes

\[
 AB=(x-\lambda)\frac{x^p-x}{x(x-1)}
   =x^{p-1}+a(x^{p-2}+\cdots+x)-\lambda,\qquad
 a=1-\lambda.                                      \tag{2}
\]

Because `A-B` is constant, `T=(A+B)/2` is monic and

\[
 T(x)^2=x^{p-1}+a(x^{p-2}+\cdots+x)+b              \tag{3}
\]

for some `b`. Both `a` and `c=a-1=-lambda` are nonzero.

## The coefficient-gap contradiction

Reverse `T`:

\[
 U(y)=y^mT(1/y),\qquad U(0)=1.
\]

Equation (3) says that `U(y)^2` agrees through degree `p-2=2m-1`
with

\[
 \frac{1+cy}{1-y}=1+a(y+y^2+\cdots).
\]

The square root with constant term one is unique, so `U` agrees through
that degree with

\[
 R(y)=\sqrt{\frac{1+cy}{1-y}}.                     \tag{4}
\]

As `deg U=m`, the coefficients of `R` in degrees `m+1,...,p-2` vanish.
For degrees below `p`, binomial coefficients depend only on their exponent
modulo `p`. Since `m+1=1/2` and `m=-1/2` in `F_p`, those coefficients are
the coefficients of

\[
 K(y)=(1+cy)^{m+1}(1-y)^m.                         \tag{5}
\]

Therefore

\[
 K(y)=U(y)+q y^{p-1}+r y^p.                        \tag{6}
\]

The next-to-leading coefficient is

\[
 q=(-1)^m c^{m+1}\frac{a}{2c},
\]

so `q` is nonzero. Let `D^[k]` denote the `k`th Hasse derivative. For
`1<=j<=m`, (6), Lucas' congruence, and `binom(p-1,k)=(-1)^k` give

\[
 D^{[m+j]}K=q(-1)^{m+j}y^{m-j}.                    \tag{7}
\]

Evaluate (7) at the two roots `1` and `-1/c` of the factors in (5). The
local coefficients of order `m+j` are respectively

\[
 (-1)^m\binom{m+1}{j}c^j a^{m+1-j},\qquad
 (-1)^{j-1}\binom m{j-1}c^j a^{m+1-j}.
\]

Their ratio, compared with the monomial ratio in (7), yields

\[
 \boxed{c^{m-j}=-\frac1{2j}}\qquad(1\le j\le m).  \tag{8}
\]

Take `j=m-1` and `j=m-2`. In `F_p`, equation (8) says

\[
 c=\frac13,\qquad c^2=\frac15.
\]

But then `1/9=1/5`, so `p` divides four. This is impossible for every
odd `p>=11`. Hence the extra one-profile equality case does not exist.

Combining this contradiction with Proposition 15.635 proves

\[
 \boxed{\{u\in L^*:2p\|u\|^2=2(p+1)\}
   =\{\pm P(e_i+C_{ij}e_j):i\ne j\}.}              \tag{9}
\]

The complete signed count is `p^2(p^2+1)`, and the complete degree-four
harmonic shell operator is

\[
 -\frac{p^2+4p-3}{4(p^2+5)}\|W\|_F^2.             \tag{10}
\]

The exact NUKA audit found no exceptional profiles after `27,720`,
`180,180`, and `6,126,120` configurations at `p=11,13,17`; this checks the
reduction but is not used by the proof.

This is an ideal Prouhet--Tarry--Escott configuration over a prime field
with a prescribed near-partition and one repeated root. Rédei's lacunary
polynomial theory and Biró's two-value polynomial theorem are the closest
literature found. They explain the coefficient-gap setting but the searched
statements do not directly give (8) or the Paley-lattice shell theorem.
This is a search record, not an unqualified priority claim.
