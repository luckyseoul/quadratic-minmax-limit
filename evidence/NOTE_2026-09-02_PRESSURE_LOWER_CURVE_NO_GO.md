# Exact order-five obstruction to the product pressure lower curve

**Status:** proved all-temperature counterexample.  This note disproves one
specific proposed lower curve; it does **not** disprove the optimized
critical-pressure reformulation and does not settle the original limit.

Write

\[
 Q_A(x)=\sum_{i<j}a_{ij}x_ix_j,
 \qquad
 Z_A(c)=\mathbb E_x\cosh\!\left({cQ_A(x)\over\sqrt n}\right).
\]

The proposed assertion was

\[
 Z_A(c)\ \ge\
 \cosh\!\left(c\sqrt{1-{1\over n}}\right)^{n/2}
 \quad\hbox{for every complete symmetric signing }A.       \tag{1}
\]

It is false for every positive temperature already at order five.

## 1. One switching-normalized witness and its exact energy law

Let `A_5` be `-1` on the edges of a five-cycle and `+1` on its five
diagonals.  A Boolean state is the same as a cut, up to global sign.  For a
cut represented by a subset of size at most two:

* the empty cut and the five singleton cuts have energy zero;
* the five two-sets that are cycle edges have energy `-4`;
* the five two-sets that are cycle diagonals have energy `+4`.

Thus, on the sixteen cuts (or equivalently on the full cube after doubling
all multiplicities),

\[
 \#\{Q=-4,0,4\}=\{5,6,5\}.                                  \tag{2}
\]

In particular `Phi(A_5)=4`; this is the same optimal order-five signing
already used for the zero-error disk obstruction.

Put `t=2c/sqrt(5)` and `u=cosh(t)`.  Equation (2) gives

\[
 Z_{A_5}(c)
 ={5\over8}\cosh(2t)+{3\over8}
 ={5u^2-1\over4},                                           \tag{3}
\]

whereas the right side of (1) is

\[
 \cosh(2c/\sqrt5)^{5/2}=u^{5/2}.                            \tag{4}
\]

Writing `v=sqrt(u)`, one has `v>1` for every `c>0`, and

\[
 u^{5/2}-{5u^2-1\over4}
 ={4v^5-5v^4+1\over4}
 ={(v-1)^2(4v^3+3v^2+2v+1)\over4}>0.                       \tag{5}
\]

Hence (1) has the wrong direction for this witness at **every** `c>0`, not
only in a finite-temperature window.  Its zero-temperature slope also makes
the failure transparent: the left side has logarithmic slope `4/sqrt(5)`,
while (4) has slope `sqrt(5)`.

## 2. The first failed coefficient is the signed four-cycle invariant

Let `N=binom(n,2)` and let

\[
 C_4(A)=\sum_{\gamma}\prod_{e\in\gamma}a_e                 \tag{6}
\]

where the sum runs over unoriented simple four-cycles of `K_n` (three on
each four-set).  The surviving even multigraphs in the fourth Walsh moment
give

\[
 \mathbb E Q_A^4=3N^2-2N+24C_4(A).                          \tag{7}
\]

On the other hand, the fourth derivative at zero of the proposed product
curve is

\[
 \left.{d^4\over dc^4}
 \cosh\!\left(c\sqrt{1-{1\over n}}\right)^{n/2}
 \right|_{c=0}
 ={(n-1)^2(3n-4)\over4n}.                                   \tag{8}
\]

The left fourth derivative in (1) is `E Q_A^4/n^2`, so their difference is

\[
 {\mathbb E Q_A^4\over n^2}
 -{(n-1)^2(3n-4)\over4n}
 ={(n-1)(n-2)\over n}+{24C_4(A)\over n^2}.                  \tag{9}
\]

For `A_5`, `C_4(A_5)=-5`, and (9) equals `-12/5`.  Consequently

\[
 Z_{A_5}(c)
 -\cosh(2c/\sqrt5)^{5/2}
 =-{c^4\over10}+O(c^6).                                     \tag{10}
\]

Equivalently, using
`8C_4(A)=tr(A^4)-n(n-1)(2n-3)`, the proposed local inequality would require

\[
 \operatorname{tr}(A^4)\ge {n(n-1)(5n-7)\over3}.           \tag{11}
\]

That requirement points in the opposite direction from spectral flatness:
a symmetric conference matrix has
`tr(A^4)=n(n-1)^2` and violates (11) for every `n>2`.  This is consistent
with Proposition 6.9, whose Brascamp--Lieb estimate makes the same product
curve an **upper** bound on conference pressure.

## 3. What unconditional all-temperature lower curve remains

There is a universal lower curve, but its zero-temperature slope is exactly
the already-known `1/pi` floor rather than the upper construction scale
`1/2`.  This makes precise why it does not close the problem.

For `0<=r<=1`, put

\[
 I(r)={1+r\over2}\log(1+r)+{1-r\over2}\log(1-r).             \tag{12}
\]

Here the endpoint uses the convention `0 log 0=0`.

Choose a state `x_*` and sign `sigma_*` with
`sigma_* Q_A(x_*)=Phi(A)`.  In the augmented partition function

\[
 Z_A(c)=\mathbb E_{\sigma\in\{\pm1\},x}
          e^{\sigma cQ_A(x)/\sqrt n},                        \tag{13}
\]

fix `sigma=sigma_*` and independently bias every coordinate of `x` toward
the corresponding coordinate of `x_*`, with mean bias `r`.  The mean signed
energy is `r^2 Phi(A)` and the relative entropy from the uniform augmented
measure is `n I(r)+log 2`.  The Gibbs variational principle, together with
the elementary lower bound `Z_A(c)>=1`, therefore gives

\[
 {1\over n}\log Z_A(c)
 \ge \max\left(0,\sup_{0\le r\le1}
 \left\{c r^2{\Phi(A)\over n^{3/2}}-I(r)-{\log2\over n}\right\}\right).
                                                                    \tag{14}
\]

Using Proposition 5.2 yields, uniformly in `A`,

\[
 {1\over n}\log Z_A(c)
 \ge \max\left(0,\sup_{0\le r\le1}
 \left\{{c\sqrt{1-1/n}\over\pi}r^2-I(r)-{\log2\over n}\right\}\right).
                                                                    \tag{15}
\]

Apart from the harmless outer cutoff and `log(2)/n` term, the variational
curve acquires a nonzero stationary optimizer when
`2c sqrt(1-1/n)/pi>1`; that optimizer solves
`atanh(r)=2c sqrt(1-1/n)r/pi`.  As `c -> infinity`, its slope is
`sqrt(1-1/n)/pi`.  Thus (15) is a genuine universal pressure theorem, but it
only repackages the existing zero-temperature lower bound.  No matching
`1/pi` upper construction is known.

Conversely, suppose one proved an asymptotic universal pressure curve
`g(c)` with

\[
 \liminf_{n\to\infty}s_n(c)\ge g(c)
 \quad\hbox{and}\quad
 \lim_{c\to\infty}{g(c)\over c}={1\over2}.                  \tag{16}
\]

The entropy sandwich `s_n(c)/c<=alpha_n` would immediately imply
`liminf alpha_n>=1/2`; the known construction upper bound would then settle
the MO problem with limit `1/2`.  In other words, obtaining the desired
large-`c` slope is not a softer substitute for the missing lower bound: it
is that lower bound after the thermodynamic sandwich.

## 4. Exact scope

The counterexample retires the literal finite-`n` curve (1), including any
proof attempt based on determinant domination, entropy, hypercontractivity,
or Lee--Yang theory that concludes that same inequality.  It does not rule
out a lower bound with an order-dependent exponential loss, nor the exact
minimum-pressure gate of Proposition 6.10.  For such a softened bound to
force the conjectural zero-temperature constant `1/2`, its loss divided by
`c n` would have to vanish as `c -> infinity`; proving that remains at least
as strong as the missing asymptotic lower bound on `m_n`.
