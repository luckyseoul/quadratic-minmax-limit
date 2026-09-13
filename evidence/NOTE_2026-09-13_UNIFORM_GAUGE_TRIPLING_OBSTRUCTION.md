# Uniform gauge averaging overshoots the tripling target

2026-09-13. An all-orders obstruction to Section 4.2 of
`NOTE_2026-09-02_UPSTREAM_RELATIVE_GAUGE_BRIDGE.md`.

The uniform reverse-KL criterion introduced there is stronger than selection
of one good gauge. In fact its order-`n` plus order-`2n` hypothesis is false
for all sufficiently large `n`, for every choice of optimal internal blocks
and every cross seed. The obstruction also applies to its positive
fractional-moment sufficient bounds. The selected-gauge certificate in
Section 4.1 and the original convergence question remain open.

## 1. An exact lower bound for the uniform average over gauges

Let `A` and `B` be complete symmetric zero-diagonal signings of orders
`k` and `m`, and let `C` be a `k`-by-`m` sign matrix. In the absorbed
coordinates of the bridge note, write

```text
D = D_alpha C D_beta,
Y_g = [[A,D],[D^T,tau B]].
```

Average over the uniform projective relative gauge, including `tau=+/-1`.
Equivalently, `alpha` and `beta` may be independent full sign vectors:
negating either whole vector changes `D` to `-D`, which leaves the Boolean
norm and symmetric partition function invariant. Thus each projective
gauge is represented with the same multiplicity.

For independent fair signs define

```text
a_j = E |epsilon_1+...+epsilon_j|,       a_0=0.
```

For every integer `0<=s<=k`,

```text
E_g Phi(Y_g) >= Phi(B) + s a_m + (k-s) a_(m+s).          (1)
```

No independence of the matrix entries or assumptions about their spectra
are used.

Choose `y_*` and `sigma in {+1,-1}` such that
`Q_B(y_*)=sigma Phi(B)`. For every `x`, the exact diamond identity gives

```text
Phi(Y_g) >= |Q_A(x)+tau Q_B(y_*)| + |x^T D y_*|
         >= Phi(B) + t Q_A(x) + x^T D y_*,
t = tau sigma.                                         (2)
```

Put `eta=D_beta y_*`, `h=C eta`, `H_i=|h_i|`, and
`xi_i=alpha_i sign(h_i)`, using `sign(0)=1` for this change of variables.
Conditional on `eta`, the `xi_i` are independent fair signs, independent
of the fair sign `t`, and `(D y_*)_i=xi_i H_i`.

Fix a subset `S` of the free block of size `s`, and let `T` be its
complement. Set

```text
x_i = xi_i                                  for i in S,
F_i = sum_(j in S) A_ij xi_j                 for i in T,
x_i = sign(xi_i H_i + t F_i)                 for i in T. (3)
```

At zero in the last line choose an independent fair sign. These are actual
Boolean spins; zero is used only for the mean of the tie coin below.
Conditional on `eta` and `xi_S`, let `sgn_0(0)=0` and put

```text
psi_i = [sgn_0(H_i+F_i)+sgn_0(-H_i+F_i)]/2.
```

Then `E[x_i | t,eta,xi_S]=t psi_i` for `i in T`. Distinct spins in `T`
are conditionally independent given these variables, because their remaining
signs `xi_i` and tie coins are independent. It follows that

```text
E[t x_i x_j | eta,xi_S]=0              for distinct i,j in T. (4)
```

This cancels the entire quadratic self-interaction of the updated block
exactly under the `t` average. The term from edges inside `S` also has
zero `t` average. The cross edges from `S` to `T` combine with the field
at each updated coordinate to give

```text
E[(xi_i H_i+t F_i)x_i | eta,xi_S]
  = (|H_i+F_i|+|H_i-F_i|)/2
  = max(H_i,|F_i|).                                     (5)
```

Thus the conditional expectation of the last expression in (2) is

```text
Phi(B) + sum_(i in S) H_i + sum_(i in T) max(H_i,|F_i|). (6)
```

Each `H_i` has the law `|epsilon_1+...+epsilon_m|`. Conditional on `eta`,
each `F_i` is a sum of `s` independent fair signs, with a law independent
of `eta`. Therefore `H_i` and `F_i` are independent for each `i`.
For independent symmetric random variables `U,V`,

```text
E max(|U|,|V|) = E (|U+V|+|U-V|)/2 = E |U+V|.
```

Apply this with independent sign sums of lengths `m` and `s` to obtain
`E max(H_i,|F_i|)=a_(m+s)`. Averaging (6) proves (1). Correlations between
different coordinates of `h` or of `F` do not enter this calculation.

## 2. The tripling constant has a strict gap

Now take `k=n`, `m=2n`, and optimal internal blocks, and put

```text
T_n = (m_n^(2/3)+m_(2n)^(2/3))^(3/2),
s_n = floor(n/2),
J_n = s_n a_(2n) + (n-s_n) a_(2n+s_n).
```

Equation (1) implies, uniformly in all such blocks and every cross seed,

```text
E_g Phi(Y_g) >= m_(2n) + J_n.                           (7)
```

For `j>=1` the exact sign-sum formula is

```text
a_j = j binom(j-1,floor((j-1)/2))/2^(j-1).
```

It follows either by conditioning one summand in
`|S_j|=sum_i epsilon_i sgn_0(S_j)`, or by pairing consecutive binomial
terms. Stirling's formula consequently gives
`a_j/sqrt(j) -> sqrt(2/pi)`. Hence

```text
J_n/n^(3/2) -> kappa = (2+sqrt(5))/(2 sqrt(pi)).          (8)
```

The existing all-orders upper bound in CORE, Section 6, is
`limsup alpha_n<=1/2`. For nonnegative `u,v`, the function
`(u^(2/3)+v^(2/3))^(3/2)-v` is increasing in both arguments.
For positive `v` its derivative in `v` is
`sqrt(1+(u/v)^(2/3))-1>=0`; the boundary follows by continuity.
The derivative in `u` is nonnegative as well. Apply the CORE bound to
both `m_n` and `m_(2n)` to conclude

```text
limsup_n (T_n-m_(2n))/n^(3/2)
 <= kappa_0 = (3 sqrt(3)-2 sqrt(2))/2.                  (9)
```

There is a rigorous constant gap

```text
kappa-kappa_0 > 1/100.                                 (10)
```

Here is an entirely rational certificate of (10). The elementary integral

```text
0 < integral_0^1 x^4(1-x)^4/(1+x^2) dx = 22/7-pi
```

gives `pi<22/7`. Polynomial division uses the quotient
`x^6-4x^5+5x^4-4x^2+4` and remainder `-4`. The following root bounds
are verified by squaring positive rationals:

```text
sqrt(pi) < 17729/10000,
sqrt(5)  > 2236/1000,
sqrt(3)  < 17321/10000,
sqrt(2)  > 14142/10000.
```

They give

```text
kappa > 21180/17729,
kappa_0 < 23679/20000,
kappa-kappa_0 > 3795009/354580000 > 1/100.
```

Combining (7)--(10), for every sufficiently large `n` and every choice
of optimal blocks and cross seed,

```text
E_g Phi(Y_g) >= T_n + n^(3/2)/200.                     (11)
```

The larger limiting gap in (10) leaves room for the asymptotic errors in
(8)--(9). No convergence of `alpha_n`, value `L=1/2`, or optimizer
classification is assumed. The known upper bound is reused; no matrix
family construction or scan is reopened.

## 3. The recently proposed reverse-KL criterion cannot hold

Use the bridge note's normalized positive partition functions, and write
`F_g=log Z_g`, `S_beta=log Z_A+log Z_B+log Z_C`. Its density is
`rho_g=Z_g/(Z_A Z_B Z_C)`. The reverse KL is exactly

```text
I_beta = -E_g log rho_g = S_beta-E_g F_g.                (12)
```

Consequently the proposed inequality `I_beta>=S_beta-beta T_n` is
equivalent to `E_g F_g<=beta T_n`. It constrains the uniform average
pressure, whereas existence of one low-pressure gauge only constrains
the minimum.

The pressure/norm sandwich and (11), with `N=3n`, imply

```text
E_g F_g - beta T_n
 >= beta n^(3/2)/200 - 3n log 2.                       (13)
```

At the proposed `beta=log(n+1)^2/sqrt(3n)`, this is

```text
n [log(n+1)^2/(200 sqrt(3)) - 3 log 2],
```

which is positive for all sufficiently large `n` and grows faster than
`n`. Therefore Section 4.2's tripling hypothesis is false for every cross
seed and every pair of optimal internal blocks once `n` is sufficiently
large. An additional fixed `O(N)` slack in that pressure inequality cannot
repair it. The formal sufficient implication remains valid, but its
hypothesis is unattainable in this setting and is retired as a target.

The same obstruction applies to the positive fractional-moment certificate

```text
(1/q) log E_g rho_g^q <= beta T_n-S_beta,       q>0,
```

including `q` depending on `n`: adding `S_beta` to its left side gives
`(1/q)log E_g exp(q F_g)>=E_g F_g` by Jensen. Equation (13) excludes the
required upper bound. This does not address negative moments, nonuniform
gauge laws, or a carefully selected individual gauge. In particular it
does not contradict the original one-good-gauge composition target.

## 4. Verification and resulting scope

Root derived and reviewed the analytic proof (1)--(13). NUKA ran one
new serial exact-arithmetic check, recorded in
`uniform_gauge_tripling_obstruction_20260913/`. On one fixed `4+4` input,
all 64 conditional contexts had exact cancellation in (4); 24 had nonzero
self-interaction at fixed `tau`, so the cancellation check was nonvacuous.
The 832 strategy outcomes were valid actual-diamond lower bounds, and the
exact mean score was `43/4`, agreeing with (1). Rational checks verified
the polynomial identity and all constant bounds used in (10).

The finite fixture corroborates the implementation of the proof's strategy;
it is not an extrapolation to arbitrary orders or an independent human
review. The all-orders obstruction follows from the analytic calculation
and the existing CORE upper bound. This result retires the uniform
reverse-KL and stated positive-moment selection targets, not the broader
relative-gauge method or the original convergence problem.
