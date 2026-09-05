# Exact half-product minimizers have subcritical operator norm

2026-09-05. **All-orders theorem for actual global minimizers.** At every
fixed critical-scale inverse temperature, every sequence of exact global
half-product minimizers satisfies `||A_N||_op=o(N^(3/4))`. The proof uses
global optimality, a hereditary energy bound, and sparse pinning into an
actual complement Gibbs law. No census, solver, or numerical experiment
is used. This is a structural theorem, not a proof of cross-order
convergence or of the original MO limit.

The one prerequisite is the fully proved actual-Gibbs field-response
lemma in `evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md`, specifically its
moderate-coordinate version (13)--(14). Its statement is reproduced
below. At the time of this proof its SHA256 is
`46f6465c9a889dc485b9c24dac6f7fef8849d27271cc86df11b94ab732ed52dd`.

## 1. Definitions and theorem

For a complete signed host A of order N, let A be symmetric, have zero
diagonal, and have entries in {-1,1} off its diagonal. Set

\[
 Q_A(x)=\frac12x^TAx,\qquad
 \Phi(A)=\max_{x\in\{-1,1\}^N}|Q_A(x)|,
\]
\[
 Z_A(\beta)=\mathbb E_x e^{\beta Q_A(x)},\qquad
 a_A(\beta)=\frac12\bigl(\log Z_A(\beta)+\log Z_A(-\beta)\bigr),
\]

where every spin expectation without another indicated law is uniform.
The same definitions apply to principal restrictions and to arbitrary
real symmetric zero-diagonal interactions. The order-zero partition
function is one.

**Theorem.** Fix c>0. For each N let A_N be any exact global minimizer of
`a_A(c/sqrt(N))` among complete signed hosts of order N. Then

\[
 \boxed{\qquad \|A_N\|_{\rm op}=o(N^{3/4}).\qquad}             \tag{1}
\]

No uniqueness, conference structure, prescribed eigenvalue profile, or
assumption about the value of the unknown optimum is required.

## 2. Actual global optimality gives every induced subgraph a norm cap

Fix N, put beta=c/sqrt(N), and write A=A_N. If S is a vertex set of
size k and T its complement, global optimality gives

\[
 0\le a_A(\beta)-a_{A_T}(\beta)
 \le\left[\binom N2-\binom{N-k}2\right]\log\cosh\beta
 \le\frac{c^2k}{2}.                                         \tag{2}
\]

For the lower bound, fix the T spins and average over the uniform S
spins. Both their internal quadratic form and the cross term have mean
zero, so Jensen gives `Z_A(+/- beta)>=Z_(A_T)(+/- beta)`.

For the upper bound, retain A_T and independently randomize every edge
incident to S. For each such completed host A', exact optimality gives
`a_A(beta)<=a_(A')(beta)`. For each phase separately, the expected
partition function of A' over these random edges is exactly

\[
 Z_{A_T}(\pm\beta)(\cosh\beta)^e,
 \qquad e=\binom N2-\binom{N-k}2.
\]

Average the optimality inequality and apply concavity of log to each
phase. Finally use `e<=Nk` and `log cosh beta<=beta^2/2`.

In addition, independently reversing the entire S spin block leaves
both internal Hamiltonians unchanged and reverses the cross term.
Averaging this reversal in the product of the two internal Gibbs laws
replaces the cross exponential by a cosh, so for either phase

\[
 Z_A(\pm\beta)\ge Z_{A_S}(\pm\beta)Z_{A_T}(\pm\beta).
\]

Consequently

\[
 a_{A_S}(\beta)\le a_A(\beta)-a_{A_T}(\beta)
                       \le\frac{c^2|S|}{2}.                 \tag{3}
\]

Both phase log partition functions are nonnegative by Jensen. One
phase contains a maximizing configuration of absolute energy Phi(A_S),
with uniform probability at least 2^(-|S|). It follows that

\[
 \beta\Phi(A_S)-|S|\log2\le2a_{A_S}(\beta).
\]

Thus, writing C=c+(log 2)/c, we obtain the hereditary bound

\[
 \boxed{\quad \Phi(A_S)\le C\sqrt N\,|S|
                     \quad\hbox{for every }S\subseteq[N].\quad} \tag{4}
\]

The denominator sqrt(N), including for small induced sets, is the
original host scale. No optimality of the individual induced hosts is
being asserted.

## 3. The hereditary bound delocalizes every large eigenvector

If z is in [-1,1]^N, form independent ternary variables W_i with value
`sgn(z_i)` with probability |z_i| and value zero otherwise. Zero diagonal
and independence give `E Q_A(W)=Q_A(z)`. Conditional on the support of W,
its nonzero coordinates are a spin configuration on that induced host.
Equation (4) therefore implies

\[
 |Q_A(z)|\le\mathbb E|Q_A(W)|
       \le C\sqrt N\,\mathbb E|\operatorname{supp}W|
       = C\sqrt N\,\|z\|_1.                                 \tag{5}
\]

Let Av=lambda v, with ||v||_2=1 and M=||v||_infinity. Substitution of
z=v/M gives

\[
 |\lambda|\le2C\sqrt N\,M\|v\|_1.                           \tag{6}
\]

Every row of A has Euclidean norm sqrt(N-1). Hence

\[
 |\lambda|M\le\sqrt{N-1},\qquad
 \|v\|_1\ge\frac{|\lambda|^2}{2CN},\qquad
 |\lambda|\le\sqrt{2C}\,N^{3/4}.                             \tag{7}
\]

The latter two bounds follow from (6) and ||v||_1<=sqrt(N). In
particular, an eigenvalue at a fixed positive fraction of N^(3/4)
forces a linear number of moderately sized eigenvector coordinates.
Precisely, if

\[
 |\lambda|\ge\varepsilon N^{3/4},\qquad
 \delta=\frac{\varepsilon^2}{2C},\qquad Y_i=\sqrt N\,|v_i|,
\]

then `N^(-1) sum Y_i>=delta`, `N^(-1) sum Y_i^2=1`, and delta<=1 by
(7). Define

\[
 I=\{i:\delta/2\le Y_i\le4/\delta\},\qquad
 \eta=\delta^2/16.
\]

The coordinates below delta/2 contribute at most delta N/2 to sum Y_i.
Those above 4/delta contribute at most
`(delta/4) sum Y_i^2=delta N/4`. Thus the contribution on I is at
least delta N/4. As every Y_i on I is at most 4/delta,

\[
 |I|\ge\eta N.                                               \tag{8}
\]

## 4. A sparse pinned set generates a dense moderate field

Suppose, for a contradiction, that (1) fails. There are a fixed
epsilon>0 and an unbounded sequence of orders with a real eigenpair
(lambda,v) as in the preceding section and
`|lambda|>=epsilon N^(3/4)`. All following constants are independent of
N on this sequence. Put

\[
 \rho_0=\varepsilon/2,\qquad r=\rho_0N^{1/4}.
\]

Take independent ternary Z_i, equal to sgn(v_i) with probability
r|v_i| and zero otherwise. These probabilities are valid: the first
bound in (7) gives `r M<=rho_0/epsilon=1/2`. Then

\[
 \mathbb E Z=rv,\qquad
 \mathbb E|\operatorname{supp}Z|=r\|v\|_1
                                      \le\rho_0N^{3/4}.       \tag{9}
\]

Independence of the centered variables and the exact squared column
norm N-1 of A give

\[
 \begin{aligned}
 \mathbb E\|\beta A(Z-rv)\|_2^2
 &=\beta^2(N-1)\sum_i\operatorname{Var}(Z_i)\\
 &\le c^2r\|v\|_1\le c^2\rho_0N^{3/4}.
 \end{aligned}                                                \tag{10}
\]

Markov's inequality applied to (9) and (10), each with factor four,
shows that a deterministic realization exists with

\[
 k:=|\operatorname{supp}Z|\le4\rho_0N^{3/4},\qquad
 \|\beta A(Z-rv)\|_2^2\le4c^2\rho_0N^{3/4}.                  \tag{11}
\]

This is a proof of existence, not a sampling or computation step.
Fix such a realization, and let S be its support, T its complement,
and q=|T|. Thus k=o(N) and q/N tends to one. On T define the actual
pinning field

\[
 w=\beta A_{T,S}Z_S.
\]

Its comparison field on all coordinates is beta r lambda v, whose
coordinate magnitudes on I equal

\[
 c\rho_0\frac{|\lambda|}{N^{3/4}}Y_i.
\]

Using both eigenvalue bounds, these lie between the fixed positive
constants

\[
 a_0=\frac{c\rho_0\varepsilon\delta}{2},\qquad
 H_0=\frac{4c\rho_0\sqrt{2C}}{\delta}.                        \tag{12}
\]

Deleting the k coordinates in S loses o(N) of the coordinates of I.
By (11), the number of coordinates where the field error exceeds a_0/2
is at most `4||beta A(Z-rv)||_2^2/a_0^2=o(N)`. Therefore, for every
sufficiently large order in the subsequence,

\[
 \#\{i\in T:a_0/2\le|w_i|\le H_0+a_0/2\}
       \ge\frac\eta2 N\ge\frac\eta2 q.                       \tag{13}
\]

There is no assertion that the remaining coordinates of w are bounded.

## 5. Actual complement Gibbs response defeats the deletion budget

Here is the prerequisite in the precise form used. For a real symmetric
zero-diagonal interaction J on q spins, suppose Phi(J)<=bq. If a>0,
H>=a, and a fraction at least d>0 of the coordinates of a real field w
have `a<=|w_i|<=H`, with all other coordinates unrestricted, then

\[
 \Psi_J(w):=\log\mathbb E_{\mu_{J,0}}e^{w^TX}
       \ge K(b,a,H,d)q,                                      \tag{14}
\]
\[
 K(b,a,H,d)=\frac{a^2d}{4}
    \exp\left[-2\sqrt{2H^2+\frac{8B(1+4B)}{d^2}}\right]>0,
 \quad B=1+\frac{4\pi b}{\log(1+\sqrt2)}.
\]

The law mu_(J,0) is the actual zero-field Gibbs law proportional to
`exp(Q_J(x))`. The lemma applies equally to J and -J. In particular,
it does not replace either law by independent uniform spins or clip
any coordinate of the actual field.

Apply (4) to S=T to obtain

\[
 \Phi(\pm\beta A_T)=\beta\Phi(A_T)\le cCq.
\]

Thus (14) applies to both complement phases with the fixed parameters

\[
 b=cC,\qquad a=a_0/2,\qquad H=H_0+a_0/2,\qquad d=\eta/2.
\]

Writing their positive response constant as K, (13)--(14) give

\[
 \Psi_{\beta A_T}(w)\ge Kq,
 \qquad \Psi_{-\beta A_T}(w)\ge Kq.                           \tag{15}
\]

Now retain in the positive phase partition function only configurations
with X_S=Z_S. Retain in the negative phase only configurations with
X_S=-Z_S. Each pinning event has uniform probability 2^(-k). Both
complements see the SAME external field w: the reversal of the pinned
spins in the negative phase cancels that phase's minus sign in the
cross term. The internal energy of S is unchanged by its reversal,
so its two signed contributions cancel when the log partition
functions are averaged. Explicitly,

\[
 Z_A(\beta)\ge2^{-k}e^{\beta Q_{A_S}(Z_S)}
                  Z_{A_T}(\beta)e^{\Psi_{\beta A_T}(w)},
\]
\[
 Z_A(-\beta)\ge2^{-k}e^{-\beta Q_{A_S}(Z_S)}
                  Z_{A_T}(-\beta)e^{\Psi_{-\beta A_T}(w)}.
\]

Hence (15) yields

\[
 a_A(\beta)-a_{A_T}(\beta)
 \ge-k\log2+\frac12\left[\Psi_{\beta A_T}(w)
                          +\Psi_{-\beta A_T}(w)\right]
 \ge Kq-k\log2.                                             \tag{16}
\]

The right side is `K N-o(N)`, whereas exact optimality's upper bound
(2) is `c^2 k/2=o(N)`. This contradiction proves (1).

## 6. What is and is not established

The theorem applies to every choice of exact half-product minimizers,
at each fixed c>0. Global optimality is used quantitatively in (2)
and then again after pinning; a generic norm cap alone would not give
this proof. The first step also proves the hereditary cap (4), but
that cap alone is not asserted to imply the little-o conclusion.

The theorem does not say `||A_N||_op=O(sqrt(N))`, does not identify
minimizers, and does not prove a useful rate in the little-o bound.
It does not assert the conclusion for arbitrary leading near-minimizers
or for minimizers of the Boolean absolute norm Phi.

It does not prove a cross-order pressure comparison. In particular,
the zero-temperature slope of the half-product pressure is HALF THE
ENERGY WIDTH, `(max Q_A-min Q_A)/2`, and is not in general Phi(A).
No passage from this structural theorem to the original MO convergence
problem is being claimed.
