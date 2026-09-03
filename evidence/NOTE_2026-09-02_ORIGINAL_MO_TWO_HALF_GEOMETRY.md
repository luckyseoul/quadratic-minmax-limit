# Original MO problem: exact two-half geometry calculator

**Status:** proved algebraic reduction plus exhaustive finite certificates at
`n=5,6,7,8`.  This does **not** prove convergence of `m_n/n^(3/2)`.

## 1. The two halves

Fix a symmetric signing `A`, put `M=Phi(A)`, and let `R` be a skew signing.
For any Boolean pair `x,y`, set

\[
 u={x+y\over2},\qquad v={x-y\over2}.
\]

Thus `u,v` are signed indicators of complementary coordinate sets.  With

\[
 I=Q_A(u)+Q_A(v),\quad X=u^TAv,\quad C=-u^TRv,
\]

direct expansion gives

\[
 Q_A(x)=I+X,\qquad Q_A(y)=I-X,\qquad x^TRy=2C.       \tag{1}
\]

The two full Boolean endpoint states imply the exact hereditary diamond

\[
 |I|+|X|=\max(|I+X|,|I-X|)\le M.                   \tag{2}
\]

Define the internal slack and skew interaction

\[
 D=2M-|Q_A(x)+Q_A(y)|=2(M-|I|),\qquad W=|x^TRy|=2|C|.
\]

Then the equal-endpoint doubling target is not merely analogous to a planar
picture; it is exactly the diagonal-envelope inequality

\[
 \boxed{\quad \max_{x,y}(W-D)\le2(\sqrt2-1)M+o_{\rm Dini}(n^{3/2}).\quad} \tag{3}
\]

Equivalently, if

\[
 B(A,R)=\max_{x,y}\{|Q_A(x)+Q_A(y)|+|x^TRy|\},
\]

then

\[
 B(A,R)-2M=\max_{x,y}(W-D).                         \tag{4}
\]

Writing `H=M-|I|-|X|>=0`, the half-scale form of (3) is

\[
 |C|-|X|-H\le(\sqrt2-1)M+o_{\rm Dini}(n^{3/2})/2.  \tag{5}
\]

This identifies the missing correlation precisely: a large oriented cut
`|C|` must be paid for by endpoint cross-energy `|X|` or unused endpoint
slack `H`.  A separate norm bound on `R` throws away exactly those two
payments.

The same statement has a useful endpoint-only form.  The identity
`|a+b|+|a-b|=2 max(|a|,|b|)` gives, with
`epsilon(x,y)=M-max(|Q_A(x)|,|Q_A(y)|)`,

\[
 \boxed{D=|Q_A(x)-Q_A(y)|+2\epsilon(x,y),}          \tag{5a}
\]

and hence

\[
 B(A,R)-2M
 =\max_{x,y}\bigl(|x^TRy|-|Q_A(x)-Q_A(y)|-2\epsilon(x,y)\bigr). \tag{5b}
\]

Thus both endpoint-energy separation and individual distance from the norm
boundary are exact payments.  Dropping either changes the problem at leading
order.

There is an equally exact directed-half interpretation.  For
`U={i:x_i=-y_i}`, split the `y`-energy on the cut into `F` on arcs from `U`
to `U^c` and `G` on the reverse arcs.  Then

\[
 Q_A(y)=I+F+G,\quad Q_A(x)=I-F-G,\quad x^TRy=2(G-F),
\]

so

\[
 |x^TRy|-|Q_A(x)-Q_A(y)|
 =\begin{cases}
 4\min(|F|,|G|),&FG<0,\\
 -4\min(|F|,|G|),&FG\ge0.
 \end{cases}                                                \tag{5c}
\]

Only cancellation between oppositely signed directed halves can be
dangerous.  For `FG<0`, the complete multiplier-two condition is exactly

\[
 2\min(|F|,|G|)
 \le(\sqrt2-1)M+\epsilon(x,y)                               \tag{5d}
\]

up to the required Dini error; pairs with `FG>=0` satisfy it automatically.
This is the sharp sign-sensitive form of the plotted geometry.

The same boundary can be read directly from the energy layers.  Put
`e(z)=M-|Q_A(z)|`.  Since
`epsilon(x,y)=min(e(x),e(y))`, the free payment is

\[
 |Q_A(x)-Q_A(y)|+2\epsilon(x,y)
 =\begin{cases}
 e(x)+e(y),&Q_A(x)Q_A(y)\ge0,\\
 2M-|e(x)-e(y)|,&Q_A(x)Q_A(y)\le0.
 \end{cases}                                                \tag{5e}
\]

Thus the narrowest constraints are pairs in the same positive or negative
near-maximizer layer.  At the critical half-excess
`t=(sqrt(2)-1)M`, success is exactly

\[
 |x^TRy|\le2t+e(x)+e(y)                                    \tag{5f}
\]

on every same-sign pair, together with the second line of (5e) on
opposite-sign pairs.  Since every skew signing has
`max_(x,y)|x^TRy|>=n mu_(n-1)`, a skew-extremizing same-sign pair must obey

\[
 e(x)+e(y)\ge
 \bigl(\sqrt{2/\pi}-(\sqrt2-1)-o(1)\bigr)n^{3/2}
 =(0.383670998\ldots-o(1))n^{3/2}.                          \tag{5g}
\]

The exact nonlinear covering dual and the proof that LP, SDP, and every
fixed-degree moment relaxation miss this layer constraint are recorded in
`NOTE_2026-09-02_BIVECTOR_ENERGY_LAYER_MINIMAX.md`.

## 2. The same geometry as one complex quadratic form

There is an exact sphere representation, not merely a change of variables.
Put

\[
 w={x+iy\over\sqrt2}\in\{(\pm1\pm i)/\sqrt2\}^n,
 \qquad H=A+iR.
\]

The matrix `H` is Hermitian, and direct expansion gives

\[
 w^*Hw=Q_A(x)+Q_A(y)-x^TRy.                        \tag{6}
\]

Replacing `y` by `-y` reverses only the last term, so

\[
 \max_w|w^*Hw|
 =\max_{x,y}\{|Q_A(x)+Q_A(y)|+|x^TRy|\}=B(A,R).   \tag{7}
\]

Consequently the desired diamond follows from the deterministic spectral
condition

\[
 \|A+iR\|_{\rm op}\le {2\sqrt2 M\over n}
\]

(and the corresponding condition with a Dini error).  This spectral
sufficient condition cannot be the general argument: Frobenius norm forces
`||A+iR||op>=sqrt(2(n-1))`, so it can be sharp only at the conference
scale `M>=n sqrt(n-1)/2`.  It does expose a tempting algebraic special case.
If

\[
 A^2=(n-1)I,\quad R^2=-(n-1)I,\quad AR+RA=0,
\]

then `(A+iR)^2=2(n-1)I`; whenever `M=n sqrt(n-1)/2`, this would prove the
zero-error doubling diamond exactly.  But the last equality is impossible
for sign-valued `A,R` at every `n>=3`.  Indeed, putting
`s_ij=A_ij R_ij` and `p_i=product_(k!=i) s_ik`, exact anticommutation would
force `p_i p_j=-1` for every pair, which is inconsistent on any three
indices.  The full proof and a quantitative parity floor are in the
[fourth-phase/Clifford audit](NOTE_2026-09-02_FOURTH_PHASE_CLIFFORD_AUDIT.md).
Approximate flattening on the discrete fourth-phase torus is not excluded.

### Exact directed-half-cut form

The fourth-phase picture has an exact real combinatorial slice.  Put
`S=A circ R`, viewed as a tournament, and for `U subset [n]` let
`F_S(U)` be the cut edges oriented from `U` to its complement.  If
`A^F` means reversing the signs on `F`, then Proposition 6.5a proves

\[
 {1\over2}B(A,R)
 =\max_{U\subseteq[n]}\Phi\left(A^{F_S(U)}\right).         \tag{7a}
\]

Indeed, for `U={i:x_i=-y_i}`, split the `y`-energy into `I` off the
cut and `F,G` on its two directed halves.  Then
`Q_A(x)+Q_A(y)=2I`, `x^TRy=2(G-F)`, and
`|I|+|G-F|=max(|Q_(A^F)(y)|,|Q_(A^G)(y)|)`; the two half-cut flips are
switching-equivalent.  Equivalently, if `D_U` is `-1` on `U` and `+1`
off it,

\[
 A^{F_S(U)}
 ={A+D_UAD_U+D_UR-RD_U\over2}
 =\operatorname{Re}\left(Z_U^*(A+iR)Z_U\right),
 \quad Z_U={I+D_U\over2}-i{I-D_U\over2}.          \tag{7b}
\]

Thus the calculator is optimizing a precise shape: one orientation must keep
the norm of **every** outward half-cut neighbor below the `sqrt(2)` line.
Global minimality of `A` gives only the reverse lower bound.  The complete
proof, ordered-prefix total-variation specialization, and regression checks
are in the [orientation audit](NOTE_2026-09-02_ORIENTATION_STRUCTURE_AUDIT.md)
and [ordered-skew audit](NOTE_2026-09-02_ORDERED_SKEW_PREFIX_HALF_CUT.md).

## 3. What the calculator optimizes

For fixed `A`, every upper-triangular entry of `R` is a binary variable.
For each Boolean pair, `x^TRy` is a linear form in those variables.  Introducing
one integer variable `B` and imposing

\[
 B\ge |Q_A(x)+Q_A(y)|\mathbin{\pm}x^TRy
\]

for every pair gives a complete finite integer minimax model.  Global signs
allow `x_1=y_1=1`, and `R` and `-R` allow one orientation bit to be fixed;
neither symmetry reduction changes the optimum.

The replay uses literal enumeration for `n=5,6` and requires CP-SAT status
`OPTIMAL` for the complete models at `n=7,8`.  The symmetric matrices have
`Phi(A)=m_n` using the independently recorded exact `m_n` table.

| `n` | `M=m_n` | `min_R B` | `min_S max_U Phi(A^F)` | `min_R max(W-D)` | target `2(sqrt(2)-1)M` | zero-error target |
|---:|---:|---:|---:|---:|---:|:---:|
| 5 | 4 | 16 | 8 | 8 | 3.3137 | fails |
| 6 | 5 | 18 | 9 | 8 | 4.1421 | fails |
| 7 | 9 | 22 | 11 | 4 | 7.4558 | passes |
| 8 | 10 | 28 | 14 | 8 | 8.2843 | passes |

The exact upper envelopes for the stored optimal orientations are:

```text
n=5: (D,Wmax) = (0,8),(4,12),(8,12)
n=6: (D,Wmax) = (0,8),(2,10),(4,12),(8,16),(10,14)
n=7: (D,Wmax) = (0,4),(4,8),(8,12),(12,16),(16,20)
n=8: (D,Wmax) = (0,4),(2,10),(4,12),(6,14),(8,16),(10,18),
                  (12,16),(14,18),(16,16),(18,18),(20,16)
```

At `n=7` the active boundary is literally `W=D+4`.  At `n=8` it is
`W=D+8` on `D=2,4,6,8,10`, missing the sharp real-valued allowance by only
`0.28427...`.  The visible slope one is forced by (4); the substantive
quantity is its optimized intercept.

## 4. An exact calculus invariant from the same data

For the discrete derivative

\[
 d_iQ_A(x)={Q_A(x)-Q_A(x^{\oplus i})\over2}
           =x_i\sum_{j\ne i}A_{ij}x_j,
\]

the random variable on the right is a sum of `n-1` independent signs,
regardless of `A`.  If

\[
 \mu_k=\mathbb E|\varepsilon_1+\cdots+\varepsilon_k|,
 \qquad
 K_n=\max_A\operatorname{Inf}_1(Q_A/\Phi(A)),
\]

then every signing has total `L^1` influence
`Inf_1(Q_A)=n mu_(n-1)`, and therefore

\[
 \boxed{\quad m_n={n\mu_{n-1}\over K_n},\qquad
 {m_n\over n^{3/2}}={\mu_{n-1}\over\sqrt n}{1\over K_n}.\quad}       \tag{8}
\]

Since `mu_(n-1)/sqrt(n) -> sqrt(2/pi)`, the original MO sequence converges
if and only if the sharp equimodular quadratic influence constants `K_n`
converge.  The plot includes their exact values for `n=2,...,15`.  Order
eleven uses an independently replayed exact certificate; order twelve then
follows from its checked witness, monotonicity, and parity.  Orders thirteen
through fifteen use published certificate packages with explicitly separated
local checks and external stream-completeness trust boundaries; see the
[order-eleven audit](NOTE_2026-09-02_EXTERNAL_M11_CERTIFICATE_AUDIT.md) and
[orders twelve--fifteen audit](NOTE_2026-09-02_EXTERNAL_N12_N15_CERTIFICATE_AUDIT.md).
They are
not monotone at these orders, so (8) is an exact analytic target rather than
a visual extrapolation.  The exact one-vertex increment identity, cut-code
form, and failure of naive tensorization are proved in the
[cut-code/spherical limit audit](NOTE_2026-09-02_CUT_CODE_SPHERICAL_LIMIT_AUDIT.md).

## 5. Honest conclusion and next analytic target

The finite transition between orders six and seven is a pattern, not an
asymptotic theorem.  The durable advance is the exact reformulation (3)--(5):
the open multiplier-two step is an `A`-weighted skew discrepancy problem on
complementary signed supports.  Any scalable proof must exploit the payment
`|X|+H=M-|I|`; controlling `|C|` independently is known to be too strong.

A useful probabilistic sufficient condition follows immediately.  For fixed
complementary signed supports `u,v` of sizes `a,b`, a uniformly random skew
signing makes `C=u^TRv` a sum of `ab` independent Rademacher signs.  Writing
`Tail_h(t)=P(|S_h|>t)` for a sum of `h` independent signs, the exact
orbit-counted first-moment criterion for existence of a suitable `R` is

\[
 {1\over8}\sum_{(u,v):ab>0}
 \operatorname{Tail}_{ab}\!\left((\sqrt2-1)M+M-|I(u,v)|\right)<1.         \tag{9}
\]

The factor `1/8` quotients the independent sign changes of `u,v` and their
swap.  No independence between bad events is asserted or needed.  The
ordinary Chernoff bound gives the more explicit sufficient condition

\[
 {1\over8}\sum_{(u,v):ab>0}
 2\exp\!\left(-{((\sqrt2-1)M+M-|I(u,v)|)^2\over2ab}\right)<1.            \tag{10}
\]

Equation (9) is only a sufficient criterion; no claim is made here that an
optimal `A` satisfies it.  Its value is that it converts the next proof task
into a weighted near-maximizer count, with larger energy slack automatically
buying a weaker discrepancy threshold.

Moreover, the independent-random-skew route cannot close at the optimal
scale at all.  Even granting the fictitious best slack `M` to every
central cut leaves exponential first-moment rate
`log(4)-4 alpha^2`, where `M=alpha n^(3/2)`; the known upper scale
`alpha<=1/2+o(1)` leaves `log(4)-1>0`.  Thus (9) cannot certify an
optimal signing, even if every internal energy vanished.  One must go beyond
an expectation/union-bound argument or replace independent random signs by a
structured discrepancy construction.  Energy-distribution refinement cannot
overcome the ideal threshold `sqrt(log(2)/2)=0.588705...>1/2`.  The exact
orbit count, threshold, and asymptotic obstruction are proved in the
[random-skew criterion audit](NOTE_2026-09-02_TWO_HALF_RANDOM_SKEW_CRITERION.md).

The feasibility problem itself is

\[
 \min_{r\in\{\pm1\}^{\binom n2}}\max_{x,y}
 \bigl(|c_{xy}\cdot r|-D_{xy}\bigr),
 \quad c_{xy,ij}=x_i y_j-x_jy_i.                  \tag{11}
\]

Each `c_xy/2` is a signed complete-bipartite cut.  Relaxing the signs to
`[-1,1]` is identically useless (`r=0` gives the endpoint-forced optimum
zero), so the required phenomenon is genuinely integral: partial coloring,
higher-order constraints, or an `A`-dependent algebraic orientation.

## 6. Reverse-engineering the limiting shape

Proposition 6.5i turns the plotted central saddle into a necessary
all-orders statement.  The dual-Gaussian rounding used for the `1/pi` lower
bound has an exact saturation-gap identity and a Hermite-rank-two overlap
bound.  Hence, if optimal signings approach that floor, its positive and
negative outputs have energies `+M-o(n^(3/2))` and `-M+o(n^(3/2))` while
their Hamming distance is `n/2+o(n)`.

If an orientation simultaneously attains the outgoing-half target, the
same argument on `K_0=[[A,R],[-R,-A]]` produces two balanced complementary
halves `u,v` with

\[
 Q_{K_0}(u)+Q_{K_0}(v)=o(n^{3/2}),\qquad
 u^TK_0v=\Phi(K_0)-o(n^{3/2}),
\]

and each axial energy has magnitude at most
`Phi(K_0)/2+o(n^(3/2))`.  Thus the middle saddle drawn by the calculator is
not a finite-order accident: it is forced at a sharp lower-floor solution.
This does not construct the orientation.  The exact proof and variance
bound are in the
[Gaussian saturation audit](NOTE_2026-09-02_GAUSSIAN_SATURATION_CENTRAL_SADDLE.md).

## 7. Replay

```bash
python scripts/original_mo_two_half_geometry.py --solve
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. pytest -q \
  tests/test_original_mo_two_half_geometry.py \
  tests/test_orientation_structure_audit.py
```

Outputs:

- `evidence/original_mo_two_half_geometry.json`
- `evidence/original_mo_two_half_geometry.png`

The JSON records every matrix, orientation, envelope, maximizing witness,
solver status, replay statistic, the exact `K_n` sequence through `n=15`,
and the SHA-256 of the source exact-value table.  The plot shows both the
original small `alpha_n` values and the exact diagonal envelopes.  These are
exhaustive finite certificates only; the original MathOverflow question
remains open.
