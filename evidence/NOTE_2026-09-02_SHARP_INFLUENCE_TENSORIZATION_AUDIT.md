# Sharp equimodular influence: dimension and tensorization audit

**Status:** proved exact reformulation, exact one-step bounds, a probabilistic
padding lemma, and method obstructions.  This note does **not** prove that
`K_n` converges.  Its purpose is to separate the genuinely useful influence
identity from four dimension-change arguments that do not preserve the flat
coefficient class at the required constant.

## 1. The identity and the exact size of a sufficient extension theorem

Put

\[
 a_n=n\mu_{n-1},\qquad
 \mu_k={\bf E}|\varepsilon_1+\cdots+\varepsilon_k|,
 \qquad K_n={a_n\over m_n}.
\]

Proposition 5.3 gives this identity because every complete equimodular
quadratic has total `L1` influence exactly `a_n`.  In particular,

\[
 K_n\alpha_n={\mu_{n-1}\over\sqrt n}
 \longrightarrow\sqrt{2/\pi}.
 \tag{1}
\]

Thus convergence of `K_n` is *equivalent* to the original MO problem, not a
weaker consequence of it.

The elementary Rademacher recurrences are

\[
 \mu_{2r}=\mu_{2r-1},\qquad
 \mu_{2r+1}={2r+1\over2r}\mu_{2r}.
\]

Consequently, with `c_n=a_(n+1)/a_n`, one has exactly

\[
 c_n=
 \begin{cases}
 (n+1)/n,&n\text{ even},\\[2mm]
 (n+1)/(n-1),&n\text{ odd}.
 \end{cases}                                      \tag{2}
\]

Monotonicity `m_n<=m_(n+1)` and the elementary one-vertex extension
`m_(n+1)<=m_n+n` therefore give

\[
 c_nK_n{m_n\over m_n+n}\le K_{n+1}\le c_nK_n.     \tag{3}
\]

Using `m_n>=n sqrt(n-1)/pi`, the possible downward loss in (3) is still
`O(n^(-1/2))`, while the possible upward gain is `O(n^(-1))`.  Neither is
summable.  Bounds of this form permit bounded nonconvergent model sequences
such as a sufficiently small perturbation of `2+sin(log n)`.

Equation (2) also states precisely what a successful one-step theorem would
have to say.  It would suffice to prove, for some nonnegative `epsilon_n`,

\[
 m_{n+1}\le c_nm_n(1+\epsilon_n),
 \qquad \sum_n\epsilon_n<\infty.                  \tag{4}
\]

Indeed (4) gives `K_(n+1)>=K_n/(1+epsilon_n)`.  Since (1) and the known
two-sided bounds keep `K_n` bounded above and away from zero, `log K_n` would
have summable total downward variation and hence would converge.  In
additive form, an error `e_n` in (4) must satisfy
`sum e_n/m_n<infinity`; for example `e_n=O(n^(1/2-delta))`, or a Dini
logarithmic improvement on `sqrt(n)`, would suffice.  Merely proving an
`O(n)` extension is off by a factor `sqrt(n)`.

## 2. Principal restriction gives only the known one-sided inequality

Let `A` have order `N`, and let `S` be any `n`-set.  For a fixed Boolean
assignment on `S`, extend the other coordinates by independent uniform
signs.  Conditional expectation kills every edge not internal to `S`, so

\[
 Q_{A[S]}(x_S)={\bf E}[Q_A(X)\mid X_S=x_S].
\]

Jensen's inequality gives `Phi(A[S])<=Phi(A)` for every `S`.  On an optimal
order-`N` form this recovers `m_n<=m_N`, and at the influence level only

\[
 K_N\le {a_N\over a_n}K_n.                         \tag{5}
\]

For `N/n` bounded away from one, the factor in (5) is asymptotic to
`(N/n)^(3/2)`, so (5) has no almost-monotonic content.  Averaging the
restrictions does not reverse the problem: for a fixed global `x`, the
average restricted energy is a scalar multiple of `Q_A(x)`, but
`E max >= max E` is the wrong direction for finding a proportionally small
restricted supremum.

Fixing variables is worse for the proposed nesting argument.  It creates
linear terms.  Those terms belong to the unrestricted degree-at-most-two
class, not to the complete homogeneous equimodular class defining `K_n`.

## 3. The monotone unrestricted influence constant is a different object

Define

\[
 \mathcal B_{2,n}=\sup_{0\ne f\text{ homogeneous quadratic on }\{\pm1\}^n}
 {\operatorname{Inf}^{(1)}(f)\over\|f\|_\infty},                 \tag{6}
\]

with arbitrary real coefficients.  Zero-padding makes
`mathcal B_(2,n)` nondecreasing.  The dimension-free bounded-degree
influence theorem of
[Filmus--Hatami--Keller--Lifshitz](https://arxiv.org/abs/1404.3396)
makes it bounded, so this unrestricted sequence does converge.

It is not `K_n`.  Already on four variables,

\[
 f(x,y,z,w)={x(z+w)+y(z-w)\over2}
\]

is Boolean-valued and homogeneous quadratic, hence has norm one and total
`L1` influence two.  Thus `mathcal B_(2,4)>=2`.  On the other hand the
complete equimodular value is

\[
 K_4={4\mu_3\over m_4}={6\over4}={3\over2}.        \tag{7}
\]

The zeros in `f` are essential: zero-padding, restriction, and convex
symmetrization all leave the flat class.  Therefore convergence of (6)
cannot be transferred to `K_n` without a new **asymptotic flattening
theorem** proving that arbitrary near-extremizers can be completed to dense
equal-modulus coefficient arrays with `1+o(1)` loss.  None of the standard
padding or blow-up operations below has that property.

Permutation/switching averages have the same defect.  Convexity can reduce
the sup norm of the averaged polynomial, but its coefficients acquire
different magnitudes (often zero), and its `L1` influence is no longer the
universal numerator `a_n`.  Selecting one member of the orbit merely returns
the original norm.  Independent rounding back to signs introduces a dense
quadratic discrepancy term, for which the available uniform control is on
the leading `n^(3/2)` scale rather than `o(n^(3/2))`.

## 4. Exact random-padding lemma, and why it is too large

Let `N=n+k` and

\[
 E={N\choose2}-{n\choose2}={k(2n+k-1)\over2}.
\]

Start with an optimal order-`n` signing and fill the `E` new edges with
independent signs.  For each fixed Boolean state the new-edge contribution
is a sum of `E` independent Rademachers.  Hoeffding followed by a union bound
over the `2^(N-1)` projective states proves the following deterministic
existence bound:

\[
 \boxed{\quad
 m_N\le m_n+\sqrt{\,2E(N\log2+1)\,}.
 \quad}                                             \tag{8}
\]

For `k=o(n)`, the added term is `O(n sqrt(k))`, i.e. relative size
`O(sqrt(k/n))` on the `n^(3/2)` scale.  This gives only local multiplicative
continuity.  At `k=1` it is order `n`, whereas (4) needs a correctly centred
order-`sqrt(n)` increment.  For `k=theta n`, (8) has a full
`Theta(n^(3/2))` error and cannot compare distinct limit points.  Padding is
therefore quantitatively weaker than the already isolated multiplier-two
and multiplier-three diamonds.

## 5. Ordinary graph blow-up has the wrong exponent

There is also an exact obstruction to the usual lexicographic graph
blow-up.  Replace each vertex of an order-`n` signing `A` by `r` vertices,
put the constant block `a_ij J_r` between distinct vertex blocks, and put an
optimal order-`r` signing `B_i` inside each block.  Call the resulting
complete signing `C`.

If `s_i` is the signed sum in block `i`, the cross-block energy is
`Q_A(s_1,...,s_n)`.  Since the absolute value of a multi-affine polynomial
on a box is maximized at a vertex,

\[
 |Q_A(s)|\le r^2\Phi(A).
\]

The internal blocks contribute at most `n m_r`.  Conversely, take every
block constant according to a maximizing state of `A`; the cross term is
`+-r^2 Phi(A)` and the internal contribution has absolute value at most
`n m_r`.  Hence the blow-up obeys the two-sided estimate

\[
 r^2\Phi(A)-nm_r\le\Phi(C)\le r^2\Phi(A)+nm_r.     \tag{9}
\]

For fixed `n` and `r->infinity`, `m_r=O(r^(3/2))`, so (9) is asymptotic to
`r^2 Phi(A)`.  A ratio-preserving influence lift would instead need order
`r^(3/2) Phi(A)` (with the fixed factor dictated by `a_(nr)/a_n`).  The
standard graph blow-up therefore destroys the ratio by a factor of order
`sqrt(r)`; vertex switches do not change this calculation.

## 6. Balanced tensors return to the live diamond

To avoid the `r^2` constant-block obstruction, the inter-block sign matrices
must be balanced.  The resulting quantity is vector-valued:

\[
 \max_{x_i\in\{\pm1\}^r}
 \left|\sum_{i<j}a_{ij}x_i^T H_{ij}x_j
       +\text{internal terms}\right|.              \tag{10}
\]

The scalar bound `Phi(A)` does not control (10) with constant one by an
ordinary Hilbert-space or Grothendieck relaxation.  The elementary CHSH
matrix

\[
 H=\begin{pmatrix}1&1\\1&-1\end{pmatrix}
\]

has scalar bilinear norm `2`, while unit vectors
`u_1 perpendicular u_2`,
`v_1=(u_1+u_2)/sqrt(2)`, and
`v_2=(u_1-u_2)/sqrt(2)` give vector value `2sqrt(2)`.  Thus a generic
scalar-to-vector tensor step already loses a fixed factor `sqrt(2)`, not
`1+o(1)`.  This does not rule out an `A`-specific integral construction, but
it rules out the unqualified vector-relaxation shortcut.

At block size two the sharp balanced construction is exactly the
equal-endpoint symmetric/skew Hadamard frame of Propositions 6.4--6.6.  Its
required `r^(3/2)` estimate is precisely the mixed-state diamond, with the
Hamming-central/joint-energy residue (6.20) still open.  At multipliers two
and three, an almost ratio-preserving tensor inequality for `K_n` is
therefore the existing two-ray convergence gate in new notation, not a way
around it.

## 7. Audit verdict

The influence identity is exact and conceptually useful: it identifies the
MO question with asymptotic regularity of the sharp flat degree-two
influence constant.  The four natural dimension arguments have precise
failure modes:

1. principal restriction gives only (5), with a factor of order
   `(N/n)^(3/2)`;
2. unrestricted symmetrization/zero-padding converges, but for the strictly
   larger constant (6), and flattening is missing;
3. random padding has an order-`n` one-vertex error instead of the required
   Dini-improved order-`sqrt(n)` error;
4. ordinary graph blow-up scales as `r^2`, while balanced block tensors are
   exactly the unresolved vector-valued diamond.

A genuinely new theorem would have to supply either the sharp extension
(4), an asymptotic flattening theorem with `1+o(1)` loss, or the already
identified ratio-preserving estimates at two multiplicatively independent
rays.  No claim of convergence follows from the present audit.
