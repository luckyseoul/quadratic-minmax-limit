# Exact minimum of the two-block diamond at Paley order 6

2026-11-02. Finite exact result, not a convergence proof.
The original MathOverflow limit remains OPEN. Nothing in this note
is entered in `src/original_mo_status.py` as a completion claim.

## Statement

Let $C$ be the Paley conference matrix of order $6$ ($q=5$),
as in `src/minmax_quadratic.py: paley_conference_matrix(5)`.
Put $A = B = C$. Then $\Phi(C) = 5 = m_6$, so $A$ and $B$ are
exact order-$6$ minimizers.

For a sign matrix $D \in \{\pm 1\}^{6 \times 6}$ and $\tau = -1$,
define

$$
\Delta(A,B;D)
  = \max_{x,y \in \{\pm1\}^6}
    \bigl(|Q_A(x) - Q_A(y)| + |x^T D y|\bigr),
\qquad
Q_A(x) = \tfrac12 x^T A x .
$$

Then

$$
\boxed{\min_{D} \Delta(A,B;D) = 18.}
$$

The minimum is over **all** sign matrices $D$, not only over the
switching orbit of a fixed seed. Hence it is also the minimum over
every gauge orbit: no gauge, and no free seed, beats $18$ for this
block pair.

## Why this matters, and what it does not do

The Section 4 fiber statement for the multiplier-two ray, at zero
error, would require $\Delta \le 2^{3/2} m_6 = 2\sqrt2 \cdot 5
\approx 14.14$. The exact minimum $18$ misses that target by
$3.86$, about $27\%$. So the zero-error two-block diamond is false
for this optimal block pair.

This does **not** refute the Dini-error version of the fiber
statement at $n=6$ alone. With $E_0 = 0$ in Section 4.1, formula
(13) gives $\Omega_2(6) = 2^{3/2}\log 2 / \log(7)^2 \approx 0.517$,
i.e. a raw budget of about $7.6$ at order $12$, taking the allowed
total to about $21.7 > 18$. The asymptotic question is whether the
excess ratio $(18 - 2\sqrt2 \cdot 5) / 12^{3/2} \approx 0.093$
shrinks with $n$. This note answers no asymptotic question, and
it says nothing about $n = 10$, about other optimal blocks, or
about the $\tau = +1$ case (whose $Q$-value distribution is the
same here, but which was not the case run).

## Upper bound: an explicit $D$ with $\Delta = 18$

The following sign matrix attains $\Delta = 18$, verified by direct
enumeration over all $64 \times 64$ pairs $(x,y)$:

```text
D = [[ 1 -1  1  1 -1  1]
     [ 1  1 -1 -1 -1 -1]
     [-1  1  1 -1  1  1]
     [ 1  1 -1  1  1 -1]
     [ 1 -1  1 -1  1 -1]
     [-1  1  1  1  1 -1]]
```

Verification (with $C$ from the repo):

```python
import sys, numpy as np
sys.path.insert(0, 'src')
from minmax_quadratic import paley_conference_matrix
C = paley_conference_matrix(5)
n = 6
X = np.array([[1 if (m >> i) & 1 else -1 for i in range(n)]
              for m in range(1 << n)])
Q = np.array([int(0.5 * x @ C @ x) for x in X])
D = np.array([[ 1,-1, 1, 1,-1, 1],
              [ 1, 1,-1,-1,-1,-1],
              [-1, 1, 1,-1, 1, 1],
              [ 1, 1,-1, 1, 1,-1],
              [ 1,-1, 1,-1, 1,-1],
              [-1, 1, 1, 1, 1,-1]])
B = X @ D @ X.T
assert set(np.unique(D).tolist()) == {-1, 1}
assert np.max(np.abs(Q[:, None] - Q[None, :]) + np.abs(B)) == 18
```

The $Q$-value distribution for $C$ is $\{-5: 12,\ -3: 20,\ 3: 20,\
5: 12\}$ over the $64$ states, so $|Q(x)-Q(y)| \le 10$ throughout.

## Lower bound: CP-SAT infeasibility for $\Delta \le 17$ and $\le 16$

The lower bound is a finite CP-SAT proof, not a hand argument.
Formulation, for a fixed target $\Delta_{\max}$:

- Variables $d_{ij} \in \{0,1\}$, $D_{ij} = 2d_{ij} - 1$.
- $x, y$ are restricted to representatives with $x_0 = y_0 = +1$.
  This is valid: $Q_A(-x) = Q_A(x)$ and $|(-x)^T D y| = |x^T D y|$,
  and similarly for $y$, so every pair has a representative in
  this set with the same contribution.
- $d_{00} = 1$ is fixed. This is valid: $D \mapsto -D$ leaves
  $|x^T D y|$ unchanged, so every $D$ has a global-sign twin
  with $d_{00} = 1$ and the same $\Delta$.
- For every representative pair $(x,y)$, with
  $\delta = |Q_A(x) - Q_A(y)|$ and bound $b = \Delta_{\max} - \delta$,
  the two linear constraints $|x^T D y| \le b$ are imposed, with
  $x^T D y = \sum_{ij} (2d_{ij}-1) x_i y_j$ expanded linearly.

No other symmetry breaking was used. In particular, fixing a
whole row or column of $D$ to $+1$ was **not** used: that is not
a symmetry of $\Delta$ with $A$ fixed, because flipping a row of
$D$ also changes which $x$ attains a given $Q_A(x)$. An earlier
run that used that break was discarded for exactly this reason.

Results, OR-Tools CP-SAT, 4 workers:

- $\Delta \le 15$: INFEASIBLE (full $64 \times 64$ formulation,
  no representative reduction).
- $\Delta \le 17$: INFEASIBLE (representative formulation above).
- $\Delta \le 16$: INFEASIBLE (representative formulation above).
- $\Delta \le 18$: feasible, by the explicit $D$ above.

Together these give $\min_D \Delta = 18$ exactly.

Two cautions on the solver record, kept here so the proof is not
overstated:

- A separate CP-SAT *optimization* run (minimize $\Delta$) timed
  out at 300 seconds with status FEASIBLE and value $18$, not
  OPTIMAL. The exact minimum does **not** rest on that run; it
  rests on the fixed-target infeasibility runs for $17$ and $16$
  plus the feasible $D$ at $18$.
- The fixed-target runs for $17$ and $16$ in the *unreduced*
  formulation timed out (UNKNOWN) at 120 seconds each. The
  infeasibility proof uses the representative formulation, which
  is valid for the reasons above and which returned INFEASIBLE
  for both targets.

## Prior-art check

This exact minimum over all $D$ at $n=6$ was not found in the
repo's notes. The closest prior material is the two-half random
skew criterion and its first-moment barrier
(`NOTE_2026-09-02_TWO_HALF_RANDOM_SKEW_CRITERION.md`), which
rules out an *independent random* skew by expectation and
explicitly says no routine variable local-lemma repair exists.
That is a statement about random $R$, not about the exact minimum
over all $D$ computed here. The Sept 28 hybrid / balanced-core
notes are a different (one-vertex) route and were not reopened.

## Open next step

The same exact question at $n=10$ (Paley $p=3$, $\Phi(C)=15$,
$m_{10}=13$, zero-error target $2\sqrt2 \cdot 13 \approx 36.77$)
is the next test. Its CP-SAT formulation has $100$ variables and,
after the representative reduction, $512 \times 512$ pairs, so a
direct fixed-target proof there is far larger than the $n=6$ case
and was not attempted in this note.
