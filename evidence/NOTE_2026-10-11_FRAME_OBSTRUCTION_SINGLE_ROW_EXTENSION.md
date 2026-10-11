# Frame obstruction to single-row extension, block bound, and Fekete counterexamples

**Status:** proved lemmas and exact finite computations. Nothing here proves or
refutes convergence of `alpha_n`; the original limit stays OPEN, and nothing was
entered in `src/original_mo_status.py`. Every statement carries one of three
labels.

- **PROVED**: complete proof in this note.
- **COMPUTED**: exact finite computation; the script path and its output are
  given in Section 10.
- **HEURISTIC**: not proved.

## 0. Conventions and constants

`Q_A(x) = sum_{i<j} A_ij x_i x_j`, `Phi(A) = max_x |Q_A(x)|` over
`x in {+1,-1}^n`, `m_n = min_A Phi(A)`, `alpha_n = m_n / n^{3/2}`,
`r_n = m_n ((1+1/n)^{3/2} - 1)`. A class is a pair `{x,-x}`. It is *tight* for
`A` if `|Q_A(x)| = Phi(A)`. `S_A` is the set of tight classes, `N0 = |S_A|`, and
`D(x) = Phi(A) - |Q_A(x)| >= 0` is the deficit of any state. For a row
`a in {+1,-1}^n`, `A+a` is the order `n+1` matrix with last row and column `a`.
`ext(A) = min_a Phi(A+a) - Phi(A)`. Exact values `m_2,...,m_15 = 1,3,4,4,5,9,10,
12,13,17,18,20,21,27` are those of `evidence/exact_m_table.json` and
`solution.md` (Section 6); values for `n <= 9` are re-derived in Section 10.

Two lower-bound constants occur in the repository. The author-reviewed constant
is `0.3258` (`README.md`, `CORE.md` Section 4, `B_tilt = 0.3258530...`); the
constant `0.3365` (`gamma_* = 0.336493...`, field-plus-spin) appears in
`evidence/NOTE_2026-10-01_DOUBLING_WINDOW_DEBT.md`. Both are displayed wherever
a lower bound for `alpha_n` enters. Nothing below depends on either.

## 1. The one-row frame obstruction

**Row formula (PROVED).** For every `A` of order `n` and every row `a`,

    Phi(A+a) = max_x ( |Q_A(x)| + |a.x| ),

since the sign of the new coordinate may be chosen to make `Q_A(x)` and the new
terms agree in sign. In particular, for `x` tight, `Phi(A+a) >= Phi(A) + |a.x|`.
This is the elementary half of the stable-state reduction of
`evidence/NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md`.

**Proposition 1.1 (Parseval bound; PROVED).** Let `S` be a finite set of vectors
in `{+1,-1}^n`, `M = sum_{x in S} x x^T`, and `lambda = lambda_min(M)`. For every
`a in {+1,-1}^n`,

    max_{x in S} |a.x| >= sqrt( n lambda / |S| ).

*Proof.* `sum_{x in S} (a.x)^2 = a^T M a >= lambda |a|^2 = n lambda`, and the
left side is at most `|S| max_x (a.x)^2`. QED

**Corollary 1.2 (frame obstruction; PROVED).** If `sum_{x in S_A} x x^T = N0 I_n`
(the tight classes form a tight frame), then every row `a` satisfies
`max_{x in S_A} |a.x| >= sqrt(n)`, hence

    Phi(A+a) >= Phi(A) + sqrt(n)   for every a, so ext(A) >= sqrt(n).

If moreover `Phi(A) = m_n`, then, since `r_n = (3/2) alpha_n sqrt(n) (1+O(1/n))`,
every one-row extension of `A` exceeds the neutral increment `r_n` by at least
`(1 - 1.5 alpha_n) sqrt(n) - O(n^{-1/2})`. The coefficient is

| `alpha_n` | 0.3258 | 0.3365 | 0.40 | 0.45 | 0.50 |
|---|---|---|---|---|---|
| `1 - 1.5 alpha_n` | 0.511 | 0.495 | 0.400 | 0.325 | 0.250 |

and it is at least `0.25` whenever `alpha_n <= 1/2`. The bound is a pointwise
inequality in `a`, so it holds for every distribution of rows. It concerns
extension of that one matrix `A`; it says nothing about `m_{n+1}` realised by a
different matrix.

**Relation to known frame statements.** `solution.md` Proposition 15.27 is a tight
frame identity for the Boolean eigenvector maximisers `Max_+` of a conference
matrix, a frame in the eigenspace `V_+` only. Erden (doi
10.5281/zenodo.22018437, Section 7, "the argmax set does not always span")
observes that a zero-cost extension needs `a.x = 0` on every tight class, which
rank `n` of the tight classes rules out. Corollary 1.2 is a quantitative form of
the second statement for the full tight set (both orientations), with cost
`sqrt(n)` instead of `> 0`. Proposition 1.1 gives the bound for any
tight set of full rank.

**Computation 1.3 (COMPUTED).** For the Paley conference matrices `C_{q+1}`
(`q = 5, 13, 17, 29`), `scripts/one_vertex_tight_frame_paley.c` enumerates all
`2^n` classes (`n = q+1`) exactly:

| `n` | `Phi` | `N0` | Gram of `S` = `N0 I` | min over rows of `max_{S}|a.x|` | `sqrt(n)` |
|---|---|---|---|---|---|
| 6 | 5 | 12 | yes | 4 | 2.449 |
| 14 | 21 | 156 | yes | 6 | 3.742 |
| 18 | 33 | 204 | yes | 6 | 4.243 |
| 30 | 75 | 812 | yes | 8 | 5.477 |

The Gram identity holds as an exact integer identity. The values `N0 = 12, 156,
204` agree with `solution.md` (Paley maximiser counts); `N0 = 812` at `n = 30`
was not recorded before. `Phi = 5, 21` equal `m_6, m_14`; the order 6 and 14
minimisers are unique up to switching (Erden, Table 1). Consequently every
order-14 matrix with `Phi = 21` has `ext = 6` (also computed over all rows in
Section 5), so `m_15 - m_14 = 6 > r_14 = 2.290` (with `m_15 = 27`). At `n = 30` every row costs at
least `8 > r_30 = 3.78` (computed with `m_30` replaced by `Phi(C_30) = 75`, an
upper bound for `m_30`).

The frame property is verified for these four orders only. A group-theoretic
explanation (the tight set is a union of orbits of the automorphism group of the
two-graph, which acts on both eigenspaces of `C`) is HEURISTIC; no proof for
general `q` is given, and no frame is claimed for any other matrix.

**Limitation to a single added row (PROVED).** The bound `sqrt(n)` is a bound for
one added row. For `d` added vertices the correct lower bound is of order
`d beta sqrt(n)` with `beta < 1` (Section 2); Corollary 1.2 does not extend to
blocks.

## 2. Block extensions

Let `B` have order `n`, `Phi(B) = F`, `D(x) = F - |Q_B(x)|`. For
`C in {+1,-1}^{n x d}` and a symmetric zero-diagonal `+-1` matrix `E` of order
`d`, put `ext = [[B, C],[C^T, E]]` and `g(C) = max_x ( ||C^T x||_1 - D(x) )`.

**Proposition 2.1 (block formula; PROVED).**
`g(C) - Phi(E) <= Phi(ext) - F <= g(C) + Phi(E)`, and `Phi(E) <= d(d-1)/2`.

*Proof.* `Phi(ext) = max_{x,y} |Q_B(x) + x^T C y + Q_E(y)|`. Upper bound: the
triangle inequality gives `|Q_B(x)| + ||C^T x||_1 + Phi(E)`, and
`|Q_B(x)| + ||C^T x||_1 = F + (||C^T x||_1 - D(x))`. Lower bound: for a given `x`
let `s = sign Q_B(x)` (`+1` if zero) and `y_j = sign(s x.C_j)`; then
`s (Q_B(x) + x^T C y + Q_E(y)) >= |Q_B(x)| + ||C^T x||_1 - Phi(E)`. QED

For `d >= 2` no reduction to stable states is available (flipping one
coordinate of `x` changes `||C^T x||_1` by up to `2d` but raises `s Q_B` by only
`2` in general), so all states enter `g(C)`.

**Proposition 2.2 (averaged bound; PROVED).** Let `S` be the tight set of `B`
(`D = 0` on `S`) and

    beta(S) = min_{c in {+-1}^n} avg_{x in S} |x.c| / sqrt(n).

Then `g(C) >= d beta(S) sqrt(n)` for every `C`, so
`Phi(ext) - F >= d beta(S) sqrt(n) - d(d-1)/2`.

*Proof.* `g(C) >= avg_{x in S} ||C^T x||_1 = sum_j avg_{x in S} |x.C_j|
>= d beta sqrt(n)`. Apply Proposition 2.1. QED

**Computation 2.3 (COMPUTED).** Same script, exact minimum over all `2^{n-1}`
columns `c`:

| `n` | `F` | `N0` | `min_c sum_S |x.c|` | `min_c avg` | `beta` | `r_n / sqrt(n)` |
|---|---|---|---|---|---|---|
| 6 | 5 | 12 | 20 | 1.6667 | 0.6804 | 0.531 |
| 14 | 21 | 156 | 392 | 2.5128 (= 98/39) | 0.6716 | 0.612 |
| 18 | 33 | 204 | 576 | 2.8235 (= 48/17) | 0.6655 | 0.657 |
| 30 | 75 | 812 | 2940 | 3.6207 (= 105/29) | 0.6610 | 0.690 |

For `n = 6, 14` the value `F` is `m_n`. For `n = 18, 30`, `F = Phi(C_n)` is only an
upper bound for `m_n` (the repository bank value is `m_18 <= 33`), so
`r_n / sqrt(n)` in those rows is an upper estimate of the true allowance. At
`n = 14, 18` the per-vertex block bound `beta sqrt(n)` already exceeds `r_n`
(`beta` exceeds `r_n / sqrt(n)` by 0.060 and 0.008); at `n = 30` it falls below it
(by 0.029), whereas the one-row bound `sqrt(30) = 5.48` exceeds `r_30 = 3.78`.
The correction `d(d-1)/2` in Proposition 2.2 is small against `d beta sqrt(n)`
only for `d = o(sqrt(n))`; the table compares the per-vertex coefficient
`beta sqrt(n)` with `r_n`.

**Proposition 2.4 (fourth-moment form; PROVED).** If `sum_{x in S} x x^T = N0 I`
and `avg_{x in S} (x.c)^4 <= K n^2` for all `c in {+-1}^n`, then
`avg_{x in S} |x.c| >= sqrt(n/K)` for all `c`, so `beta(S) >= K^{-1/2}`.

*Proof.* With `t = x.c`, the average of `t^2` over `S` is `c^T (M/N0) c = n`.
Holder with exponents `3/2` and `3` applied to `t^2 = |t|^{2/3} |t|^{4/3}` gives
`E t^2 <= (E|t|)^{2/3} (E t^4)^{1/3}`, so
`E|t| >= (E t^2)^{3/2} / (E t^4)^{1/2} >= n^{3/2} / (sqrt(K) n)`. QED

For Gaussian kurtosis (`K = 3`) the per-vertex block cost would be at least
`0.577 sqrt(n)`, above the allowance `1.5 alpha sqrt(n)` exactly when
`alpha < 0.385` (so for `alpha = 0.3258` and `0.3365`, and not for
`alpha >= 0.40`). Whether tight sets satisfy a fourth-moment bound with `K = 3` is
HEURISTIC; the computed `beta >= 0.66` forces `K >= 1/beta^2 >= 2.2` for the
extremal column.

**Reachability (HEURISTIC).** The bound of Proposition 2.2 is an average over
tight classes; a block `C` with `d = o(n)` columns must control the maximum over
all `N0` classes, and a Gaussian-fluctuation estimate puts that maximum above the
mean by of order `sqrt(d n log N0)`. Reaching `d beta sqrt(n)` therefore plausibly
needs `d` large compared with `log N0`, and for `d` of order `n` the comparison
with a linear allowance `d r_n` is no longer the right one. No block of cost near
`d beta sqrt(n)` is exhibited here.

**Modifying rows is a block extension (PROVED).** Replacing the rows and columns
of a set `K` of `k` vertices of `A` and adding one vertex gives exactly a block
extension of `B = A - K` (order `n-k`) with `d = k+1`, and conversely every block
extension of `B` arises. Hence Propositions 2.1 and 2.2 apply with `F`
replaced by `Phi(B) <= Phi(A)`.

## 3. Chain criterion, and what summability of the excess requires

Write `e_{n,d} = m_{n+d} - m_n (1 + d/n)^{3/2}`. Then, exactly,

    alpha_{n+d} - alpha_n = e_{n,d} / (n+d)^{3/2}.                      (3.1)

With `d = 1`, `e_{n,1} = m_{n+1} - m_n - r_n`, which is at most the excess `rho_n`
of `evidence/SCOPE_FREEZE_2026-09-25_ONE_VERTEX_ONLY.md` for any admissible row.

**Proposition 3.2 (chain criterion; PROVED).** Suppose `theta(n) -> 0` and that
for every large `n_0` there is an unbounded chain
`n_0 = n^(0) < n^(1) < ...` with `n^(i+1) - n^(i) <= theta(n^(i)) n^(i)` and

    sum_i e_{n^(i), d_i} / (n^(i+1))^{3/2} <= eta(n_0),   d_i = n^(i+1) - n^(i),

where `eta(n_0) -> 0` as `n_0 -> infinity`. Then `alpha_n` converges.

*Proof.* Let `L = liminf alpha_n`, `U = limsup alpha_n`, `eps > 0`. Choose `n_0`
with `alpha_{n_0} <= L + eps`, `eta(n_0) <= eps`, and `theta(n) <= eps` for
`n >= n_0`. Telescoping (3.1) along the chain gives
`alpha_{n^(i)} <= alpha_{n_0} + eta(n_0) <= L + 2 eps` for every `i`. For
`N >= n_0` pick `i` with `n^(i) <= N < n^(i+1)`. Since `m` is nondecreasing
(`CORE.md` Section 2), `alpha_N <= alpha_{n^(i+1)} (n^(i+1)/N)^{3/2}
<= (L + 2 eps)(1 + eps)^{3/2}`. Hence `U <= (L + 2 eps)(1 + eps)^{3/2}`, and
`eps -> 0` gives `U <= L`. QED

The sum is signed (negative `e` helps) and blocks are allowed. For `d_i = 1` it is
the Dini-summable excess criterion of `CORE.md` Section 4 (summability of
`rho_n / (n+1)^{3/2}`).

**Corollary 3.3 (PROVED).** `rho_n = O(sqrt(n) / (log n)^{1+eps})` for some
`eps > 0` suffices, since `sum 1/(n (log n)^{1+eps}) < infinity`. For blocks: if
`e_{n^(i), d_i} <= sum_{n^(i) <= k < n^(i+1)} eps_k` with a per-vertex excess
`eps_k = O(sqrt(k) / (log k)^{1+eps})`, then, because
`(n^(i+1))^{-3/2} <= k^{-3/2}` for these `k`, the chain sum in Proposition 3.2 is at
most `sum_{k >= n_0} eps_k k^{-3/2}`, which tends to `0`.

**Proposition 3.4 (`rho_n = c sqrt(n)` is not enough; PROVED).** Let
`g(t) = a + b sin(omega t)` with `a = 0.42`, `b = 0.03`, `omega = 0.1`, and
`m_n = round(n^{3/2} g(ln n))`. Then `m` is nondecreasing, `m_{n+1} - m_n <= n-2`
for `n >= 9`, `alpha_n in [0.39, 0.45]` for `n >= 8` (inside the proved window
`0.3365 <= alpha_n <= 1/2`), and
`e_{n,1} = b omega sqrt(n) cos(omega ln n) + O(1) <= 0.003 sqrt(n) + O(1)`, yet
`alpha_n = g(ln n) + O(n^{-3/2})` oscillates forever.

*Proof.* `h(t) = n^{3/2} g(ln n)` has derivative
`n^{1/2}(1.5 g + g') >= n^{1/2}(1.5 * 0.39 - 0.003) > 0`, so `m` is nondecreasing
(rounding preserves monotonicity), and its increments are at most
`n^{1/2}(1.5 * 0.45 + 0.003) + 1 <= n - 2` for `n >= 9`. Rounding changes
`alpha_n` by `O(n^{-3/2})` and `e_{n,1}` by `O(1)`. Finally
`m_n (1+1/n)^{3/2} = (n+1)^{3/2} g(ln n) + O(1)`, so
`e_{n,1} = (n+1)^{3/2} [g(ln(n+1)) - g(ln n)] + O(1) = b omega sqrt(n)
cos(omega ln n) + O(1)`. QED

The model is not a statement about the true `m_n`; it shows that the scalar
inequalities monotonicity, `m_{n+1} - m_n <= n - 2`, the window, and
`e_{n,1} <= c sqrt(n)` do not imply convergence. Related records: the
2026-10-11 `HANDOFF.md` entry (an integer sequence inside the scalar axioms that
oscillates in `(0.42, 0.46]`) and the public claim register of Marcus Miccelli
(entry `superadditivity`). **COMPUTED sanity check:**
`scripts/log_periodic_alpha_model_check.py` (Section 10).

**Consequence for frame minimisers (PROVED from Corollary 1.2).** If `A` is a
minimiser of order `n` whose tight classes form a tight frame, then every
one-row extension has `rho_n >= (1 - 1.5 alpha_n) sqrt(n) - O(n^{-1/2})`, at least
`0.25 sqrt(n) - o(sqrt(n))` when `alpha_n <= 1/2 + o(1)`; this is larger than the
sufficient scale of Corollary 3.3 by a factor of at least `(log n)^{1+eps}`.
If frame minimisers existed for a set `N` of orders with `sum_{n in N} 1/n =
infinity`, the one-row inequality with summable `rho_n` would fail for `A`
chosen among those minimisers at orders in `N`, because
`sum_{n in N} rho_n / (n+1)^{3/2} >= c sum_{n in N} 1/n`. This does not exclude
other matrices at those orders, nor near-minimisers whose tight set is not a
frame.

## 4. First-order optimality: at least three near-tight classes

**Proposition 4.1 (PROVED).** If `Phi(A) = m_n > 2`, then at least three classes
satisfy `|Q_A(x)| >= Phi(A) - 2`.

*Proof.* Fix a pair `{i,j}` and let `A'` be `A` with `A_ij` negated. Minimality
gives `Phi(A') >= Phi(A)`, so there is `x` and a sign `s` with
`s Q_{A'}(x) >= Phi(A)`, where `Q_{A'}(x) = Q_A(x) - 2 A_ij x_i x_j`. Since
`s Q_A(x) <= Phi(A)`, this forces `s A_ij x_i x_j = -1` and
`s Q_A(x) >= Phi(A) - 2`. As `Phi(A) - 2 > 0`, `s = sign Q_A(x)`. Say that
`(x, s)` *covers* the pairs `{i,j}` with `s A_ij x_i x_j = -1`. Every pair is
covered by a class with `|Q_A| >= Phi(A) - 2`. A class with `s Q_A(x) = q` covers
exactly `(C(n,2) - q)/2 <= (C(n,2) - Phi(A) + 2)/2 < C(n,2)/2` pairs. Two classes
therefore cannot cover all `C(n,2)` pairs. QED

This is a lower bound of `3` only; it gives no upper bound `N0 <= eps n` and no
lower bound `N0 >= c n`. It is distinct from the `m_{n+1} <= m_n + n - 2` argument
of `CORE.md`, which uses a vector with `|Q| <= m_n - 2`.

## 5. One-row extension of every known minimiser, n = 5,...,14

**COMPUTED** by `scripts/one_vertex_minimiser_extension_table.py`, exact
enumeration of all rows, for one representative of every switching class of
minimisers. For `n <= 9` the lists come from the exhaustive enumeration of
`scripts/normalized_exhaustive_phi_counts.c`. For `n = 10,...,13` they were
obtained by annealing, and the class counts (2, 4, 1, 1) agree with Erden,
Table 1 (classes for `n = 4,...,14`: 1, 1, 1, 6, 2, 15, 2, 4, 1, 1, 1). Order 14
is the Paley matrix.

| `n` | `m_n` | classes | `N0` range | `ext` (multiplicity) | `r_n` | `ext - r_n` |
|---|---|---|---|---|---|---|
| 5 | 4 | 1 | 10 | 1 (1) | 1.258 | -0.26 |
| 6 | 5 | 1 | 12 | 4 (1) | 1.301 | +2.70 |
| 7 | 9 | 6 | 3-7 | 1 (4), 3 (2) | 1.996 | -1.00, +1.00 |
| 8 | 10 | 2 | 8 | 2 (2) | 1.932 | +0.07 |
| 9 | 12 | 15 | 12-28 | 1 (4), 3 (11) | 2.055 | -1.05, +0.95 |
| 10 | 13 | 2 | 40 | 4 (1), 6 (1) | 1.998 | +2.00, +4.00 |
| 11 | 17 | 4 | 5-24 | 1 (2), 3 (2) | 2.370 | -1.37, +0.63 |
| 12 | 18 | 1 | 20 | 6 (1) | 2.296 | +3.70 |
| 13 | 20 | 1 | 78 | 1 (1) | 2.352 | -1.35 |
| 14 | 21 | 1 | 156 | 6 (1) | 2.290 | +3.71 |

Observations (COMPUTED, from the table and the exact values of Section 0).

1. For every minimiser at `n = 6, 8, 10, 12, 14` every one-row extension exceeds
   `r_n`. At `n = 8` the margin is `+0.068`, so the excess
   `m_{n+1} - m_n - r_n` is positive at `n = 6, 8, 10, 14` (values `+2.70, +0.07,
   +2.00, +3.71`) and negative at `n = 5, 7, 9, 11, 12, 13` (values `-0.26,
   -1.00, -1.05, -1.37, -0.30, -1.35`). The positive excess is therefore not
   confined to `n = 2 mod 4`.
2. `m_{n+1} = m_n + min ext(A)`, the minimum over minimisers `A` of order `n`,
   for `n = 5,...,11, 13, 14` (for example `9 + 1 = 10`, `17 + 1 = 18`,
   `21 + 6 = 27`). It fails at `n = 12`: the unique order-12 minimiser has
   `ext = 6`, so `18 + 6 = 24 > 20 = m_13`, and every order-13 matrix with
   `Phi = 20` restricts to an order-12 matrix that is not a minimiser.
3. `ext = n - 2` at `n = 6` (value 4) is the equality case of the step bound
   `m_{n+1} <= m_n + n - 2` recorded in `CORE.md`.
4. Three or four points with positive excess do not distinguish `O(1)` from
   `c sqrt(n)` growth of the excess along `n = 2 mod 4` (HEURISTIC remark). The
   one-row costs of the Paley matrices at `n = 14, 18, 30` (`6, 6, >= 8`, against
   `r_n = 2.29, 2.79, 3.78`) are compatible with either.

## 6. The plus-I lift in quartic phase, and exhaustive minima

For symmetric zero-diagonal `+-1` `B` of order `m`, let
`K(B) = [[B, B+I],[B+I, -B]]` (`evidence/NOTE_2026-09-19_PLUS_I_COHERENT_LIFT.md`),
so `Q_K(u,v) = Q_B(u) - Q_B(v) + u^T B v + u.v` for `u, v in {+-1}^m`.

**Proposition 6.1 (quartic-phase identity; PROVED).** Define `z_k in mu_4 =
{1, i, -1, -i}` by `u_k + i v_k = (1+i) z_k` (a bijection `(u,v) <-> z`), put
`Q_B(z) = sum_{i<j} B_ij z_i z_j`, and let `r(z)` be the number of `k` with `z_k`
real. Then

    Q_K(u, v) = 2 ( Re Q_B(z) - Im Q_B(z) ) + ( 2 r(z) - m ).

*Proof.* With `w = u + i v`, symmetry and zero diagonal of `B` give
`Q_B(w) = Q_B(u) - Q_B(v) + i u^T B v`. Also `w_i w_j = (1+i)^2 z_i z_j =
2i z_i z_j`, so `Q_B(w) = 2i Q_B(z)`, i.e. `Q_B(u) - Q_B(v) = -2 Im Q_B(z)` and
`u^T B v = 2 Re Q_B(z)`. Finally `w_k^2 = 2 i u_k v_k` and `w_k^2 = 2i z_k^2`, so
`u_k v_k = z_k^2`, which is `+1` for real `z_k` and `-1` otherwise; thus
`u.v = 2r - m`. QED

**Corollary 6.2 (PROVED).** Put `Psi_1(B) = max_{z in mu_4^m} |Re Q_B(z) - Im
Q_B(z)|`. Since `|2r(z) - m| <= m`, Proposition 6.1 gives
`|Phi(K(B)) - 2 Psi_1(B)| <= m`. (The identity is exact for the plus-I block
itself; the real-valued bound `m_{2n} <= Psi(A) + n` of `HANDOFF.md` concerns the
block `[[A,A],[A,-A]]`.) **COMPUTED check:** the identity holds for every
`(u,v)` on 176 matrices `B` (all `B` of order 4, random `B` of orders 2 to 6),
`scripts/plus_i_quartic_phase_identity_check.py`.

**Exhaustive minima over all `B` (COMPUTED).** `K(DBD) = (D (+) D) K(B) (D (+) D)`
for diagonal `+-1` matrices `D`, so `Phi(K(B))` is invariant under switching of `B`
and `B` may be normalised to have first row and column `+1`. Exhaustive
enumeration of all `2^{(m-1)(m-2)/2}` normalised `B` with exact `Phi(K(B))`
(`scripts/plus_i_lift_all_b_exhaustive.c`), using that `Phi(K(B))` is congruent
to `C(2m,2)` modulo 2:

| `m` | `n = 2m` | parity of `Phi` | no `B` with `Phi(K(B)) <=` | minimum over all `B` | normalised `B` attaining it | `m_n` |
|---|---|---|---|---|---|---|
| 5 | 10 | odd | 11 | 13 | 12 | 13 |
| 6 | 12 | even | 16 | 18 | 12 | 18 |
| 7 | 14 | odd | 23 | 25 | 2292 | 21 |
| 8 | 16 | even | 30 | 32 | 3780 | `<= 30` |
| 9 | 18 | odd | 31 | 33 | 5040 | `<= 33` |

(Column six counts the normalised `B` with `Phi(K(B)) <=` the minimum; none lies
below it.)
The minima 13 and 18 at `n = 10, 12` equal `m_n`; the minima `25` and `32` exceed
`m_14 = 21` and the repository value `m_16 <= 30` (`HANDOFF.md`;
`evidence/NOTE_2026-09-11_GROWTH_CHAIN_M16.md`), so the plus-I lift over all `B`
is not optimal at `n = 14, 16`; at `n = 18` the minimum `33` equals the
repository bank value. Earlier records are single-source or restricted families:
`K(C_5) = 13`, the all-diagonal scans of
`NOTE_2026-09-11_ONE_VERTEX_EXTENSION_SCANS.md`, and `K(C_9) = 33`. The value 32
at `n = 16` is a minimum over lifts and is not an upper bound for `m_16`
(`m_16 <= 30` supersedes it).

## 7. The one-sided problem

**Proposition 7.1 (PROVED).** `min_A max_x Q_A(x) = floor(n/2)`.

*Proof.* Upper bound: for `A_ij = -1` off the diagonal, `Q_A(x) = -(s^2 - n)/2`
with `s = sum x_i`, whose maximum is `n/2` (`n` even, `s = 0`) or `(n-1)/2`
(`n` odd, `s = +-1`). Lower bound: for any `A` choose `x_1 = 1` and, for
`i = 2,...,n` in turn, `x_i = sign h_i` with `h_i = sum_{j<i} A_ij x_j` (any sign
if `h_i = 0`). Then `Q_A(x) = sum_i x_i h_i = sum_i |h_i|`. The sum `h_i` has
`i - 1` terms `+-1`, so `h_i` is odd for even `i` and `|h_i| >= 1`; hence
`Q_A(x) >= floor(n/2)`. QED

So the one-sided quantity is `Theta(n)`, and only the two-sided `Phi` has order
`n^{3/2}`.

## 8. Gaussian superadditivity and Fekete counterexamples

Let `g_ij` (`i < j`) be i.i.d. `N(0,1)`, `X_n(sigma) = n^{-1/2} sum_{i<j} g_ij
sigma_i sigma_j` on `{+-1}^n`, and `E_n = E max_sigma X_n(sigma)`.

**Proposition 8.1 (PROVED).** `E_{a+b} >= E_a + E_b` for all `a, b >= 1`. Hence
`lim E_n / n` exists, equals `sup E_n / n`, and is at most `sqrt(ln 2)`.

*Proof.* Let `n = a + b` and write `sigma = (s, t)` with `s in {+-1}^a`,
`t in {+-1}^b`. Put `Y(s,t) = X_a(s) + X_b'(t)` with `X_b'` an independent copy of
`X_b`; then `E max Y = E_a + E_b`, since the index set is the product of the two
index sets and the maximum of a sum of a function of `s` and a function of `t`
is the sum of the maxima. For `sigma, sigma'` put
`eps_i = sigma_i sigma'_i` and `R = sigma.sigma'/n`. Then
`Var(X_n(sigma) - X_n(sigma')) = n^{-1} sum_{i<j} (sigma_i sigma_j - sigma'_i
sigma'_j)^2 = n^{-1} [ n(n-1) - 2 sum_{i<j} eps_i eps_j ] = n^{-1} [ n^2 -
(sigma.sigma')^2 ] = n (1 - R^2)`. In the same way
`Var(Y(s,t) - Y(s',t')) = a (1 - R_1^2) + b (1 - R_2^2)` with `R_1 = s.s'/a`,
`R_2 = t.t'/b`, and `n R = a R_1 + b R_2`. By Cauchy-Schwarz
`(a R_1 + b R_2)^2 <= (a+b)(a R_1^2 + b R_2^2)`, i.e. `n R^2 <= a R_1^2 + b R_2^2`,
so `Var(X_n(sigma) - X_n(sigma')) >= Var(Y(s,t) - Y(s',t'))`. The Sudakov-Fernique
inequality for centred Gaussian processes on a finite set gives `E max X_n >=
E max Y`. Fekete's lemma applies to the superadditive sequence `E_n`; each
`X_n(sigma)` has variance `(n-1)/2`, so `E_n <= sqrt(2 n ln 2 (n-1)/2) <= n sqrt(ln 2)`.
QED

This is the zero-temperature Gaussian case of the Guerra-Toninelli argument (the
repository records the common-temperature interpolation in
`evidence/NOTE_2026-09-02_THERMODYNAMIC_INTERPOLATION_GATE.md`); the limit is the
Sherrington-Kirkpatrick ground-state energy `0.7632...` (Parisi formula,
Talagrand). The proof above compares increments for a fixed Gaussian disorder. For `m_n`
the coupling matrix is the minimisation variable with entries fixed to `+-1`, so
no disorder is averaged and no rescaling of couplings is available; that is the
mechanism recorded in the interpolation-gate note.

**Proposition 8.2 (Fekete failure for `f(n) = m_n / sqrt(n)`; PROVED from exact
values).** Both superadditivity and subadditivity of `f` fail:

    f(6) = 5/sqrt(6) = 2.0412 < 2 f(3) = 6/sqrt(3) = 3.4641,
    f(10) = 13/sqrt(10) = 4.1110 > 2 f(5) = 8/sqrt(5) = 3.5777.

*Proof.* Arithmetic with `m_3 = 3`, `m_5 = 4`, `m_6 = 5`, `m_10 = 13`
(`evidence/exact_m_table.json`, `N10_STRUCTURE.md`; `m_5, m_6` are also re-derived
in Section 10). QED

`f(n) ~ alpha_n n` has the same linear scale as `E_n`, so the contrast is
between a superadditive Gaussian sequence and a sequence for which neither
inequality holds. For `H(n) = m_n^{2/3}` the inequality `H(m+n) <= H(m) + H(n)`
fails at `2+2` and at `3+5` (also recorded in the public Marcus Miccelli claim
register, entry `composition.exact_subadditivity`). Any Fekete argument for `m_n`
requires an error term; the block cross term of `solution.md` Section 10 is of the
same order as `m_n`.

## 9. Saw-tooth in the exact data and the Paley-30 core

**Data (COMPUTED from the exact values).** For `n = 2,...,15`, `alpha_n` is
0.354, 0.577, 0.500, 0.358, 0.340, 0.486, 0.442, 0.444, 0.411, 0.466, 0.433,
0.427, 0.401, 0.465. The increments `m_{n+1} - m_n` for `n = 2,...,14` are 2, 1,
0, 1, 4, 1, 2, 1, 4, 1, 2, 1, 6. Local minima of `alpha_n` occur at `n = 5, 6, 10,
14` and local maxima at `n = 3, 7, 11, 15`: `n = 2 mod 4` (conference orders 6,
10, 14) is low and `n = 3 mod 4` high, with the large increments at
`n = 6, 10, 14`. By `solution.md` Proposition 6.1, applied to `n_k = a + rk`, every residue
class modulo `r` has the same `liminf` and `limsup` as the full sequence. Hence
`alpha_{4k+3} - alpha_{4k+2} >= A > 0` cannot hold for all large `k`: with
`L = liminf alpha_n`, one would have `liminf_k alpha_{4k+2} >= L` and so
`liminf_k alpha_{4k+3} >= L + A`, contradicting that the class `3 mod 4` also has
`liminf = L`. The pattern at `n <= 15` is therefore transient, and three data
points do not determine its trend.

**The Paley-30 core (COMPUTED).** Deleting one vertex from `C_30` leaves an order
29 matrix `B` (the Seidel matrix of the Paley graph on 29 points). Script
`scripts/one_vertex_tight_frame_paley.c 29 core`:

| `B` | `Phi` | `N0` | Gram of tight classes | Gram eigenvalues | one-row cost on `S` |
|---|---|---|---|---|---|
| `C_30` | 75 | 812 | `812 I` | all 812 | 8 |
| `B = C_30` minus a vertex | 72 | 406 | not a frame | 126 to 416 | 3 |

For `B` the Parseval bound of Proposition 1.1 is
`sqrt(29 * 126 / 406) = 3.000`, so it is attained: restoring the deleted vertex
gives `Phi = 75 = 72 + 3`. Hence the step `29 -> 30` through the core costs `3`
(against `r_29 = 3.756` if `m_29 = 72`, which is not known; `Phi(B) = 72` only
bounds `m_29` from above), while any extension of `C_30` to order 31 costs at
least `8`. The cost is low before the frame order and at least `sqrt(n)` at it, so
a signed sum of such steps can be small. The formulation (1)-(2) of
`evidence/SCOPE_FREEZE_2026-09-25_ONE_VERTEX_ONLY.md`, which extends a minimiser of
each order separately, does not see a step through a non-minimising core. Convergence needs only a bounded signed
chain sum (Proposition 3.2). Turning this compensation into a theorem requires a
handle on `m_{n+1}` at frame orders that is not available here (HEURISTIC).

## 10. Reproduction

All programs use only the C standard library or the Python standard library,
apart from `one_vertex_minimiser_extension_table.py`, which needs `numpy`. C
programs build with `gcc -O2 -fopenmp` (OpenMP optional); timings are for 8 cores.

| Script | Command | Content |
|---|---|---|
| `scripts/one_vertex_tight_frame_paley.c` | `./a.out 5`, `13`, `17`, `29`, `29 core` | Sections 1, 2, 9; `q = 29` takes about 20 s |
| `scripts/one_vertex_minimiser_extension_table.py` | `python3 ...` | Section 5, under 1 s |
| `scripts/normalized_exhaustive_phi_counts.c` | `./a.out n T` | exact `m_5,...,m_9` and counts, Section 0 |
| `scripts/plus_i_lift_all_b_exhaustive.c` | `./a.out m T` | Section 6; `m = 9` takes about 2.5 min |
| `scripts/plus_i_quartic_phase_identity_check.py` | `python3 ...` | Proposition 6.1 |
| `scripts/log_periodic_alpha_model_check.py` | `python3 ...` | Proposition 3.4 |

Outputs.

`one_vertex_tight_frame_paley` (`n = 6, 14, 18, 30, 29 core`):

    q=5 n=6 Phi=5 N0=12
    Gram of tight classes equals N0*I exactly: yes
    eigenvalues of Gram: min 12.0000 max 12.0000 (N0 = 12);  sqrt(n*min/N0) = 2.4495
    min_c max_x |x.c| = 4   (one-row cost over the tight classes)
    min_c sum_x |x.c| = 20   avg = 1.66667   avg/sqrt(n) = 0.6804
    q=13 n=14 Phi=21 N0=156
    Gram of tight classes equals N0*I exactly: yes
    eigenvalues of Gram: min 156.0000 max 156.0000 (N0 = 156);  sqrt(n*min/N0) = 3.7417
    min_c max_x |x.c| = 6   (one-row cost over the tight classes)
    min_c sum_x |x.c| = 392   avg = 2.51282   avg/sqrt(n) = 0.6716
    q=17 n=18 Phi=33 N0=204
    Gram of tight classes equals N0*I exactly: yes
    eigenvalues of Gram: min 204.0000 max 204.0000 (N0 = 204);  sqrt(n*min/N0) = 4.2426
    min_c max_x |x.c| = 6   (one-row cost over the tight classes)
    min_c sum_x |x.c| = 576   avg = 2.82353   avg/sqrt(n) = 0.6655
    q=29 n=30 Phi=75 N0=812
    Gram of tight classes equals N0*I exactly: yes
    eigenvalues of Gram: min 812.0000 max 812.0000 (N0 = 812);  sqrt(n*min/N0) = 5.4772
    min_c max_x |x.c| = 8   (one-row cost over the tight classes)
    min_c sum_x |x.c| = 2940   avg = 3.62069   avg/sqrt(n) = 0.6610
    q=29 n=29 (core) Phi=72 N0=406
    Gram of tight classes equals N0*I exactly: no
    eigenvalues of Gram: min 126.0000 max 416.0000 (N0 = 406);  sqrt(n*min/N0) = 3.0000
    min_c max_x |x.c| = 3   (one-row cost over the tight classes)
    min_c sum_x |x.c| = 1134   avg = 2.79310   avg/sqrt(n) = 0.5187

`normalized_exhaustive_phi_counts` (`n T`: `5 4`, `6 5`, `7 9`, `8 10`, `9 12`);
no value below `m_n` occurs in any run:

    n=5 Phi=4 normalised matrices: 12
    n=6 Phi=5 normalised matrices: 12
    n=7 Phi=9 normalised matrices: 3240
    n=8 Phi=10 normalised matrices: 4200
    n=9 Phi=12 normalised matrices: 1607760

`plus_i_lift_all_b_exhaustive` (`m T`):

    m=5 (n=10) T=11 count(Phi(K(B))<=T)=0 min Phi(K(B)) among them=-1
    m=5 (n=10) T=13 count(Phi(K(B))<=T)=12 min Phi(K(B)) among them=13
    m=6 (n=12) T=16 count(Phi(K(B))<=T)=0 min Phi(K(B)) among them=-1
    m=6 (n=12) T=18 count(Phi(K(B))<=T)=12 min Phi(K(B)) among them=18
    m=7 (n=14) T=23 count(Phi(K(B))<=T)=0 min Phi(K(B)) among them=-1
    m=7 (n=14) T=25 count(Phi(K(B))<=T)=2292 min Phi(K(B)) among them=25
    m=8 (n=16) T=30 count(Phi(K(B))<=T)=0 min Phi(K(B)) among them=-1
    m=8 (n=16) T=32 count(Phi(K(B))<=T)=3780 min Phi(K(B)) among them=32
    m=9 (n=18) T=31 count(Phi(K(B))<=T)=0 min Phi(K(B)) among them=-1
    m=9 (n=18) T=33 count(Phi(K(B))<=T)=5040 min Phi(K(B)) among them=33

`one_vertex_minimiser_extension_table.py`: the table of Section 5, preceded by one
line per class (`n`, class, `Phi`, `N0`, `ext`, `r_n`, `ext - r_n`).

`plus_i_quartic_phase_identity_check.py`:

    identity holds for all (u,v) on 176 matrices B: True

`log_periodic_alpha_model_check.py`:

    8 <= n < 2000000: m nondecreasing: True; m_(n+1)-m_n <= n-2 for n>=9: True
    alpha_n in [0.4111, 0.4498]
    max over n of (e_n - b*omega*sqrt(n)) = 0.993, so e_n <= 0.003 sqrt(n) + 1 on this range

## 11. Open points

The frame property for general `q`, for near-minimisers, and for any matrix not
of Paley type is unproved. The reach of block extensions with `d = o(n)` is
unresolved. Nothing here bears on the existence of `lim alpha_n`.
