# The Max_+ frame identity for every Paley order n = p^2 + 1, and the unconditional Max-Lipschitz bound

**Status:** one proved theorem with its corollaries, and exact finite computations.
Nothing here proves or refutes convergence of `alpha_n`; the original limit stays
OPEN and nothing was entered in `src/original_mo_status.py`. Every statement carries one
of three labels.

- **PROVED**: complete proof in this note.
- **COMPUTED**: exact finite computation; the script and its output are listed in Section 6.
- **HEURISTIC**: not proved.

## 0. Context and what is new

`solution.md` Proposition 15.27 states the identity

    (1/|Max_+|) sum_{y in Max_+} y y^T = I + C/p = 2 P_+                        (F)

for the Boolean eigenvectors `Max_+ = {y in {+1,-1}^n : C y = p y}` of the Paley conference
matrix `C` of order `n = p^2 + 1`, `p = sqrt(n-1)`, as an *assumption*, "certified for
Paley `n = p^2+1` at `p = 3, 5` (and holds at `p = 7`)". The Max-Lipschitz bound of the same
proposition, Proposition 15.20d (`k_star = O(n^{3/2})` implies `lim alpha_n = 1/2`) and the
`k_star = o(n^2)` form recorded in `evidence/E1_RIGIDITY_ATTACK.md` depend on (F). The note
`evidence/NOTE_2026-09-19_PALEY_GAP_ATTACK.md`, Section 3, and `HANDOFF.md` (2026-09-19)
record that the automorphism group `AGammaL` with square multiplications has a
2-dimensional commutant on `V_+` (`p = 3, 5, 7, 11`), so that Schur's lemma alone does not
give (F), which then "remains a finite certificate at `p = 3, 5, 7`".

This note proves (F) for every odd prime `p` (Theorem 2). The two free scalars of the
2-dimensional commutant are fixed by two linear constraints that every Boolean
distribution satisfies: `y_inf^2 = 1` and `trace = n`. Consequently the Max-Lipschitz
bound and Proposition 15.20d no longer depend on the unproved frame hypothesis. The
remaining hypothesis of those statements (a bound on the switching distance `k_star` of
minimisers from `C`) is untouched.

## 1. Setup

Let `p` be an odd prime, `q = p^2`, `F = GF(q)`, `chi` its quadratic character (`chi(-1) = 1`
since `q = 1 mod 4`), `Tr : F -> F_p` the absolute trace, `w = exp(2 pi i / p)` and
`psi_a(x) = w^{Tr(ax)}`. Index `R^n`, `n = q+1`, by `F` together with a point `inf`. The
symmetric Paley conference matrix is `C_{xy} = chi(x-y)` (`x != y` in `F`), `C_{xx} = 0`,
`C_{inf,x} = C_{x,inf} = 1`, `C_{inf,inf} = 0`; it satisfies `C^2 = q I`. Put `p = sqrt(q)`,
`P_+ = (I + C/p)/2`, `V_+ = range(P_+)` (the `+p`-eigenspace, dimension `n/2`), and
`Q_A(x) = sum_{i<j} A_ij x_i x_j`.

`H` denotes the group of permutations `x -> s x + b` of `F` (`s` a nonzero square, `b in F`)
extended by `inf -> inf`, of order `q(q-1)/2`. Each `g in H` satisfies `C_{g x, g y} = C_{xy}`,
so the permutation matrix `P_g` commutes with `C`, preserves `V_+` and preserves
`Max_+`.

**Gauss sum.** `G = sum_z chi(z) psi_1(z)` satisfies `G^2 = chi(-1) q = q` and `G` is real
(substitute `z -> -z`), so `G = eps p` with `eps in {+1,-1}`. For `a != 0`,
`(C (psi_a, 0))_x = sum_z chi(z) psi_a(x-z) = chi(a) G psi_a(x)` and `(C (psi_a,0))_inf = sum_y psi_a(y) = 0`.
Hence `(psi_a, 0)` is an eigenvector with eigenvalue `chi(a) G = chi(a) eps p`. Also
`u_0 = (1,...,1; p)` satisfies `C u_0 = p u_0`. Let `S = {a != 0 : chi(a) = eps}`, a coset of the
squares, `|S| = (q-1)/2`.

**Lemma 1 (PROVED).** `V_+ = span{u_0} (+) span{(psi_a, 0) : a in S}`.
*Proof.* The listed vectors are `+p`-eigenvectors, linearly independent, and there are
`1 + (q-1)/2 = n/2 = dim V_+` of them (`C^2 = qI` and `tr C = 0` give both eigenvalues
multiplicity `n/2`). QED

**Lemma 2 (PROVED; a Boolean eigenvector exists).** Choose `lambda in F^*` with
`chi(lambda) = eps` and `h : F_p -> {+1,-1}` taking the value `+1` exactly `(p+1)/2` times.
Put `f(x) = h(Tr(lambda x))` and `y = (f ; sum_t h(t))`. Then `y in {+1,-1}^n` and `C y = p y`.
*Proof.* `sum_t h(t) = +-1` since `p` is odd, so `y` is Boolean. Expanding `h` on `F_p` in
additive characters, `f = sum_{c in F_p} hhat(c) psi_{c lambda}` with `hhat(0) = (sum_t h(t))/p`. Every
element of `F_p` is a square in `F`, so `chi(c lambda) = chi(lambda) = eps` for `c in F_p^*`,
i.e. `c lambda in S`. Therefore
`y = hhat(0) (1,...,1; p) + sum_{c != 0} hhat(c) (psi_{c lambda}, 0)`
because `p hhat(0) = sum_t h(t)`; by Lemma 1 `y in V_+`. QED

(Lemma 2 reproves `Phi(C) = n p / 2` of `CORE.md` Section 8.)

## 2. The frame identity

**Theorem 2 (PROVED).** Let `mu` be any probability measure on `Max_+` that is invariant
under `H`, for instance the uniform measure on `Max_+`, or the uniform measure on the `H`-orbit
of the vector `y` of Lemma 2. Then

    E_mu [ y y^T ] = I + C/p = 2 P_+ .

*Proof.* Let `M = E_mu[y y^T]`. It is real symmetric and positive semidefinite, `P_g M P_g^T = M`
for `g in H`, `M = P_+ M P_+` (every `y` lies in `V_+`), and `M_ii = 1` for all `i` (Boolean `y`).
Let `P_0` be the orthogonal projection onto `span{u_0}` and `P_W` the orthogonal projection onto
`W = span{(psi_a, 0) : a in S}`, so `P_+ = P_0 + P_W` by Lemma 1. The vectors `u_0` and
`(psi_a, 0)` (`a in S`) are pairwise orthogonal, `|u_0|^2 = q + p^2 = 2q`, `|(psi_a,0)|^2 = q`.

1. *Translations.* For `tau_c : x -> x + c`, `P_{tau_c} (psi_b, 0) = w^{-Tr(bc)} (psi_b, 0)` and
   `P_{tau_c} u_0 = u_0`. Invariance of `M` gives, for `a, b in S`,
   `<psi_a, M psi_b> = w^{Tr((a-b)c)} <psi_a, M psi_b>` for every `c`. The trace form is nondegenerate, so for
   `a != b` some `c` has `Tr((a-b)c) != 0` and the entry vanishes. The same argument with the fixed vector `u_0`
   gives `<u_0, M psi_b> = 0` for `b in S`. Hence `M` is diagonal in the orthogonal basis
   `{u_0} u {(psi_a,0) : a in S}` of `V_+`.
2. *Square multiplications.* `P_s (psi_a, 0) = (psi_{a/s}, 0)` and `P_s u_0 = u_0`. The squares act transitively on
   the coset `S`, so the diagonal entries `<psi_a, M psi_a>/q`, `a in S`, are one constant `lambda_2`. Writing
   `lambda_1 = <u_0, M u_0>/(2q)`, steps 1 and 2 together with `M = P_+ M P_+` give
   `M = lambda_1 P_0 + lambda_2 P_W`.
3. *Two constraints.* The vectors `(psi_a, 0)` vanish at `inf`, while `u_0(inf) = p`, so
   `(P_W)_{inf,inf} = 0` and `(P_0)_{inf,inf} = p^2/(2q) = 1/2`. Since `M_{inf,inf} = 1`,
   `1 = lambda_1 / 2`, i.e. `lambda_1 = 2`. Then `n = tr M = lambda_1 + lambda_2 dim W = 2 + lambda_2 (q-1)/2`, and
   `n = q + 1` gives `lambda_2 = 2`.

Hence `M = 2 (P_0 + P_W) = 2 P_+ = I + C/p`. QED

**Remark.** The commutant of `H` on `V_+` is exactly the two scalars `(lambda_1, lambda_2)`
(dimension 2, as recorded in `NOTE_2026-09-19_PALEY_GAP_ATTACK.md`, Section 3). Booleanity
enters only through `M_ii = 1`, and the single coordinate `inf`, which `u_0` sees and the
`psi_a` do not, separates the two scalars.

**Corollary 3 (Max_- , switching; PROVED).** For `g` a non-square, `D = diag(-1 at inf, +1 else)`:
`(D P_g) C (D P_g)^T = -C`, so `D P_g` maps `Max_+` bijectively onto `Max_- = {y : C y = -p y}`; the
uniform measure on `Max_-` has second-moment matrix `I - C/p`. For a diagonal sign matrix
`D'`, the matrix `C' = D' C D'` has `Max_+(C') = D' Max_+(C)` and second-moment matrix
`I + C'/p`.

**Corollary 4 (cut polytope; PROVED).** The matrix `I + C/p` is the second-moment matrix of a probability
distribution on `{+1,-1}^n`, i.e. it lies in the cut polytope (correlation polytope of
Boolean vectors); with the orbit measure the support has at most `q(q-1)/2` points
(`p(p^2-1)/4` for `p = 3, 5, 7, 11, 13`, COMPUTED, Section 4).

## 3. Consequences

**Corollary 5 (correlation lower bound; PROVED).** For every `A in S_n`, `n = p^2+1`, and every
diagonal sign matrix `D'`,

    Phi(A) >= | sum_{i<j} A_ij (D' C D')_ij | / p .                            (1)

*Proof.* Let `y` be uniform on `Max_+(C)`. `Q_A(D' y) = sum_{i<j} A_ij D'_i D'_j y_i y_j`, and
`E y_i y_j = C_ij / p` by Theorem 2. So `E Q_A(D' y) = <A, D'CD'>/p`, and `Phi(A) >= |E Q_A(D' y)|`. QED

**Corollary 6 (Max-Lipschitz, unconditional; PROVED).** Let `N = n(n-1)/2 = n p^2/2` and let
`k` be the Hamming distance (number of differing pairs `i<j`) between `A` and some
`D' C D'`. Then `<A, D'CD'> = N - 2k` and `Phi(C) = n p / 2 = N/p`, so

    Phi(A) >= Phi(C) - 2k/p .                                                 (2)

Every value of `Q` has the parity of `N`, so `Phi(A) = Phi(C) (mod 2)` and

    Phi(A) >= Phi(C) - 2 floor(k/p) .                                         (3)

This is the statement of `solution.md` Proposition 15.27 (2) with the frame hypothesis
discharged. In particular a perturbation of `C` by `k < p` entries, in any switching class,
does not lower `Phi` at all, and `k = O(sqrt n)` flips lower it by at most `O(1)`.
For comparison the edge-Lipschitz bound `Phi(A) >= Phi(C) - 2k` (Proposition 15.20b) is weaker by the
factor `p`. Values for the known undercuts:

| `n` | `p` | `Phi(C)` | best known `Phi` | gain | bound (3): `k >=` | known switching distance `k` |
|---|---|---|---|---|---|---|
| 10 | 3 | 15 | 13 (exact `m_10`) | 2 | 3 | 5 (`solution.md`, N10 analysis) |
| 26 | 5 | 65 | 61 | 4 | 10 | 122 (`CORE.md` Section 8) |
| 50 | 7 | 175 | 169 | 6 | 21 | not recorded |

**Corollary 7 (normalised distance; PROVED).** Write `gamma_p = (Phi(C) - m_n)/n^{3/2}`. Every
minimiser `A` of order `n = p^2+1` is at normalised switching distance
`k_star/N >= gamma_p sqrt(n/(n-1)) >= gamma_p` from `C`. Equivalently `m_n >= Phi(C) - 2 k_star/p` gives
`alpha_n >= (1/2) sqrt(1 - 1/n) - 2 k_star/(p n^{3/2})`, and `2 k_star/(p n^{3/2}) = (2 k_star/n^2)(1+O(1/n))`.

**Corollary 8 (conditional settlement; PROVED).** If `k_star(n_k) = o(n_k^2)` along `n_k = p_k^2+1`,
then `alpha_{n_k} -> 1/2` and, by `CORE.md` Section 8 (ratio-dense Paley orders), `lim alpha_n = 1/2`.
This removes the frame hypothesis from `solution.md` Proposition 15.20d and weakens its
hypothesis `k_star = O(n^{3/2})` to `k_star = o(n^2)`. The hypothesis on `k_star` is not
proved; the undercut at `n = 26` has `k/N = 122/325 = 0.375`, and (2) shows that any
undercut of size `c n^{3/2}` has `k/N >= c`.

## 4. Computations (COMPUTED)

Script `scripts/paley_square_max_plus_frame_check.py` (numpy only). It builds `GF(p^2)`, `C`, the vector `y`
of Lemma 2, the `H`-orbit of `y`, and compares the orbit average of `y y^T` with `I + C/p`
(exact; maximal entrywise deviation printed):

| `p` | `n` | `H`-orbit size of `y` | max abs deviation of orbit average from `I + C/p` |
|---|---|---|---|
| 3 | 10 | 6 | 0 |
| 5 | 26 | 30 | 0 |
| 7 | 50 | 84 | 0 |
| 11 | 122 | 330 | 0 |
| 13 | 170 | 546 | 0 |

The orbit size equals `p(p^2-1)/4` in all five cases (observation, not proved). With `--full`, all
Boolean vectors of `V_+` are enumerated through a pivot basis: `p = 3`: 12 vectors, `p = 5`: 260
vectors (130 classes), `p = 7`: 11452 vectors (5726 classes); the average of `y y^T` over each full set equals `I + C/p`
(deviation 0), in agreement with the earlier certificates and with Theorem 2. For
`p = 5` the tight classes of `C` are the 260 classes `Max_+ u Max_-` (130 each), each class
satisfying `C y = +-p y`. Each check of `E Q_A = <A,C>/p` for random `A` passes.

## 5. Prime orders: the signed tight covariance (COMPUTED, not proved)

For prime `q = 1 mod 4` there is no Boolean eigenvector and `Phi(C) < n p / 2`. Let `S` be the set of tight
classes (`|Q_C| = Phi(C)`), `sigma(x) = sign Q_C(x)`, `N0 = |S|`. The program
`scripts/paley_signed_tight_covariance.c` computes in integer arithmetic

    sum_{x in S} sigma(x) x x^T = (N0 Phi(C) / N) C ,       sum_{x in S} x x^T = N0 I .

| `q` | `n` | `Phi(C)` | `N0` | Gram `= N0 I` | signed matrix `= (N0 Phi/N) C` |
|---|---|---|---|---|---|
| 5 | 6 | 5 | 12 | yes | yes |
| 13 | 14 | 21 | 156 | yes | yes |
| 17 | 18 | 33 | 204 | yes | yes |
| 29 | 30 | 75 | 812 | yes | yes |
| 37 | 38 | 109 | 2812 | yes | yes |

Each row is the output of one run of the program (Gray-code enumeration of all `2^(n-1)` classes; `q = 37` takes about
6 minutes on 8 cores). `Phi(C)` is the exact maximum of `|Q_C|` over all classes, so `Phi(C_38) = 109` is exact, with
`N0 = 2812` tight classes (1406 with each sign; every tight set listed is split evenly between the two signs).

**Consequence (PROVED given the computed identity at these orders).** For `n in {6, 14, 18, 30, 38}` and every
`A in S_n`, with `k` the switching distance of `A` from `C`, `N = n(n-1)/2`:
`Phi(A) >= |E_S sigma Q_A| = (Phi(C)/N) |<A,C>| = Phi(C)(1 - 2k/N)`. The gain is therefore at most
`2 Phi(C) k / N`, with the parity refinement `Phi(A) = Phi(C) (mod 2)`:

| `n` | `Phi(C)` | `2 Phi/N` (gain per flip, upper bound) | least `k` that the bound allows for a gain of 2 |
|---|---|---|---|
| 14 | 21 | 0.462 | 5 |
| 18 | 33 | 0.431 | 5 |
| 30 | 75 | 0.345 | 6 |
| 38 | 109 | 0.310 | 7 |

(`Phi(C_14) = 21 = m_14`, so no gain exists at `n = 14`.) At squares the same quantity is `2 Phi/N = 2/p`.

**HEURISTIC (open).** The identity `E_S[sigma x x^T] = (Phi/N) C` for all prime `q = 1 mod 4` would follow if the
signed automorphism group of `C` (including the signed inversion `x -> -1/x`) acts irreducibly on each of `V_+`,
`V_-` with `V_+` not isomorphic to `V_-` (dimension `(q+1)/2` principal-series type); the argument of Theorem 2 does not
apply because no Boolean vector lies in `V_+`. Only the five computed orders are claimed.

## 6. Reproduction

| Script | Command | Content |
|---|---|---|
| `scripts/paley_square_max_plus_frame_check.py` | `python3 ... 3 5 7 11 13` ; `python3 ... --full 3 5 7` | Section 4 (numpy only; `--full` for `p=7` takes about 1 minute) |
| `scripts/paley_signed_tight_covariance.c` | `gcc -O2 -fopenmp -o a.out ...; ./a.out q` for `q = 5, 13, 17, 29, 37` | Section 5 (`29` about 30 s on 8 cores; `37` about 6 minutes) |

Output of `paley_square_max_plus_frame_check.py 3 5 7 11 13`:

    p=3 n=10 orbit size=6 (p(p^2-1)/4=6) max|M-(I+C/p)|=0.0e+00 y_inf=1 chi(lambda)=1
    p=5 n=26 orbit size=30 (p(p^2-1)/4=30) max|M-(I+C/p)|=0.0e+00 y_inf=1 chi(lambda)=-1
    p=7 n=50 orbit size=84 (p(p^2-1)/4=84) max|M-(I+C/p)|=0.0e+00 y_inf=1 chi(lambda)=1
    p=11 n=122 orbit size=330 (p(p^2-1)/4=330) max|M-(I+C/p)|=0.0e+00 y_inf=1 chi(lambda)=1
    p=13 n=170 orbit size=546 (p(p^2-1)/4=546) max|M-(I+C/p)|=0.0e+00 y_inf=1 chi(lambda)=-1

Output with `--full 3 5` and `--full 7` (appended to the corresponding line): `all Boolean vectors in V_+: 12`,
`260`, `11452` (including sign pairs), deviation `0.0e+00`.

Output of `paley_signed_tight_covariance`:

    q=5 n=6 Phi=5 N0=12  Gram==N0*I: yes   signed M==(N0*Phi/N)*C: yes
    q=13 n=14 Phi=21 N0=156  Gram==N0*I: yes   signed M==(N0*Phi/N)*C: yes
    q=17 n=18 Phi=33 N0=204  Gram==N0*I: yes   signed M==(N0*Phi/N)*C: yes
    q=29 n=30 Phi=75 N0=812  Gram==N0*I: yes   signed M==(N0*Phi/N)*C: yes
    q=37 n=38 Phi=109 N0=2812  Gram==N0*I: yes   signed M==(N0*Phi/N)*C: yes

## 7. Open points

- The hypothesis of Corollaries 7 and 8 (`k_star = o(n^2)` for minimisers at `n = p^2+1`) is not proved; the known undercut
  at `n = 26` has `k/N = 0.375`.
- Theorem 2 gives `E y y^T` for `H`-invariant measures on `Max_+`. It carries no information on higher moments of `Max_+`
  (for instance `E (y.c)^4`), which control the block-extension bound `beta` of `NOTE_2026-10-11_FRAME_OBSTRUCTION_SINGLE_ROW_EXTENSION.md`.
- The signed tight covariance at prime orders (Section 5) is a finite observation.
- The bound (1) is sharp only at `A = D' C D'`; it gives `Phi(A) >= 0` and no information for `A` at normalised distance `1/2` from every
  switched Paley matrix.
