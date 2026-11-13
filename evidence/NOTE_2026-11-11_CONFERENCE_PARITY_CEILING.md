# Conference parity ceiling: Φ(C54) = 189, Φ(C102) = 501, Φ(C62) = 229, and no uniform Paley gap

**Status:** draft session note (2026-11-11), independent review pending.
Lemma 1, Theorem 2 and Theorem 4 are proved below (elementary). §3 holds
exact finite certificates. §5 is witnesses only. The original MO limit is
OPEN. Nothing here bounds any `m_n` from below, nothing is a cross-order
statement about `m_n`, and nothing is entered in `src/original_mo_status.py`.

**Implications established** (the reason this note exists):

1. `Φ(C54) = 189`. This closes the bracket `{189, 191, 193, 195}` left open in
   `NOTE_2026-11-11_C54_LOWER_BOUND.md`, with no enumeration.
2. `Φ(C102) = 501` and `Φ(C62) = 229`, both beyond exhaustive reach. `C62` is
   the first Paley matrix proved to sit strictly below the ceiling of
   Theorem 2 (235).
3. The prime Paley family has **no** uniform gap below the spectral ceiling:
   `ρ(C_{q+1}) ≥ b/√q` for `q = a² + b²`, `b` odd (Theorem 4). The target
   "`ρ ≤ 1 − η` uniformly on the ratio-dense prime subfamily" in
   `NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md` is false
   for the full prime family, and statement (U) of
   `NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md` is false for Paley.

**Reused results.** The switched degree-parity lemma of Proposition 15.762
(`NOTE_2026-09-03_CONFERENCE_CUBE_GAP.md`), proved there for order `p² + 1`;
the Momihara–Suda excess bound quoted in the 08-29 note; the order-54
witness; `c42_result.json`. **New here:** the lemma for every `q ≡ 1 (mod 4)`,
the resulting ceiling at nonsquare `q`, its attainment, the torus witness for
the quadratic form, and the defect-vector decision procedure.

Notation: `C` a symmetric conference matrix of order `n = q + 1`
(`C = Cᵀ`, zero diagonal, `±1` elsewhere, `C² = qI`), `Q_C(x) = ½ xᵀCx`,
`S = ½ n√q`, `ρ(C) = Φ(C)/S`. For Boolean `x` the **local fields** are
`h_i = x_i (Cx)_i`, so `Σ h_i = 2 Q_C(x)`.

## 1. Parity lemma

**Lemma 1.** For every symmetric conference matrix `C` and every Boolean `x`,
all `h_i` are odd and congruent to each other modulo 4. Flipping one
coordinate of `x` moves every field to the other class.

*Proof.* Put `D = diag(x) C diag(x)`; it is a symmetric conference matrix
with row sums `h_i`. For `i ≠ j`, `0 = (D²)_ij = Σ_{k≠i,j} D_ik D_jk`. For
`u, v ∈ {±1}`, `uv ≡ u + v − 1 (mod 4)`. Summing over the `n − 2` values of
`k`,

    0 ≡ (h_i − D_ij) + (h_j − D_ji) − (n − 2) = h_i + h_j − 2 D_ij − (n − 2)   (mod 4).

`2 D_ij ≡ 2` and `n ≡ 2 (mod 4)`, so `h_i + h_j ≡ 2 (mod 4)`; both are odd,
hence `h_i ≡ h_j (mod 4)`. Flipping `x_k` changes `h_i` by `∓2 D_ik` for
`i ≠ k` and negates `h_k`; each is a shift by `2 (mod 4)`. ∎

Proposition 15.762 states this for `n = p² + 1` using `8 | n − 2`; only
`4 | n − 2` is needed.

## 2. The ceiling

Two more identities hold for every `x`:

    (B)  Σ h_i² = |Cx|² = q n,          (R)  D h = D² 1 = q 1.

**Theorem 2.** For odd `κ > 0` put `F(κ) = n (q + κ² − 4) / (4κ)`. If the
fields of `x` are `≡ κ + 2 (mod 4)` then `Q_C(x) ≤ F(κ)`. Hence, for every
symmetric conference matrix of order `q + 1`,

    Φ(C)  ≤  max_{r ∈ {1,3}}  min_{κ ≡ r (mod 4)}  F(κ)   <   S      (q not an odd square).

Equality `Q_C(x) = F(κ)` holds iff every field is `κ − 2` or `κ + 2`, i.e.
`Cx = κx + 2y` with `y` Boolean; then `Cy = ½(q − κ²) x − κ y`.

*Proof.* With `c = κ − 2`, every `h_i ≡ c (mod 4)`, so
`(h_i − c)(h_i − c − 4) ≥ 0`. Sum over `i` and use (B):
`qn − (2c + 4)·2Q + c(c + 4) n ≥ 0`, which is `Q ≤ F(κ)`. Apply the same to
`−C` for the minimum. `F(κ) < S` unless `κ = √q ± 2`. ∎

**Integral refinement.** Write `h_i = c + 4 t_i`, `T₁ = Σ t_i`, `T₂ = Σ t_i²`.
Then (B) reads `16 T₂ + 8c T₁ = n (q − c²)`, with `T₂ ≥ |T₁|` and
`T₂ ≡ T₁ (mod 2)`, and

    Q = (n c + 4 T₁)/2 = n (q + c²)/(4c) − 4 T₂ / c .

So `Q` takes few values near the top, each with a short list of admissible
field multisets, and (R) removes some of those: for each `i` there must be
`(q − h_i)/2` indices `j ≠ i` with `Σ h_j = (Σh − h_i − q)/2`.
`scripts/parity_ceiling/ceiling.py` computes the closed form, the integral
ceiling and the ceiling with (R). All three bound `Φ(C)` for **every**
symmetric conference matrix of that order.

**Relation to the excess bound.** The Momihara–Suda bound allows two adjacent
odd row sums `k, k + 2` because it treats `yᵀCx` with independent row and
column switching. For `y = x` Lemma 1 forbids adjacent odd values, which is
where the `4` in `F` comes from: `S − F = n (4 − (κ − √q)²)/(4κ)` can reach
`3n/(4√q)`, against at most `n/(4√q)` for the excess bound. The ceiling still
satisfies `1 − F/S = O(1/q)`, so, as the 08-29 note says of integrality in
general, it cannot produce a fixed gap.

## 3. Table and exact values

`spectral` is `S`; `closed` is Theorem 2 rounded down to an odd integer;
`ceiling` uses the integral refinement and (R).

      q    n  spectral  closed  ceiling  Paley value
      5    6     6.71       5       5    5
      9   10    15.00      15      15    15
     13   14    25.24      23      21    21
     17   18    37.11      33      33    33
     25   26    65.00      65      65    65
     29   30    80.78      79      75    75
     37   38   115.57     111     109    109
     41   42   134.47     129     129    129      (also the 2^41 enumeration)
     49   50   175.00     175     175    175
     53   54   196.56     195     189    189      NEW exact
     61   62   242.12     237     235    229      NEW exact, below the ceiling
     73   74   316.13     311     305    >= 297
     89   90   424.53     421     417    >= 401
     97   98   482.59     475     473    >= 457
    101  102   512.54     505     501    501      NEW exact
    109  110   574.22     567     565    >= 555
    113  114   605.92     601     599    >= 591
    173  174  1144.31    1141    1133    >= 1131  (Theorem 4)
    229  230  1740.27    1737    1731    >= 1725  (Theorem 4)
    293  294  2516.23    2513    2509    >= 2499  (Theorem 4)

At every Paley order up to 54, and at 102, the Paley matrix attains the
ceiling. So at those orders **no symmetric conference matrix has larger Φ than
Paley**.

**Order 54 by hand** (`q n = 2862`).
Class `h ≡ 1`: `h = 5 + 4t` gives `5T₁ + 2T₂ = 189`; `T₂ ≥ T₁` forces
`T₁ ≤ 27`, so `Q = 135 + 2T₁ ≤ 189`.
Class `h ≡ 3`: `h = 7 + 4t` gives `7T₁ + 2T₂ = 27`, `T₁` odd, `T₁ ≤ 3`,
`Q = 189 + 2T₁`. `T₁ = 1` needs `T₂ = 10`, violating `T₂ ≡ T₁ (mod 2)`.
`T₁ = 3` needs `T₂ = 3`: field 11 on three points `I`, 7 elsewhere; then (R)
at `i ∈ I` reads `53 = 7·11 + 4 Σ_{j∈I∖i} D_ij ∈ {69, 77, 85}`, impossible. So
`Q ≤ 187` in this class.
Hence `Q_C(x) ≤ 189` for every `x` and every `C` of order 54, likewise for
`−C`, and the stored witness gives `Φ(C54) = 189`. In particular 191 and 193
are not admissible values at all.

**Order 102 by hand** (`q n = 10302`).
Class `h ≡ 1`: `h = 9 + 4t` gives `9T₁ + 2T₂ = 255`, `T₁` odd, `T₁ ≤ 23`;
`T₁ = 23` needs `T₂ = 24` (parity), so `T₁ ≤ 21` and `Q = 459 + 2T₁ ≤ 501`.
Class `h ≡ 3`: `h = 11 + 4t` gives `11T₁ + 2T₂ = −255`, so `T₁ ≤ −29`;
`T₁ = −29` needs `T₂ = 32` (parity), so `Q = 561 + 2T₁ ≤ 499`.
The witness in `witnesses.json` has `Q = 501`. Hence `Φ(C102) = 501`.

**Order 62** (`q = 61`). The admissible values from the top are
`235, 229, 227, 221, …`; 231 and 233 do not occur. At 235 integrality leaves
six field multisets. Four violate (R); the other two are

    {3: 2, 7: 50, 11: 9, 15: 1}      and      {3: 3, 7: 47, 11: 12}.

*Reduction.* For modal value `c` put `e = (Cx − cx)/4`, so `e_i = t_i x_i`,
and `m = (q − c²)/4`. Then `(C + c) e = m x`; conversely an integer `e` with
`(C + c) e = m x`, `x` Boolean, `e_i = t_i x_i`, gives `Cx = cx + 4e`. So the
pattern occurs iff there are a support `S`, values `t_i` and signs `x_i` on
`S` with

    (Ce)_i = x_i (m − c t_i)  on S,          |(Ce)_i| = m  off S.

Here `c = 7`, `m = 3`, `|S| = 12` or `15`: a search over sparse integer
vectors, not over `2^61` states. `PSL(2,61)` acts 2-transitively by signed
permutations preserving `C` and `x ↦ −x` is a symmetry, so the first support
vertex is pinned to `∞` with `x = +1`, the second to `0`, and the third to
`1` or the least nonsquare (the two orbits of the split torus).
`scripts/parity_ceiling/defect.c` runs the depth-first search with exact
reachability pruning on every row.

    pattern {3:2, 7:50, 11:9, 15:1}      224,280 nodes        0 solutions
    pattern {3:3, 7:47, 11:12}       520,107,727 nodes        0 solutions
    the four patterns failing (R)    rejected at the root  0 solutions

Since `C62` is equivalent to `−C62`, `|Q_{C62}(x)| ≤ 229` for all `x`; the
witness gives **`Φ(C62) = 229`**, hence `m_62 ≤ 229` and
`229 / 62^{3/2} = 0.4691`.

*Checks on the solver.* Against brute force over all states at orders 14 and
18 it was asked about every multiset admissible by Lemma 1 and (B): 245
multisets, 89 of which occur, zero disagreements. It finds the known
maximiser patterns at orders 14, 18, 30, 38, 42. Receipt:
`c62_exclusion.json`. Independent evidence at order 62: every state with a
nontrivial stabiliser in `PSL(2,61)` has `|Q| ≤ 229` (the four orbit quotients
for the primes 2, 3, 5, 31 enumerated exactly; the unipotent case gives 61),
and 1,030,000 annealing restarts reach 229 in 252,869 runs and nothing higher.

**Explicit upper bounds this gives** (draft, per-order): `m_54 ≤ 189`,
`m_62 ≤ 229`, `m_102 ≤ 501`, and from the ceiling alone `m_74 ≤ 305`,
`m_90 ≤ 417`, `m_98 ≤ 473`, `m_110 ≤ 565`, `m_114 ≤ 599`.

## 4. The torus witness: no uniform gap for prime Paley

**Theorem 4.** Let `q ≡ 1 (mod 4)` be prime, `q = a² + b²` with `b` odd, and
`C` the Paley matrix of order `n = q + 1`. There is a Boolean `x` with

    xᵀ C x = (q + 1) b,        fields in {b − a, b + a}.

Hence `Φ(C_{q+1}) ≥ (q + 1) b / 2`, i.e. `ρ(C_{q+1}) ≥ b/√q`.

*Proof.* With `v_∞ = (1,0)`, `v_j = (j,1)` we have
`C_PQ = χ(det(v_P, v_Q))`. For `g ∈ SL(2,q)` write `g v_P = λ_P(g) v_{gP}`;
then `(g·x)(gP) = χ(λ_P(g)) x(P)` is an action by signed permutations that
preserves `xᵀCx`, and `−I` acts trivially. Let `T` be the image of the
norm-one subgroup of `GF(q²)*` acting on `GF(q)² ≅ GF(q²)`. It is cyclic of
order `n/2` and acts freely on the `n` points, so there are two orbits and
each carries one invariant `±1` vector `u₁, u₂` up to sign. `C` preserves
`V = span(u₁, u₂)`, and in this basis `C|_V = [[α, β], [β, δ]]` with integer
entries, symmetric because `|u₁| = |u₂|`. `C² = q` gives `α² + β² = q`,
`β(α + δ) = 0`; `q` is not a square, so `β ≠ 0` and `δ = −α`. `α` is a sum of
`(q − 1)/2` signs, hence even, so `|β| = b`. Take `x = u₁ ± u₂` with the sign
of `β`: `xᵀCx = (n/2)(α + 2|β| + δ) = n b`, and the fields are `b ± α`. ∎

Explicitly `x(∞) = 1`, `x(j) = φ(j² − d)`, `d` a nonsquare, `g` a primitive
root, `φ(g^k) = +1` for `k mod 4` in `{0,1}` (or in `{0,3}`; one choice gives
`+(q+1)b`, the other `−(q+1)b`). `scripts/parity_ceiling/torus.py`
verifies `xᵀCx = (q+1)b` in exact arithmetic for all 80 primes
`q ≡ 1 (mod 4)` below 1000.

**Consequences.**

- For `q = b² + 4` the witness is the equality case `κ = b` of Theorem 2, and
  `S − b(q + 1)/2 = 2(q + 1)/(√q + b) = √q + O(q^{-1/2})`. Proved values: `Φ(C174) ≥ 1131`
  (`ρ ≥ 0.9884`), `Φ(C230) ≥ 1725`, `Φ(C294) ≥ 2499`, `Φ(C734) ≥ 9909`
  (`ρ ≥ 0.9973`). With the ceiling, `Φ(C174) ∈ {1131, 1133}`.
- 11 of the 80 primes below 1000 have `b/√q` above `√15/4 = 0.9682`, the
  first-moment cap of the Haar model. The random-orthogonal reading of Paley
  fails at those primes by proof, not by parity step.
- By Hecke's equidistribution of the arguments of Gaussian primes there are
  infinitely many primes with `a/b` arbitrarily small, so
  `lim sup ρ(C_{q+1}) = 1` over primes, unconditionally. No constant `η > 0`
  works for the whole prime Paley family.
- Theorem 4 is the quadratic-form counterpart of Momihara–Suda's Theorem 1.2
  (bilinear excess at `q = 4m² + 1`, row sums `2m ± 1`): the bilinear problem
  sees the even component `a`, the quadratic form sees the odd component `b`.
  (Their statements were read from the arXiv listing and a machine summary of
  the PDF, not line by line; recheck before citing in a paper.)

## 5. Where Paley leaves the ceiling (order 62 exact, the rest witnesses)

`best` is the largest value found by annealing on the full cube
(`sa.c`, receipts in `witnesses.json`). `sym` is the best state with a
nontrivial stabiliser in `PSL(2,q)`, from the orbit quotients of every prime
order (`sym.py`; exact Gray code up to 34 orbits, annealing on larger
quotients). `|Stab|` is the exact stabiliser order of the stored best state.

      q    n  ceiling   best    sym   |Stab|  Phi/n^1.5   restarts (hits at best)
     61   62    235     229    229      6     0.4691     1,030,000 (252,869)   exact, §3
     73   74    305     297    297      2     0.4666       180,000 (22,599)
     89   90    417     401    401      2     0.4697       120,000 (6,233)
     97   98    473     457    457      2     0.4711        20,000 (465)
    101  102    501     501    501     50     0.4863        20,000 (66)        exact, §3
    109  110    565     555    555     10     0.4811        20,000 (187)
    113  114    599     591    591     24     0.4855        20,000 (51)

For comparison, the maximisers at the ceiling orders have stabilisers of order
14 (`q = 13`), 24 (17), 30 (29), 18 (37), 60 (41), 54 (53): the torus
normaliser, `S₄`, `A₅`, or a dihedral group, with at most four orbits.

Readings, none of them a theorem beyond order 62:

- Up to order 54 every Paley value is the ceiling, reached by a state with a
  large stabiliser. That is why the Haar band of the ROM note failed: those
  values are arithmetic.
- From order 62 the two separate. At 62 (exact), and at 74, 90, 98 on the
  evidence above, the best states have stabilisers of order 6, 2, 2, 2 and sit
  6, 8, 16, 16 below the ceiling, with `ρ` between 0.9395 and 0.9470. Orders
  102 and 114 carry large-stabiliser states again.
- Annealing is not a certificate at larger orders: at `q = 173` it reaches
  1085 and never finds the torus state at 1131. At orders 74 to 114 the
  quotient search covers every state with a nontrivial stabiliser, with the
  involution quotient annealed rather than enumerated.

## 6. What this does and does not do

- It does not decide convergence, does not bound any `m_n` from below, and
  does not show any signing family has a permanent gap.
- It retires one target and sharpens another. A permanent gap
  (`lim sup α_n < 1/2`) would follow from an `η > 0` and primes
  `q_k ≡ 1 (mod 4)` with `q_{k+1}/q_k → 1` and `ρ(C_{q_k+1}) ≤ 1 − η`.
  Theorem 4 shows such a set must avoid small `a/b`; §5 shows it must also
  avoid the sporadic large-stabiliser states; Theorem 2 shows arithmetic alone
  gives only `1 − ρ ≥ c/q`. Order 62 is the first proved instance of a Paley
  matrix off its ceiling.
- Equivalent form of the target (Paley, `q` prime). For `f = x|_{GF(q)}` with
  periodic autocorrelation `A_f(d) = Σ_j f(j) f(j+d)`, `s = Σ_j f(j)` and
  `M_N(f) = Σ_{k nonsquare} |f̂(k)|²`,

      xᵀ C x = 2 x_∞ s + Σ_{d≠0} χ(d) A_f(d),
      S − Q_C(x) = ( M_N(f) + ½ (s − x_∞ √q)² ) / √q        (exact).

  So `ρ ≤ 1 − η` says: every `±1` sequence of length `q` keeps about `η q²/2`
  of its Fourier mass on the nonsquares; equivalently no binary sequence has
  periodic autocorrelation close to `√q χ(d)`. The congruence
  `A_f(d) ≡ q (mod 4)` for binary sequences is the classical relative of
  Lemma 1.
- `defect` is an exact tool at orders far beyond enumeration whenever the
  window between a witness and the ceiling holds few admissible patterns.
  Order 74 has 26 patterns in its window (values 303 and 305, supports 14 to
  26, mostly loose `|t| = 1` vertices); the first four tried did not finish in
  a minute each, so `Φ(C74)` stays in `{297, 303, 305}`.

**Next unresolved target.** A proof, for some explicit ratio-dense set of
primes, that `M_N(f) + ½ (s ∓ √q)² ≥ η q²` for every `±1` sequence `f` of
length `q`. Order 62 (`(1 − ρ)/2 = 2.7%` of `q²`) is the base case on file.

**Update 2026-11-13.** No such set of primes exists with a fixed `η > 0`:
`ρ(C_{q+1}) → 1` along all primes `q ≡ 1 (mod 4)` (Theorem S of
`NOTE_2026-11-13_PALEY_PHASE_MODEL.md`, proved modulo review of one cited
equidistribution theorem). The order-62 value stands as a finite fact; the
target above is withdrawn, not proved.

## Sources

- K. Momihara, S. Suda, *Conference matrices with maximum excess and
  two-intersection sets*, arXiv:1611.01305 (Proposition 1.1, Theorem 1.2).
- E. Hecke, *Eine neue Art von Zetafunktionen und ihre Beziehungen zur
  Verteilung der Primzahlen II*, Math. Z. 6 (1920).
- `solution.md` Proposition 15.762; `evidence/NOTE_2026-09-03_CONFERENCE_CUBE_GAP.md`.

## Reproduction

See `scripts/parity_ceiling/README.md`. Receipts:
`evidence/parity_ceiling_20261111/{ceilings, torus, witnesses, small_maximisers, symmetric_states, stabilisers, c62_exclusion}.json`.
Tests: `python3 -m pytest tests/test_parity_ceiling_note.py -q`.
