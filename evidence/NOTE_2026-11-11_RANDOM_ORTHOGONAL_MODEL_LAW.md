# The flat-spectrum minima track the random orthogonal model

2026-11-11. Claude Fable 5.2 (`claude-fable-5-2`), single session, asked to
attack the limit by a route not already in the repository.

**Status.** The original MO limit is still OPEN. Nothing below was entered
in `src/original_mo_status.py`. This note contains one proved lemma (§1),
one proved vector inequality (§6), exact finite evaluations (§2, §4), and an
empirical law with predictions that can fail (§3, §5). It contains no
convergence proof and no cross-order bound for the actual minima.

Dedup check before writing: `random orthogonal`, `Marinari`, `R-transform`,
`free probability`, `sqrt(15)`, `0.4841`, `n^2(n-1)/4` have no hits in
`*.md`/`*.py`; `Haar` appears only in an August counterexample built from
random equal-diagonal projectors (`solution.md`) and in two August GPU logs;
none of them compares ground states. The cut-code,
Gale--Berlekamp, fourth-phase, two-block, tensorization and nuclear-budget
material that I rederived on the way is already in the September notes and
is not repeated here.

Notation: `Q_A(x) = sum_{i<j} A_ij x_i x_j`, `Phi(A) = max_x |Q_A(x)|`,
`C_{q+1}` the Paley conference matrix of `F_q`, and

    S_n = (1/2) n sqrt(n-1)   (n = 2 mod 4, conference scale)
    S_n = (1/2) n^(3/2)       (n = 0 mod 4, symmetric Hadamard of trace 0)

the spectral ceiling for a signing whose spectrum is flat.

## 1. Lemma (first moment for a Haar flat spectrum)

Let `n` be even, `D = diag(I_{n/2}, -I_{n/2})`, `O` Haar on `O(n)`,
`J = O D O^T`. Then for every `delta > 0`

    P( max_{x in {+-1}^n} |x^T J x| >= (sqrt(15)/4 + delta) n ) -> 0

exponentially fast. Since `tr J = 0`, `x^T J x` is the zero-diagonal form
of `J`. In the repository's normalisation (`A = sqrt(n) J`, all eigenvalues
`+-sqrt(n)`, off-diagonal entries of variance about one) this reads

    Phi(sqrt(n) J) <= (sqrt(15)/8 + delta) n^(3/2),   sqrt(15)/8 = 0.48412...

**Proof.** Let `P_-` be the spectral projector of `J` on eigenvalue `-1`.
For fixed `x`, `x^T J x = n - 2|P_- x|^2 = n(1 - 2B)` with
`B = |P_- x|^2 / n`. `O^T x / sqrt(n)` is uniform on the sphere, so `B` is
`Beta(a, a)` with `a = n/4`. For `0 < eps < 1/2` the density is increasing
on `[0, eps]`, hence

    P(B <= eps) <= eps (eps(1-eps))^(a-1) / B(a,a).

The duplication formula gives
`B(a,a) = sqrt(pi) Gamma(a) / (2^(2a-1) Gamma(a+1/2)) >= sqrt(pi) / (2^(2a-1) sqrt(a))`,
so `P(B <= eps) <= sqrt(a) (4 eps (1-eps))^a`. A union bound over the
`2^(n-1)` antipodal classes gives

    P( max_x x^T J x >= (1-2 eps) n ) <= (sqrt(n)/4) (64 eps (1-eps))^(n/4).

The base is below one exactly when `eps(1-eps) < 1/64`, i.e. when
`1 - 2 eps > sqrt(15)/4`. `-J` has the same law. QED.

This is the entropy-crisis point of the random orthogonal model (ROM) of
Marinari, Parisi and Ritort. With the R-transform `R(t) = (sqrt(1+4t^2)-1)/(2t)`
of the symmetric Bernoulli law, the annealed entropy
`log 2 - (1/4) log((1 + sqrt(1+4 beta^2))/2)` vanishes at `beta^2 = 240`
(`T = 0.0645`, their `T_S = 0.065`), where the energy is `sqrt(15)/8`.
Their statement that deterministic orthogonal couplings share the ROM
thermodynamics "but for special properties that the deterministic model can
have for special N" is the phenomenon in §3--§4.

## 2. Exact finite-size values of the Haar model

`scripts/rom_law/rom.py` + `gs2.c`: exhaustive Gray-code enumeration of
`max_x x^T J x` and `max_x -x^T J x` for independent Haar draws. Column
"mean" is the mean of `max / n` over draws and both signs.

    n   draws  mean    sd        n   draws  mean    sd
    4   400    0.8209  0.146     18  300    0.9449  0.0155
    6   400    0.8745  0.085     20  200    0.9472  0.0129
    8   400    0.9093  0.050     22  200    0.9498  0.0112
    10  400    0.9210  0.039     24  100    0.9526  0.0103
    12  400    0.9284  0.029     26  60     0.9554  0.0086
    14  400    0.9375  0.022     28  40     0.9536  0.0074
    16  300    0.9417  0.019     30  24     0.9551  0.0061

Annealing estimates (lower estimates; the same annealer reproduces the
exhaustive `n = 30` mean to four digits): `0.9563` (34), `0.9579` (38),
`0.9585` (42), `0.9569` (50). In units of `S_n` the draw-to-draw standard
deviation is about `0.5` to `0.6` at every `n` from 14 to 30.

## 3. The law

For every even order from 6 to 38 at which a flat-spectrum signing exists
and a minimum or record is on file, the best known value of `Phi` is the
Haar mean times `S_n`, rounded down to the admissible parity, or one parity
step below that.

    n   S_n      Haar*S_n  admissible  best  gap    source of "best"
    6   6.71     5.87      5           5     0.87   m_6 exact
    8   11.31    10.29     10          10    0.29   m_8 exact
    10  15.00    13.82     13          13    0.82   m_10 exact
    12  20.78    19.30     18          18    1.30   m_12 exact
    14  25.24    23.66     23          21    2.66   m_14 exact = Phi(C14)
    16  32.00    30.14     30          30    0.14   record, m_16 in [27,30]
    18  37.11    35.06     35          33    2.06   record = Phi(C18)
    20  44.72    42.36     42          40    2.36   record
    24  58.79    56.00     56          56    0.00   record
    26  65.00    62.10     61          61    1.10   record
    28  74.08    70.65     70          70    0.65   record
    30  80.78    77.15     77          75    2.15   record = Phi(C30)
    38  115.57   110.71    109         109   1.71   record = Phi(C38)

Thirteen rows, every gap in `[0, 2.7]`, nine exact matches with the
admissible integer and four one step lower. The gap does not grow with `n`.
`n = 4` is the one small-order exception (`m_4 = 4 = S_4`, Haar `3.3`).
The Haar numbers for 24 and 28 were computed before I read the repository's
records 56 and 70.

Read in ratio form: `Phi(C_{q+1})/S` at the prime orders is `0.832`,
`0.889`, `0.928`, `0.943` (q = 13, 17, 29, 37) against Haar means `0.9375`,
`0.9449`, `0.9551`, `0.9579`; the difference is the constant two or so
units, not a different limit.

## 4. New exact evaluations in this session

All by `gs2.c` over the full projective cube; `C6 -> 5`, `C14 -> 21`,
`C18 -> 33` reproduce the shipped values.

- Flat Hadamard lift `K(C, C-I) = [[C, C-I],[C-I, -C]]` (the Seidel part of
  the trace-zero symmetric Hadamard matrix `[[C+I, C-I],[C-I, -C-I]]`):
  order 12 from C6: `18`; order 28 from C14: `70`; **order 36 from C18:
  `max Q = 108`, `min Q = -108`**, i.e. `Phi = S_36` exactly (2^35 states,
  185 s). The `+I` variants give `18` and `70` at orders 12 and 28.
- Principal restrictions (delete the first d points, infinity included):
  C14: `20, 20, 19, 17` (d = 1..4); C18: `32, 32, 29, 27`;
  C30: `72, 70, 69, 65, 64, 58` (d = 1..6); **C38 minus two points: `108`**
  (hence C38 minus one point is `108` by parity and monotonicity).

So the order-36 Hadamard lift is crystalline: it has Boolean eigenvectors on
both sides, like C10 (`15 = S_10`) and C26 (`65 = S_26`). At 10 and 26 the
minimum or record (`13`, `61`) sits at the Haar level and is attained by
other signings. At 36 nobody has searched; both algebraic candidates above
give `108`. A greedy pool descent with an annealing oracle
(not shipped) reproduced the floor `70` at order 28 and did
not move the order-36 lift off `108`; single-edge moves cannot leave a
crystalline maximiser set.

## 5. Predictions that can fail

Using the Haar means of §2 and the band "admissible integer below
Haar*S_n, or one step lower":

1. `Phi(C42)` (q = 41) is `127` or `125`. `S = 134.47`, Haar*S = `128.9`.
2. Some signing of order 36 has `Phi <= 104`; expected `102`.
   `S = 108`, Haar*S about `103.4` (mean interpolated between 34 and 38).
   Current bound: `108`.
3. Some signing of order 32 has `Phi <= 86`. `S = 90.51`, Haar*S about
   `86.5` (mean interpolated between 30 and 34).
4. Some signing of order 50 has `Phi <= 167`. `S = 175`, Haar*S >= `167.5`.
   On file: `169` (plus-I lift), `175` (C50).
5. `Phi(C54)` (q = 53) is `187` or `185`. `S = 196.56`, Haar*S about
   `188.5` (mean extrapolated from 42 and 50).

`Phi(C42) >= 129`, or a certified `m_36 >= 106`, would count against the
law. Item 1 is one enumeration of 2^41 states.

## 6. A proved by-product: the vector trade-off is an exact hyperbola

Let `V^+(A) = max (1/2)<A,X>` and `V^-(A) = max -(1/2)<A,X>` over
`X >= 0` with unit diagonal. For every complete signing of order `n >= 2`

    V^+(A) V^-(A) >= n^2 (n-1) / 4,

with equality for every conference matrix and for the all-negative signing.

**Proof.** By SDP duality `V^+ = (1/2) min{sum d_i : D - A >= 0}` and
`V^- = (1/2) min{sum d'_i : D' + A >= 0}` with `D, D'` diagonal. The 2 by 2
minors give `d_i d_j >= 1`, so all `d_i, d'_i > 0`. For any real diagonal
`Delta = diag(delta)`, since `A_ii = 0` and `A_ij^2 = 1`,

    0 <= tr[(D - A) Delta (D' + A) Delta]
       = sum_i (1 + d_i d'_i) delta_i^2 - (sum_i delta_i)^2.

Taking `delta_i = 1/(1 + d_i d'_i)` gives `sum_i 1/(1 + d_i d'_i) <= 1`.
Put `t_i = 1/(1 + d_i d'_i)` and `phi(t) = sqrt(1/t - 1)`, which is
decreasing and convex on `(0, 3/4]`. If all `t_i <= 3/4`, Jensen gives
`sum phi(t_i) >= n phi(1/n) = n sqrt(n-1)`. Otherwise exactly one
`t_i > 3/4`, the others are below `1/4` with mean below `1/(4(n-1))`, and
Jensen on those gives `sum >= (n-1) sqrt(4n-5) >= n sqrt(n-1)` for
`n >= 3`; `n = 2` is AM--GM. Cauchy--Schwarz finishes:
`(sum d_i)(sum d'_i) >= (sum sqrt(d_i d'_i))^2 >= n^2 (n-1)`. QED.

Checked with SCS on 60 random signings of orders 4 to 14 (minimum ratio
`1.0000`, maximum harmonic sum `1.0000`, both attained in the all-positive
class at small orders), C10, C14, and the all-negative
signing at orders 8 and 9 (ratio `1.00000` each).

The Boolean analogue is false at every tested order: `m_10 = 13` gives
`max Q * max(-Q) <= 169 < 225`. So the vector problem has no room below the
hyperbola and the Boolean problem does, which is the gap §1--§3 quantify.
This flips nothing.

## 7. What this does and does not do

It does not prove convergence, and it does not bound `m_N` by `m_n`.

It replaces the object of study. Every route on file tries to carry a
minimiser from one order to another. The table in §3 says the minima at flat
orders are already described, to within one parity step, by an ensemble
whose limit is a number: `e_H = lim E max_x x^T J x / n <= sqrt(15)/4`.
On that reading

    lim alpha_n = e_H / 2,   0.4776 <= e_H / 2 <= 0.48412,

the limit exists because both sides of the game are pinned to the same
constant, and `1/2` is not the value. The lower end is the exhaustive
`n = 30` mean; it is an observation, not a bound.

Two statements would turn the reading into a proof. Neither is proved here.

- **(U)** A ratio-dense family of signings with
  `Phi <= (e_H/2 + o(1)) n^(3/2)`. Lemma 1 is this statement for Haar
  conjugates, which are not sign matrices. For prime-order Paley it is the
  counting bound
  `#{x : |P_- x|^2 <= eps n} = 0` for `eps < (1 - sqrt(15)/4)/2 - o(1)`,
  equivalently: no binary sequence of prime length `q` has periodic
  autocorrelations `C_x(d)` with
  `sum_d chi(d) C_x(d) > (sqrt(15)/4 + o(1)) q^(3/2)`; equivalently the
  annealed bound on the Paley Ising partition function up to
  `beta = sqrt(240)`. (U) alone gives `limsup alpha_n <= 0.48412 < 1/2`.
- **(L)** Every signing has `Phi >= (e_H/2 - o(1)) n^(3/2)`. The proved
  floor is `B_infty`, about `0.3259`. I have no route to (L), and the
  constant two-unit deficit in §3 leaves open that sign matrices beat the
  Haar model by a vanishing or by a positive margin.

Do not record the law as a theorem. Do not cite §3 as a bound on any `m_n`.

## Sources

- E. Marinari, G. Parisi, F. Ritort, Replica field theory for deterministic
  models (II): a non-random spin glass with glassy behaviour,
  [arXiv:cond-mat/9406074](https://arxiv.org/abs/cond-mat/9406074).
- E. Marinari, G. Parisi, F. Ritort, The fully frustrated hypercubic model
  is glassy and aging at large D,
  [arXiv:cond-mat/9410089](https://arxiv.org/abs/cond-mat/9410089).
- G. Parisi, M. Potters, Mean-field equations for spin models with
  orthogonal interaction matrices,
  [arXiv:cond-mat/9503009](https://arxiv.org/abs/cond-mat/9503009).
- M. Degli Esposti, C. Giardinà, S. Graffi, Energy landscape statistics of
  the random orthogonal model,
  [arXiv:cond-mat/0207681](https://arxiv.org/abs/cond-mat/0207681).

## Repro

    cd scripts/rom_law && gcc -O3 -march=native -funroll-loops -o gs2 gs2.c
    python3 paley.py                 # C6, C14, C18 -> 5, 21, 33
    python3 rom.py 26 60             # one row of the Haar table
    python3 rom_sa.py                # annealing rows 30..50
    python3 lift36.py -1             # order-36 Hadamard lift, 2^35 states
    python3 restrict.py 29 1 2 3 4 5 6

Receipts: `evidence/rom_law_20261111/` (`rom_results.jsonl`,
`rom_sa_results.json`, `comparison_table.json`, `exact_evaluations.json`).
Host: 2 cores, no GPU. Checks: `tests/test_rom_law_note.py`.

The GPU dispatch suggestion and its stop rules are in the 2026-11-11 entry
of `HANDOFF.md`.
