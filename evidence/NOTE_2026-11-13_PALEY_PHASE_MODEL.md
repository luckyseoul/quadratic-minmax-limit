# Paley phase model: Jacobi-sum diagonalisation, exact square-wave formula, and `rho -> 1` along prime orders

**Status:** draft session note (2026-11-13). The original MO limit is OPEN;
nothing here is entered in `src/original_mo_status.py`.

- Sections 1-2 are proved below (elementary, self-contained) and every
  identity is re-checked in exact integer arithmetic.
- Section 3 is a proof that uses one external theorem (Rojas-León,
  arXiv:2207.12439v3, Theorem 4). Its hypotheses are verified in the text, but
  the application was made from the statement of that theorem, not from its
  proof, and has not been independently reviewed. Treat Theorem S as
  **proved modulo that review**.
- Sections 4-6 are certificates and numerics. Section 7 is an exact value at
  order 46 unrelated to Paley.

**Implication targeted** (AGENTS rule). Theorem E(2) of `solution.md` asks
whether `rho(C_k)` converges along Paley orders. This note gives
`rho(C_{q+1}) -> 1` along all primes `q = 1 (mod 4)`. Prior art searched:
the two-block decomposition (`NOTE_2026-09-11`), the projector form
(Prop 15.18), the torus and quartic witnesses (`NOTE_2026-11-11`,
`NOTE_2026-11-12`). None diagonalises the form by Jacobi sums or treats a
frequency other than the quartic one.

Notation: `q = 1 (mod 4)` prime, `chi` the quadratic character, `g` a
primitive root, `n' = q-1`, `M = n'/2` (even), `e(t) = exp(2 pi i t)`,
`psi(g^d) = e(d/n')`, `C` the Paley conference matrix on `{inf} u GF(q)`
(`C_inf,u = 1`, `C_uv = chi(u-v)`), `S = (q+1) sqrt(q) / 2` the spectral
ceiling, `rho(C) = Phi(C)/S`. Sections 1-2 hold verbatim for prime powers
`q = 1 (mod 4)`; the code is for primes.

## 1. The phase model (proved)

Write a sign vector on `GF(q)*` as `z_s = x(g^s)`, `s in Z/n'`.

**Lemma 1 (twisted circulant).** `chi(g^s - g^t) = (-1)^s h(t-s)` with
`h(d) = chi(1 - g^d)`, `h(0) = 0`, and `h(-d) = (-1)^d h(d)`.

*Proof.* `chi(g^s - g^t) = chi(g^s) chi(1 - g^(t-s))` and `chi(g^s) = (-1)^s`.
Also `chi(1 - g^-d) = chi(g^-d) chi(-1) chi(1 - g^d)`. ∎

**Lemma 2 (Jacobi sums).** Put
`J_k = sum_d h(d) e(kd/n') = sum_{x != 0} chi(1-x) psi^k(x) = J(chi, psi^k)`.
Then `J_{-k} = conj(J_k)`, `J_{M-k} = J_k`, `J_0 = J_M = -1`, and
`|J_k|^2 = q` for `k` not in `{0, M}`.

*Proof.* `h` is real. `J_{M-k} = sum_d h(d) (-1)^d e(-kd/n') = sum_d h(-d)
e(k(-d)/n') = J_k` by Lemma 1. `J_0 = sum_{y != 1} chi(y) = -1`. The modulus is
the classical value of a Jacobi sum of two non-trivial characters with
non-trivial product; `psi^k` is trivial only for `k = 0` and `chi psi^k` only
for `k = M`. ∎

**Proposition 3 (phase model).** With `zhat(k) = (1/n') sum_s z_s e(-ks/n')`,

    sum_{s,t} chi(g^s - g^t) z_s z_t = n' sum_k J_k zhat(k) zhat(M-k),

and for the conference form

    2 Q_C(x) = n' sum_k J_k zhat(k) zhat(M-k)
               + 2 x_inf x_0 + 2 x_inf sum_s z_s + 2 x_0 sum_s (-1)^s z_s.

*Proof.* By Lemma 1 the left side is `sum_s w_s (h*z)_s` with
`w_s = (-1)^s z_s`, `(h*z)_s = sum_d h(d) z_{s+d} = sum_k zhat(k) J_k e(ks/n')`,
and `what(-k) = zhat(M-k)`. The border terms are `C_inf,0 = 1`,
`C_inf,g^s = 1`, `C_0,g^s = chi(-g^s) = (-1)^s`. ∎

So the form couples frequency `k` only to `M-k`, through one unit complex
number `J_k / sqrt(q) = e(phi_k / 2 pi)`. The self-paired frequencies are
`k = +-M/2` (the quartic characters); there `|Im J| = a`, `|Re J| = b` for
`q = a^2 + b^2`, `a` even, which is Theorem Q of `NOTE_2026-11-12`. Every
other frequency was untouched until now.

## 2. Square-wave states: an exact formula (proved)

Let `k in Z/n'` with `M' = M / gcd(k, M)` even, `k' = k / gcd(k, M)` (then
`k'` is odd), `delta in Z/M'`, and `sq(v) = +1` if `v mod M'` lies in
`[0, M'/2)`, `-1` otherwise. Define the state

    x(g^(2u)) = sq(k' u),   x(g^(2u+1)) = sq(k' u - delta),   x(0) = x(inf) = 1.

**Theorem R.** With `phi_j = arg J_j`,

    Q_C(x) = 1 + (1/2)(q-1) sqrt(q) * rho*(k, delta),

    rho*(k, delta) = (4 / M'^2) * sum_{m odd, 0 < m < M'}
                     sin(phi_{mk}) sin(pi m (2 delta + k') / M') / sin^2(pi m / M').

The weights `4 / (M'^2 sin^2(pi m / M'))` are positive and sum to 1.

*Proof.* Put `E(u) = z_{2u}`, `O(u) = z_{2u+1}` on `Z/M`; then
`zhat(j) = (1/2)[Ehat(j) + e(-j/n') Ohat(j)]` with transforms on `Z/M`. Both
are `M'`-periodic, and the substitution `v = k'u` gives `Ehat(mk) = c_m`,
`Ohat(mk) = c_m e(-m delta / M')` for odd `m mod M'`, zero at all other
frequencies, where `c_m = (4/M') / (1 - e(-m/M'))` is the transform of `sq`.
The two frequencies of `Z/n'` above `mk mod M` are `mk` and `mk + M`, and

    zhat(mk + eps M) = (1/2) c_m [1 + (-1)^eps e(-theta_m / 2 pi)],
    theta_m = pi m (2 delta + k') / M'.

The partner `M - (mk + eps M)` is `(-m)k + (1-eps) M`, so
`zhat(j) zhat(M-j) = -(i/2) (-1)^eps |c_m|^2 sin(theta_m)`. Summing over `eps`
with `J_{mk+M} = conj(J_{mk})` (Lemma 2) gives

    sum_j J_j zhat(j) zhat(M-j) = sum_{m odd mod M'} |c_m|^2 sin(theta_m) Im J_{mk}.

`mk` is never `0` or `M` modulo `n'` for `0 < m < M'`, so
`Im J_{mk} = sqrt(q) sin(phi_{mk})`, and `|c_m|^2 = 4 / (M'^2 sin^2(pi m/M'))`.
Both `sum_s z_s` and `sum_s (-1)^s z_s` vanish because `sq` is balanced, so the
border contributes `x_inf x_0 = 1`. Parseval for `sq` gives the weight sum. ∎

**Corollary R1.** `Phi(C_{q+1}) >= 1 + (1/2)(q-1) sqrt(q) rho*(k, delta)` for
every such `(k, delta)`.

**Corollary R2 (harmonic bound).** Let `eps_m = (-1)^((m-1)/2)`, `r >= 1`,
`M' > 4r`, `sigma = +-1`, and suppose

    sigma eps_m sin(phi_{mk}) >= 1 - eta     for m = 1, 3, ..., 2r-1.

Then for a suitable `delta`,

    rho*(k, delta) >= 1 - eta - 4/(pi^2 r) - 2 pi^2 r^2 / M'^2.

*Proof.* Take `2 delta + k' = sigma tau (mod 2M')` with `tau = M'/2` if that
is odd and `tau = M'/2 + 1` otherwise. Then
`sin(pi m (2 delta + k')/M') = sigma eps_m cos(pi m / M')^[tau != M'/2]`. The
summand of Theorem R is invariant under `m -> M' - m`, so the terms with
`m <= 2r-1` or `m >= M'-2r+1` are each at least
`(1 - eta) cos((2r-1) pi / M')` and carry total weight
`W >= (8/pi^2) sum_{m <= 2r-1 odd} 1/m^2 > 1 - 2/(pi^2 r)`; all other terms are
at least `-1` with total weight `1 - W`. ∎

In words: a square wave of frequency `k` is nearly a `+sqrt(q)` eigenvector as
soon as the Jacobi phases at `k, 3k, 5k, ...` sit near `-+pi/2` with
alternating signs. The quartic states of the two previous notes are the
degenerate frequency where all odd multiples coincide.

## 3. Saturation along prime orders (proved modulo review of one citation)

**Theorem S.** `Phi(C_{q+1}) / S -> 1` as `q -> infinity` through primes
`q = 1 (mod 4)`. Equivalently `Phi(C_n) / n^(3/2) -> 1/2` along prime Paley
orders.

*Proof.* Fix `r` and `eta > 0`. It suffices to find, for all large `q`, an odd
`k` with `gcd(k, M) <= sqrt(M)` satisfying the hypothesis of Corollary R2;
then `M' >= sqrt(M)`, and Corollaries R1-R2 give
`liminf Phi/S >= 1 - eta - 4/(pi^2 r)` for every `r, eta`, while `Phi <= S`.

*Equidistribution.* For `j != 0, M`: `J_j = G(chi) G(psi^j) / G(chi psi^j)` and
`G(chi) = sqrt(q)`, so `e^(i phi_j) = gamma(psi^j) / gamma(chi psi^j)` with
`gamma = G / sqrt(q)` on the unit circle. Odd `k` are `k = 1 + 2j`, i.e.
`psi^k = psi Lambda^2` with `Lambda = psi^j` arbitrary (two-to-one). Hence

    gamma(psi^(mk)) = gamma(psi^m Lambda^(2m)),
    gamma(chi psi^(mk)) = gamma(chi psi^m Lambda^(2m)),    m = 1, 3, ..., 2r-1,

are `2r` Gauss sums of monomials in one variable character `Lambda`, with
exponents `2m <= 4r` and fixed twists `psi^m`, `chi psi^m`. Theorem 4 of
arXiv:2207.12439v3 (finite fields of varying characteristic, bounded
exponents, arbitrary twists) gives joint equidistribution of these `2r`
numbers in `(S^1)^(2r)` provided the vectors

    v_{m,eta} = sum_{xi : xi^(2m) = eta} xi,   eta in {psi^m, chi psi^m},

are linearly independent in the rational vector space on characters. Fix
`xi_0` with `xi_0^2 = psi` and `kappa` of order 4. For `q > 4r` the solution
sets are `xi_0 mu_{2m}` and `xi_0 kappa mu_{2m}` (`mu_d` = characters of order
dividing `d`). Every element of `kappa mu_{2m}` has order divisible by 4 and
no element of `mu_{2m'}` has, so the two families have disjoint supports; and
within a family the vector for `m` contains a character of exact order `2m`,
which lies in no `mu_{2m''}` with `m'' < m`. So the hypothesis holds. The map
`(s, s') -> s/s'` on each pair sends Haar measure to Haar measure, hence the
points `(phi_k, phi_{3k}, ..., phi_{(2r-1)k})`, `k` odd, equidistribute on the
`r`-torus.

*Conclusion.* The box `{sigma eps_m sin(theta_m) > 1 - eta for all m}` is open
with positive Haar measure and null boundary, so it contains at least `c M`
odd `k` for large `q`, `c = c(r, eta) > 0`. At most
`sum_{d | M, d > sqrt(M)} n'/d <= 2 sqrt(M) tau(M) = o(M)` odd `k` have
`gcd(k, M) > sqrt(M)`. ∎

**Review items for Theorem S.** (i) The exact definition of the vector space
and of the solution sets in the cited hypothesis; (ii) that Theorem 4 there
allows the twists to depend on `q` with no bound on their order, as its
statement indicates; (iii) the two-to-one parametrisation of odd characters.
The numerics of Section 5 test the conclusion of the equidistribution step
directly and agree with it.

**Remark (small characteristic; not used).** For `q = p^e` the identity
`G(Lambda^p) = G(Lambda)` forces `phi_{pk} = phi_k`. The target pattern
`sigma eps_m` is multiplicative in `m`, so this is consistent when
`p = 1 (mod 4)` and conflicts at harmonic `p` when `p = 3 (mod 4)`; square
waves then lose about `2/(p^2+1)`. No statement is made for prime powers.

## 4. What this settles and what it does not

- **E(2) on prime orders, with `rho_* = 1`.** The 11-11 note had
  `limsup rho = 1`; the "uniform bound off resonance" named there and on 11-12
  as the next target does not exist. The deficits seen below `q = 1000`
  (`rho` between 0.93 and 1) are a finite-size effect: among `M` frequencies
  one can align about `log q / log(pi / theta)` harmonics to tolerance
  `theta`, so the deficit decays roughly like `1 / log q`.
- **No Paley matrix of prime order can certify a constant below 1/2.** Any
  proof of `liminf alpha_n < 1/2` has to use other signings.
- **Nothing here bounds `m_n` from below.** E(1) is untouched.
- **Heuristic, not a theorem.** The mechanism needs only a flat spectrum on
  paired frequency blocks and about `n` frequencies with generic phases. A
  signing that is a twisted circulant over a cyclic group of order comparable
  to `n` with a flat spectrum should saturate the same way. Candidates for a
  constant below 1/2 should therefore have no cyclic symmetry of order
  comparable to `n`.

## 5. Certificates

`scripts/phase_model/verify.py` -> `evidence/phase_model_20261113/identities.json`
(`q` = 13, 29, 61, 101, 181, 197, 1009, 2521, 10009): Lemma 1 entrywise, Lemma 2,
`|Im J_{M/2}| = a`, Proposition 3 against exact integer `Q` (and against the
dense matrix for `q <= 2521`), and Theorem R on 40 random `(k, delta)` per `q`
including even `k` and `gcd(k, M) > 1`; largest deviation of the closed formula
from the exact integer is `9.1e-7` (at `q = 10009`, `Q ~ 5e5`).

`scripts/phase_model/scan.py` -> `evidence/phase_model_20261113/scan.json`.
For each `q`: the best square-wave state (all odd `k` coprime to `M` and all
shifts for `q <= 40000`; preselected by 16 harmonics above), its exact integer
`Q` by Kronecker substitution, and the exact `Q` after a deterministic batch
1-opt ascent.

      q      | k      | delta  | Q square  | rho sq | Q polished | rho pol
    ---------|--------|--------|-----------|--------|------------|--------
    101      | 17     | 4      | 437       | 0.8526 | 437        | 0.8526
    229      | 29     | 14     | 1477      | 0.8487 | 1549       | 0.8901
    509      | 231    | 202    | 5005      | 0.8700 | 5243       | 0.9113
    1009     | 65     | 344    | 14585     | 0.9092 | 14753      | 0.9197
    2017     | 1025   | 739    | 41033     | 0.9055 | 41585      | 0.9177
    5009     | 3561   | 98     | 168473    | 0.9503 | 168535     | 0.9506
    10009    | 2495   | 3      | 472657    | 0.9439 | 473857     | 0.9463
    20021    | 1223   | 1865   | 1343749   | 0.9486 | 1347821    | 0.9515
    50021    | 19659  | 8943   | 5358669   | 0.9580 | 5362851    | 0.9587
    100049   | 22241  | 26421  | 15214481  | 0.9615 | 15276329   | 0.9654
    200009   | 175993 | 86901  | 43219001  | 0.9663 | 43236225   | 0.9667
    500009   | 97749  | 138631 | 170547233 | 0.9647 | 170547313  | 0.9647
    1000033  | 617699 | 316170 | 482674929 | 0.9653 | 482681409  | 0.9653

Each `Q` is a rigorous lower bound for `Phi(C_{q+1})`. These are not the
maxima: at `q = 101` the exact value is 501 (`rho = 0.9775`), and below 1000
the symmetric-state windows of 11-12 are higher. The table shows the square
waves alone climbing as Theorem S predicts.

## 6. Equidistribution check

Same file, per `q`: the largest of 22 Weyl sums
`|mean_k e(i(a phi_k + b phi_3k + c phi_5k))|`, `|a|,|b|,|c| <= 3`, times
`sqrt(q)`, stays between 3.3 and 7.0 for `q` from 101 to 1000033; and the number
of odd `k` whose first `r` phases lie within `theta` of the square-wave target
matches the Haar prediction `2 M (theta/pi)^r`:

      q      | r=1, 0.3      | r=2, 0.3    | r=3, 0.3  | r=4, 0.5
    ---------|---------------|-------------|-----------|----------
    10009    | 940 / 956     | 80 / 91     | 12 / 9    | 8 / 6
    100049   | 9488 / 9554   | 912 / 912   | 52 / 87   | 64 / 64
    500009   | 47252 / 47747 | 4500 / 4560 | 436 / 435 | 376 / 321
    1000033  | 95112 / 95496 | 8812 / 9119 | 952 / 871 | 600 / 642

## 7. Exact value at order 46 (Mathon), a non-Paley conference matrix

45 is not a prime power, so there is no Paley matrix of order 46. The parity
ceiling of `NOTE_2026-11-11` at `q = 45` with the row identity leaves one
admissible top level: `Q = 147` with fields `{3: 10, 7: 33, 11: 3}`. Hence
`Phi(C) <= 147` for **every** symmetric conference matrix of order 46.

`evidence/phase_model_20261113/mathon46.json` holds six graphs
`SRG(45,22,10,11)` (graph6, public data, source recorded in the file), the
bordered conference matrices (`C^2 = 45 I` checked), and for each a sign vector
with `Q = 147` and exactly that field pattern. The six matrices have the same
descendant 4-clique profile, consistent with one switching class. So

    Phi(C_46, Mathon) = 147,     m_46 <= 147,     alpha_46 <= 0.4712.

This is a non-Paley, non-geometric conference matrix sitting exactly on the
parity ceiling. No order-46 bound was on file before.

## 8. Next unresolved target

1. Independent review of the citation step in Theorem S (items (i)-(iii)).
2. With Paley prime orders closed at `rho_* = 1`, the only finite test on file
   that can still move the constant is the order-50 one: a conference matrix
   of order 50 with no Boolean `+-7` eigenvector has `Phi <= 167 < 169` by
   Prop 15.762. Every construction available offline (nets `L_4(7)`,
   Pasechnik, two-block group-developed) is regularizable. Spence's list
   `conf2graph.50` of the known regular two-graphs on 50 vertices decides it;
   the file is served as binary and could not be read from this session.
3. An effective rate in Theorem S needs a discrepancy bound in growing
   dimension; not attempted.

## Receipts

    python3 scripts/phase_model/verify.py
    python3 scripts/phase_model/scan.py 101 229 509 1009 2017 5009 10009 20021 \
        50021 100049 200009 500009 1000033        # about 5 minutes
    python3 -m pytest tests/test_paley_phase_model_note.py -q
