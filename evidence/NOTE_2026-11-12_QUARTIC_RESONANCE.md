# Quartic square-wave states: the Paley-core lift never undercuts by more than p − 1, and Φ(C198) = 1373

**Status:** draft session note (2026-11-12), independent review pending.
Theorems Q and L are proved below (elementary; one Jacobi sum). §4 is an
exact certificate. §6 is a table of rigorous two-sided windows and claims
nothing about values inside a window. The original MO limit is OPEN. Nothing
here bounds any `m_n` from below and nothing is entered in
`src/original_mo_status.py`.

**Implications established.**

1. *The lift.* For every prime `m ≡ 1 (mod 4)`, `m = a² + b²`, both coherent
   lifts `K_±(Q_m) = [[Q, Q ± I], [Q ± I, −Q]]` of the Paley core satisfy
   `Φ(K_±) ≥ (m − 1)(a + b) + 1`. At `n = 2m = p² + 1` this is
   `Φ(C) − (p − 1)`. So along the Paley-core lift family the undercut of the
   crystalline value is **at most `p − 1` for every `p`**. This is the
   inequality left open in D9 of
   `NOTE_2026-09-19_SUPERLEVEL_DISCRETE_DERIVATIVE.md` ("would make the gap
   `p − 1` a theorem along this family"), in the direction that matters for
   `γ_p`: this family can never certify more than
   `γ_p ≥ (p − 1)/n^{3/2}`. Equality holds at `p = 3, 5, 7` by the exact values
   on file (13, 61, 169).
2. *Paley.* `Φ(C_{q+1}) ≥ ½(q − 1)a + 1` for every prime `q = a² + b²`,
   `a` even. With the parity ceiling this gives the new exact value
   `Φ(C198) = 1373`, and it explains 5, 33, 109, 501 at orders 6, 18, 38, 102.
3. Rigorous windows `L ≤ Φ(C_{q+1}) ≤ U` for all 80 primes `q ≡ 1 (mod 4)`
   below 1000.

**Reused.** Lemma 1, Theorem 2, Theorem 4 of
`NOTE_2026-11-11_CONFERENCE_PARITY_CEILING.md`; the exact lift values and the
identity `Q_K(x,y) = Q_B(x) − Q_B(y) + xᵀBy + x·y` of
`NOTE_2026-09-19_PLUS_I_COHERENT_LIFT.md`; `NOTE_2026-11-11_K50_MINUS_I.md`.
Dedup: `quartic`, `cyclotomic`, `Jacobi sum` occur in the repository only in
the August dual-shell and coprimality notes, not for Boolean states of Paley
matrices or lifts.

## 1. One lemma

`m` is a prime, `m ≡ 1 (mod 4)`; `χ` the quadratic character; `Q` the Paley
core, `Q_xy = χ(x − y)` on `GF(m)`; `ψ` a character of `GF(m)*` of order 4, so
`ψ² = χ`. Put

    J = J(χ, ψ) = Σ_z χ(1 − z) ψ(z)  ∈  Z[i],        |J|² = m.

`Re J` is a sum of `(m − 3)/2` signs, hence odd; `Im J` is a sum of
`(m − 1)/2` signs, hence even. So with `m = a² + b²`, `a` even, `b` odd:
`|Im J| = a`, `|Re J| = b`.

For `α ∈ C` let `g_α(y) = Re(α ψ(y))` for `y ≠ 0` and `g_α(0) = 0`.

**Lemma.** (i) `Q g_α = g_{ᾱ J̄}`. (ii) `⟨g_α, g_β⟩ = ½(m − 1) Re(α β̄)`.
(iii) `g_α` is orthogonal to `1` and to `χ`, and `(Q g_α)(0) = 0`.

*Proof.* For `x ≠ 0` substitute `y = xz`:
`Σ_{y≠0} χ(x − y) ψ(y) = χ(x) ψ(x) J = ψ³(x) J`, and
`Re(α J ψ³(x)) = Re(ᾱ J̄ ψ(x))`. At `x = 0` the sum is `Σ ψ³ = 0`. (ii) and
(iii) are orthogonality of `ψ, ψ², ψ³`. ∎

For `θ ∈ {π/4, 3π/4, 5π/4, 7π/4}` the **square wave**
`w_θ = g_{√2 e^{iθ}}` takes values `±1` on `GF(m)*`: it is `+1` on two
consecutive quartic classes and `−1` on the other two. The Lemma gives

    w_θᵀ Q w_θ  = (m − 1) Re(e^{2iθ} J),
    w_θᵀ Q w_θ' = (m − 1) Re(e^{i(θ+θ')} J),
    w_θ · w_θ'  = (m − 1) cos(θ − θ').

## 2. Theorem Q (Paley conference matrices)

**Theorem Q.** For every prime `q ≡ 1 (mod 4)`, `q = a² + b²` with `a` even,

    Φ(C_{q+1})  ≥  ½ (q − 1) a + 1.

*Proof.* Take `x = (x_∞; u)` with `u = u₀ δ₀ + w_θ`. Then
`xᵀCx = 2 x_∞ Σ_j u(j) + uᵀQu = 2 x_∞ u₀ + (q − 1) Re(e^{2iθ} J)`, using
(iii). `e^{2iθ} = ±i`, so `Re(e^{2iθ} J) = ∓ Im J`; choose `θ` with value `+a`
and `x_∞ u₀ = 1`. ∎

The fields of this state are `1` at `∞` and `0` and, on the four quartic
classes, either `(a+b+2, a+b−2, a−b, a−b)` or `(a−b+2, a−b−2, a+b, a+b)`.

Together with Theorem 4 of the parity-ceiling note,

    Φ(C_{q+1})  ≥  max( ½ (q + 1) b ,  ½ (q − 1) a + 1 ),

so `ρ(C_{q+1}) ≥ max(cos ω, (1 − 2/(q+1)) sin ω)` with `tan ω = a/b`: the
torus state is strong when the Gaussian angle `ω` is near `0`, the square wave
when it is near `90°`. The larger of the two is at least `1/√2 − O(1/q)`, and
at least `0.745` of `S` for every prime below 1000.

## 3. Theorem L (coherent lifts of the Paley core)

For `(u, v) ∈ {±1}^m × {±1}^m`,
`Q_{K_±}(u, v) = ½ uᵀQu − ½ vᵀQv + uᵀQv ± u·v`.

**Theorem L.** For every prime `m ≡ 1 (mod 4)`, `m = a² + b²`,

    Φ(K_±(Q_m))  ≥  (m − 1)(a + b) + 1.

*Proof.* Take `u = u₀ δ₀ + w_θ`, `v = v₀ δ₀ + w_θ'`. By the Lemma the terms
containing `δ₀` vanish except `± u₀ v₀`, and

    Q_{K_±} = (m − 1) [ ½ Re(e^{2iθ}J) − ½ Re(e^{2iθ'}J) + Re(e^{i(θ+θ')}J) ± cos(θ − θ') ] ± u₀ v₀.

With `θ' = θ + π/2` the bracket is `Re((1 + i) e^{2iθ} J)`, and with
`θ' = θ − π/2` it is `Re((1 − i) e^{2iθ} J)`. Since `e^{2iθ} = ±i`, the four
values `±(Re J + Im J)`, `±(Re J − Im J)` all occur; the largest is
`|Re J| + |Im J| = a + b`. Choose `u₀ v₀` so that the last term is `+1`. ∎

**Corollary (crystalline orders).** If `2m − 1 = p²` then
`(a + b)² + (a − b)² = 2m = p² + 1`, and `2m` has one representation as a sum
of two squares, so `a + b = p`. Hence

    Φ(K_±(Q_m))  ≥  (m − 1) p + 1  =  m p − (p − 1)  =  ½ n √(n − 1) − (p − 1),     n = 2m = p² + 1.

Here `m p` is the value of the crystalline Paley matrix of order `p² + 1` when
it exists. The primes `m` of this form have `a, b` consecutive, i.e. Gaussian
angle `45°`: exactly where the two Paley families of §2 are weakest, the lift
is resonant.

The same proof works over a prime-power field. Over `GF(25)` (`p = 7`) the
square-wave pair evaluates to 169 for both signs, the exact value on file
(checked in `tests/test_quartic_resonance_note.py`).

      p     m    order    m p    lift state    status
      3     5      10      15        13        = Φ, exact (m_10)
      5    13      26      65        61        = Φ, exact
      7    25      50     175       169        = Φ, exact (both signs)
      9    41      82     369       361        lower bound
     11    61     122     671       661        lower bound
     15   113     226    1695      1681        lower bound
     19   181     362    3439      3421        lower bound
     25   313     626    7825      7801        lower bound
     29   421     842   12209     12181        lower bound
     ...  (35, 39, 45, 49, 51, 59 likewise; `quartic.json`)

Checked in exact arithmetic for all 59 primes `m < 700` and both signs
(`scripts/parity_ceiling/quartic.py`).

Two remarks for anyone running lifts.

- Full-cube annealing at order 122 reached only 645 in 6,000 restarts and
  never found the 661 state; the 10-orbit quotient of the dilation subgroup of
  order 15 finds it at once. At order 82 annealing found 361 in 18 of 6,000
  runs. Annealing is not evidence of an undercut at these orders.
- Whether `Φ(K_±(Q_m)) = m p − (p − 1)` exactly for `p ≥ 9` is open. The fields
  of a lift state are all congruent mod 4 (same argument as Lemma 1 of the
  parity-ceiling note, because `K² = I₂ ⊗ ((2m+1)I − 2J + 2Q)` is `≡ 0 (mod 4)`
  off the diagonal), but `Σ h²` is not constant, so the ceiling argument does
  not transfer as is.

## 4. Exact value at order 198

`q = 197 = 14² + 1`. Theorem Q gives `Φ(C198) ≥ 98·14 + 1 = 1373`. Ceiling
(`q n = 39006`): in the class `h ≡ 1 (mod 4)` write `h = 13 + 4t`; then
`2T₂ = 693 − 13 T₁`, so `T₁` is odd and `T₂ ≥ T₁` forces `T₁ ≤ 46`;
`T₁ = 45` needs `T₂ = 54` (parity), so `T₁ ≤ 43` and
`Q = 1287 + 2T₁ ≤ 1373`. In the class `h ≡ 3`, `h = 15 + 4t` gives
`2T₂ = −693 − 15T₁`, `T₁ ≤ −55`; `T₁ = −55` needs `T₂ = 66` (parity), so
`Q = 1485 + 2T₁ ≤ 1371`. Hence **`Φ(C198) = 1373`**, and every symmetric
conference matrix of order 198 has `Φ ≤ 1373`.

For `q = a² + 1` with `a ≡ 2 (mod 4)` the same computation gives ceiling
`= ½a³ + 1` exactly when `a ≤ 14`, i.e. `q = 5, 37, 101, 197`; from `a = 18`
on a window of 8 opens (`q = 677`: `[8789, 8797]`).

## 5. Resonance, in one picture

With `ω` the Gaussian angle of the prime (`tan ω = a/b`):

    Paley, torus state        ρ ≥ cos ω                 (strong near 0°)
    Paley, square wave        ρ ≥ sin ω (1 − 2/(q+1))   (strong near 90°)
    lift,  square-wave pair   Φ / (m √(2m−1)) ≥ cos(ω − 45°)(1 − O(1/m))   (strong near 45°)

Every explicit family on file is pinned to its flat value by a state that is
constant on the orbits of a large subgroup, and the pinning is governed by
`ω`. The lift does not escape: the orders `p² + 1` force `ω = 45°`.

## 6. Windows for all primes below 1000

`windows.py` takes `L` to be the best of: the torus state, the square wave,
and every state fixed up to sign by a subgroup with at most 26 orbits (cyclic
subgroups of both tori, `A₄`, `S₄`, `A₅`), each quotient maximised exactly and
the maximiser re-evaluated on the full matrix. `U` is the parity ceiling.
Receipt: `windows.json` (80 primes, 74 stored witnesses).

- Exact (`L = U`): `q = 5, 13, 17, 29, 37, 41, 53, 101, 197`; with the order-62
  certificate of the previous note, ten exact Paley values.
- Windows of width at most 10: `173: [1131, 1133]`, `229: [1725, 1731]`,
  `257: [2049, 2055]`, `293: [2499, 2509]`, `401: [4001, 4009]`,
  `577: [6913, 6921]`, `677: [8789, 8797]`, and `73, 109, 113`.
- Certified `L/S ≥ 0.97` at 28 of the 80 primes, `≥ 0.95` at 47. The best
  certified state comes from the square wave (20 primes), the torus (17),
  other cyclic quotients (28), `A₄/S₄/A₅` (15).
- `L` is a lower bound from small quotients only. It is weak (`< 0.90`) at
  `q = 157, 277, 317, 613, 653, 773`; that is absence of a small-quotient
  state, not a bound on `Φ`. At `q = 137` annealing has a state at 777 while
  the small-quotient `L` is 759.

## 7. What this does and does not do

- It closes the Paley-core lift as a source of a permanent gap: its undercut is
  at most `p − 1` at every crystalline order where it is defined. It says
  nothing about lifts of other blocks.
- It does not decide convergence and does not bound any `m_n` from below.
- It does not change the unresolved target of the previous note (a uniform
  bound off resonance on a ratio-dense set of primes). It adds the second
  explicit resonance (`ω` near `90°`) and the table shows that small-orbit
  states above `0.97` are common, so that set has to avoid more than small
  `a/b`.
- New exact data: `Φ(C198) = 1373`; hence `m_198 ≤ 1373`
  (`1373/198^{3/2} = 0.4928`, a resonant order, of no use as an undercut).

## Reproduction

    python3 scripts/parity_ceiling/quartic.py 1000
    gcc -O2 -o /tmp/qopt scripts/parity_ceiling/qopt.c -lm
    QOPT=/tmp/qopt python3 scripts/parity_ceiling/windows.py 5 1000     # about 25 minutes on two cores
    python3 -m pytest tests/test_quartic_resonance_note.py -q

Receipts: `evidence/quartic_resonance_20261112/{quartic, windows, lift_search}.json`.
