# Relative-gauge composition is exactly a cross-block diamond after gauge absorption

**Status.** This note is a bridge to the current public upstream framework at
[`Robby955/mo-413935-research`](https://github.com/Robby955/mo-413935-research),
main commit `61ac2268b6f9234c3ea268e7a7072ac72141f36b` (commit date
2026-08-23). It proves an exact equivalence and a two-ray settling
specialization. It does **not** prove either ray or the MathOverflow limit.

## 1. Gauge absorption

Let `A,B` be symmetric zero-diagonal signings of orders `n,k`, let `C` be
an `n`-by-`k` sign matrix, and let

```text
g = ([alpha],[beta],tau) in P_n x P_k x {+1,-1}.
```

The upstream relative-gauge block is

```text
Y_g = [[D_alpha A D_alpha, C],
       [C^T, tau D_beta B D_beta]].
```

For a rectangular sign matrix `D`, define the two-block diamond

```text
Delta_tau(A,B;D)
  = max_(x,y) (|Q_A(x)+tau Q_B(y)| + |x^T D y|).
```

Then exactly

```text
M(Y_g) = Delta_tau(A,B; D_alpha C D_beta).             (1)
```

Indeed, conjugating `Y_g` by `D_alpha direct-sum D_beta` fixes its Boolean
maximum and changes it to

```text
[[A, D_alpha C D_beta],
 [(D_alpha C D_beta)^T, tau B]].
```

For fixed `x,y`, replacing `x` by `-x` leaves both internal quadratic
energies fixed and reverses the cross energy. Therefore

```text
max(|U+V|,|U-V|)=|U|+|V|
```

gives (1). The choices of representatives of `[alpha]` and `[beta]` only
negate the cross matrix, which does not change the diamond.

Consequently, for a fixed seed `C`,

```text
min_g M(Y_g)
 = min_(tau,[alpha],[beta])
   Delta_tau(A,B;D_alpha C D_beta),                    (2)
```

and, when the cross seed is itself free,

```text
min_C min_g M(Y_g)
 = min_(tau,D) Delta_tau(A,B;D).                       (3)
```

If the admissible family for `B` is closed under coefficient negation (in
particular, if all optimal order-`k` signings are allowed), `tau` can be
absorbed into `B` as well.

Thus a relative gauge does not enlarge the one-step class of block
signings. Its possible advantage is proof-theoretic: the balanced fibers,
their exact Fourier factorization, or an abundance theorem may locate many
good row/column switchings of one structured seed and thereby retain an
iterable state.

## 2. Exact relation to Propositions 6.5c and 6.8

At `n=k` and `B=A`, the `tau=-1` case of (1) is precisely Proposition
6.5c:

```text
M([[A,D],[D^T,-A]])
 = max_(x,y)(|Q_A(x)-Q_A(y)|+|x^T D y|).              (4)
```

The `tau=+1` case is the same-diagonal diamond. If `R` is a skew signing
and `D_0` is an arbitrary diagonal signing, take the complete cross matrix
`D=-R+D_0`. Grouping the two coordinates of each cloud in Proposition 6.5
by layers gives exactly this block matrix. Removing `D_0` changes every
energy by at most `n`, while

```text
Delta_+(A,A;-R)
 = max_(x,y)(|Q_A(x)+Q_A(y)|+|x^T R y|)
 = 2 K(A,R).                                           (5)
```

For orders `n` and `2n`, the `tau=+1` case of (1) is exactly the two-state
objective (6.39) in Proposition 6.8. Hence the upstream zero-temperature
existence problem specializes to the same internal/cross anticorrelation
already exposed by the local two-half and `1:2` diamonds. The Fourier fiber
labels are additional structure for proving that objective small, not a
different construction which has already done so.

## 3. A one-character abundance lemma

The shell factorization has a precise theorem-shaped sufficient target which
is stronger than merely showing that the occupancy is nonconstant. Let
`b:G -> Z_>=0`, let `K=|G|`, and use the upstream normalization

```text
mu = hat b(1) = K^(-1) sum_g b(g).
```

For any nontrivial `{+1,-1}`-valued character `chi`, put

```text
delta_chi = mu-|hat b(chi)| >= 0.
```

Then at least

```text
(K/2)(1-delta_chi)                                  (6)
```

fibers are empty whenever `delta_chi<1` (with the evident integer rounding).
In particular, `delta_chi<=1-epsilon` gives at least `epsilon K/2` empty
fibers, while `delta_chi=0` makes a whole character half empty.

To prove this, choose `s=sign(hat b(chi))`. The total occupancy on the
opposite character half is exactly

```text
sum_(chi(g)=-s) b(g) = (K/2)(mu-|hat b(chi)|)
                     = K delta_chi/2.                 (7)
```

Every occupied fiber contributes at least one to this integer sum, so at
most `K delta_chi/2` of the `K/2` fibers in that half are occupied. This
proves (6).

Combining (7) with the upstream factorization turns good-fiber abundance
into one explicit character-sum inequality:

```text
hat b_s(I,J,epsilon)
 = K^(-1) sum_(d_A+d_B+d_C<s)
   A_dA(I,epsilon) B_dB(J,epsilon) C_dC(I,J).
```

It is enough to make the absolute value of one nontrivial factored sum lie
within `1-epsilon` of the trivial coefficient. This is an exact sufficient
condition, not a proved estimate for signing shells. When `mu` is
exponential, additive accuracy one is exponentially fine in relative terms;
the lemma therefore sharpens the target but does not evade the upstream
precision wall unless an algebraic support or sign mechanism supplies the
near-saturation.

The one stored order-14 calibration does not exhibit that mechanism. A
direct read of its already-constructed Walsh transform gives raw trivial
coefficient `304908` and largest nontrivial absolute coefficient `19568`
(transform index 12). Hence

```text
mu = 304908/8192 = 76227/2048,
max_(chi!=1)|hat b(chi)| = 19568/8192 = 1223/512,
min_(chi!=1) delta_chi = 71335/2048 = 34.83154296875.
```

This is far above one, so no single-character abundance certificate sees
the unique empty order-14 fiber. Any successful use of the lemma needs new
asymptotic algebraic concentration, not extrapolation from that finite
example.

### Mod-two Walsh collapse

There is also a precise boundary on the proposed ``parity'' version of the
character-sum route.  It prevents a useless finite census of factored Walsh
coefficients modulo two, but does not rule out a direct shell-specific
pairing or an odd-prime congruence.

Let `G` be any nontrivial elementary abelian two-group and let
`b:G -> Z`.  Write `W_b(chi)=sum_g b(g)chi(g)` for the *unnormalized*
Walsh numerator (the normalized coefficient `hat b=W_b/|G|` is not defined
modulo two).  For every Walsh character `chi:G -> {+1,-1}`, reduction modulo
two gives

```text
W_b(chi) = sum_g b(g) chi(g) = sum_g b(g)  (mod 2).         (8)
```

Thus every factored coefficient in the labeled formula has the same
mod-two residue: the Walsh transform matrix itself reduces to the all-ones
matrix.  Fourier inversion supplies no recovery of individual fiber
residues in characteristic two, because its normalizing factor `1/|G|` is
not defined there.  In particular, transform-level parity alone cannot
distinguish an everywhere-odd occupancy from a pattern containing empty
fibers.

For example, the two integer occupancies `b(g)=1` for every `g`, and any
occupancy obtained by replacing two values `1,1` by `0,2`, have the same
total parity and hence identical mod-two Walsh-transform data by (8), while
only the latter has a vacancy.  This is an algebraic obstruction to a proof
that uses only the existing factored transforms reduced modulo two.  It does
not say that an actual signing shell realizes the example, nor does it
exclude a direct involution on its product triples, a congruence modulo an
odd prime, or an exact sign identity before reduction.

## 4. The only relative-gauge estimates needed for convergence

Let `m_n` be the optimal Boolean maximum and put `H(n)=m_n^(2/3)`. Let
`Omega_2,Omega_3 >= 0` have vanishing dyadic Dini tails:

```text
sum_(j>=0) sup_(u>=2^j N) Omega_r(u) -> 0
```

for `r=2,3`. It is enough to prove the following two fiber statements for
all sufficiently large `n`.

1. Choose optimal order-`n` blocks `A_n,B_n` and a cross seed `C_n`. With
   `L_2=M(A_n)+M(B_n)+B(C_n)` and

   ```text
   T_2=2^(3/2)m_n+n^(3/2)Omega_2(n),
   ```

   the upstream occupancy has an empty fiber at deficit `L_2-T_2`.

2. Choose optimal blocks `A_n,B_(2n)` and an `n`-by-`2n` seed `E_n`. With
   `L_3=M(A_n)+M(B_(2n))+B(E_n)` and

   ```text
   T_3=(m_n^(2/3)+m_(2n)^(2/3))^(3/2)
       +n^(3/2)Omega_3(n),
   ```

   the corresponding occupancy has an empty fiber at deficit `L_3-T_3`.

The upstream max-plus identity turns these assertions into

```text
m_(2n) <= T_2,
m_(3n) <= T_3.                                         (8)
```

Since `m_n=Theta(n^(3/2))`, the mean-value theorem turns (8) into the two
Dini-summable `H`-inequalities of Proposition 6.3. Therefore these two
statements imply convergence. No uniform all-pairs good-fiber theorem is
needed.

By (1), each empty-fiber conclusion is equivalently one good cross-block
diamond in the seed's row/column switching orbit. Thus the theorem-shaped
hinge is not the already-proved shell factorization itself. It is a uniform
sign or lower-tail theorem for the factored occupancy which reaches these
two critical thresholds with a Dini-summable error.

### 4.1 Positive-temperature certificates do not need integer-level accuracy

Added 2026-09-13. This combines the gauge construction above with the
standard pressure/norm sandwich in
`NOTE_2026-09-02_THERMODYNAMIC_INTERPOLATION_GATE.md`. It corrects the
overstatement that every analytic treatment of the connected layers must
resolve an exponentially large occupancy to absolute accuracy below one.
That precision requirement belongs to the occupancy certificate in Section
3. The following alternative has a different error budget. Neither
certificate's main inequality has been proved for the required signings.

Write `N=n+k`, and use the **same raw inverse temperature** `beta>0` for
all three blocks. Define normalized partition functions

```text
Z_A(beta) = E_x cosh(beta Q_A(x)),
Z_B(beta) = E_y cosh(beta Q_B(y)),
Z_C(beta) = E_(x,y) cosh(beta x^T C y),
Z_g(beta) = E_(x,y) cosh(beta Q_(Y_g)(x,y)).
```

Expectations here are uniform over full Boolean cubes. With the uniform
gauge `g=([alpha],[beta_g],tau)` from Section 1, set

```text
rho_g(beta) = Z_g(beta)/(Z_A(beta) Z_B(beta) Z_C(beta)).   (9)
```

Then `rho_g>0` and `E_g rho_g=1` exactly. To see this, average over the
internal gauges first. For each fixed `(x,y)`, the gauged internal spins
are independent uniform projective spins, so their quadratic energies have
the same laws as `Q_A` and `Q_B`. Averaging `tau` turns the right-block
exponential factor into `cosh(beta Q_B)`. In the remaining average,
`x -> -x` reverses the cross energy and cancels the term containing its
hyperbolic sine. Expanding the last cosine therefore leaves
`Z_A Z_B Z_C`. In general omitting `tau` does not give this product.
The letter `beta_g` denotes the right gauge and is unrelated to temperature.

There are `2^(N-1)` projective spin states. Including the two exponential
terms in `cosh` gives, for every gauge,

```text
exp(beta Phi(Y_g))/2^N <= Z_g(beta) <= exp(beta Phi(Y_g)).
```

Consequently

```text
Phi(Y_g) <= [log Z_A + log Z_B + log Z_C
              + log rho_g + N log 2]/beta.              (10)
```

This gives a precise one-sided error criterion. Let `T` be the desired raw
norm bound, and suppose a gauge and a real model value `ell_g` satisfy

```text
log rho_g <= ell_g + N epsilon,       epsilon >= 0,
log Z_A + log Z_B + log Z_C + ell_g <= beta T.           (11)
```

Only an upper error bound at the selected gauge is needed; there is no
requirement to reconstruct every gauge density. Equations (10)--(11) imply

```text
Phi(Y_g) <= T + N(log 2 + epsilon)/beta.                 (12)
```

For the two rays in Section 4, take `k=(r-1)n`, `r=2,3`, and set

```text
T_(2,n) = 2^(3/2) m_n,
T_(3,n) = (m_n^(2/3)+m_(2n)^(2/3))^(3/2),
c_n = log(n+1)^2,       beta_n = c_n/sqrt(r n).
```

If (11) holds for all sufficiently large `n` on both rays with their
stated optimal internal blocks and `0<=epsilon_(r,n)<=E_0`, where `E_0`
is independent of `n`, then (12) proves the Section 4 hypotheses with

```text
Omega_r(n) = r^(3/2)(log 2 + E_0)/log(n+1)^2.           (13)
```

These errors have the required vanishing dyadic Dini tails. For `L>=2`,
`a=log L`, `b=log 2`, and `K_r=r^(3/2)(log 2+E_0)`, monotonicity gives

```text
sum_(j>=0) sup_(u>=2^j L) Omega_r(u)
 <= K_r sum_(j>=0) 1/(a+jb)^2
 <= K_r (1/a^2 + 1/(ab)) -> 0.                         (14)
```

For the second inequality, each `j>=1` term is at most
`[1/(a+(j-1)b)-1/(a+jb)]/b`; summing telescopes. Thus a one-sided
`O(N)` error in the **whole log density** is admissible at this temperature
schedule. It produces an `O(1/log(n+1)^2)` normalized norm error, not an
integer-level occupancy estimate.

The unresolved content is (11), especially its second inequality and its
validity at the growing temperatures `c_n`. A fixed-temperature remainder
bound with uncontrolled dependence on `c` does not suffice. The previously
recorded `H4/H6/H8` collisions disprove exact recovery from those truncated
states; they do not disprove a bound on the whole remainder in (11).
The normalization `E_g rho_g=1` alone also does not supply the generally
negative log-density bound needed there.

The proof above is analytic. A new, single-fixture rational-polynomial
check in `relative_gauge_temperature_20260913/` corroborates the exact
normalization, the projective normalization factor, and the necessity of
the `tau` average. It does not verify (11), the Dini hypotheses for actual
signings, or convergence. No new source-signing census is involved.

## 5. What the current upstream commit actually proves

The upstream commit proves the balanced max-plus identity, exact shell
Fourier factorization, a Parseval improvement, a complete finite
integer-moment vacancy hierarchy, and the first mixed four-cycle Hamiltonian.
It also proves that generic low-degree or polynomial-precision vacancy tests
can miss an isolated empty fiber exponentially badly. Its order-14 example
has one empty fiber, but this is a finite calibration.

It does not prove either statement in Section 4. In particular:

- the advertised good-fiber abundance theorem is an active target;
- the exact character-sum sign needed to force a critical empty fiber is an
  active target;
- the mixed four-cycle Hamiltonian is only the first nonconstant term, and
  the upstream text explicitly leaves higher connected Eulerian layers
  uncontrolled on the mean-field diagonal;
- the all-pairs power-saving composition inequality is a sufficient
  hypothesis, not an established estimate.

Accordingly, nothing in current upstream `main` closes or strictly advances
the local residual diamonds (6.20) or (6.42)--(6.43). Its genuinely extra
content is an exact general labeling language and finite vacancy machinery.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. pytest -q \
  tests/test_upstream_relative_gauge_bridge.py
```
