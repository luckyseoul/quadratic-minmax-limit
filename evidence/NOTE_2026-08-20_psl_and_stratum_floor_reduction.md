# Leftover-1 reduction: PSL constituents, quartic exception, and stratum diagnostics

Date: 2026-08-20. Branch: `codex/leftover-moment-attack`.

No settlement flag is flipped. This note records one general representation
theorem, two exact reductions, and two failed proof routes.

## 1. Exact PSL decomposition

Put `q=p^2`, `n=q+1`, `d=n/2`, and let `W_e=V_+` be the degree-`d` even
Weil constituent of `PSL(2,q)`. Direct substitution in the standard ordinary
character table gives

```
Sym^2(W_e) = 1 + St + W_e + sum_{alpha in A_e} rho(alpha),
|A_e| = (q-9)/8,
```

where the `rho(alpha)` are distinct principal-series irreducibles of degree
`q+1`. The diagonal image is the projective-line module `1+St`; hence

```
Z = W_e + sum_{alpha in A_e} rho(alpha).
```

This is multiplicity-free and has the required dimensions

```
dim Z   = d + ((q-9)/8)n = n(n-6)/8,
dim Z^U = 1 + 2(q-9)/8   = (q-5)/4.
```

GAP's independent character table and power-map calculation reproduces the
decomposition at `q=25,49,121`: one exceptional plus respectively 2, 5, and
14 principal constituents, all once.

Character-table source used for the symbolic class-family calculation:
Jeffrey Adams, *Character tables for GL(2), SL(2), PGL(2), and PSL(2) over a
finite field*, Section 6.4 (`https://math.umd.edu/~jda/characters/characters.pdf`).

Because `Phi` commutes with `PSL(2,q)`, it is scalar on each constituent.
Every eigenvalue not supported solely on `W_e` therefore has multiplicity at
least `n`. The only possible sub-`n` eigenvalue is the single scalar
`lambda_exc=Phi|W_e`, of multiplicity `d=n/2`.

Thus the old missing assertion `mult(lambda_min)>=n` is replaced exactly by
the one scalar check `lambda_exc>=6`.

## 2. The exceptional scalar is one quartic variance

The `U`-fixed line in `W_e` is the quartic multiplicative-character mode
`psi^2=chi`. If `D` is the minus set of a Max+ vector, put

```
N(a)    = |D intersect (D-a)|,
Z_psi   = sum_{a != 0} psi(a) N(a).
```

Combining the Fourier identities of Props 15.279 and 15.473 gives

```
lambda_exc = 32 E|Z_psi|^2 / (q(q-1)).
```

Consequently the entire exceptional issue is

```
E|Z_psi|^2 >= 3q(q-1)/16.                 (QVAR)
```

This is substantially narrower than proving the floor on every
multiplicative-character mode.

Exact checks:

| p | E|Z_psi|^2 | threshold | lambda_exc |
|---:|---:|---:|---:|
| 5 | `3300/13` | `225/2` | `176/13` |
| 7 | `317520/409` | `441` | `4320/409` |
| 11 | `3931.461697314` | `2722.5` | `8.664378396284` |

The p=11 value independently reproduces the unique multiplicity-61 spectral
cluster.

## 3. Variance route after the decomposition

The spectral mean is `mu=8(n-2)/(n-6)`. If the minimum is on principal
constituents, multiplicity at least `n` gives the sufficient condition

```
Var(spec Phi) <= 32(n+10)^2/(n-6)^3,
||delta||^2   <= n(n+10)^2/[6(n-6)^2].
```

If the minimum were the exceptional `d`-block, the variance room is exactly
half as large. It is preferable to prove `(QVAR)` separately and retain the
larger principal room.

## 4. Profile-stratum diagnostics

The p=11 enumeration is in profile order `k=1,3,4,5,6`. Every stratum clears
`(QVAR)` on its own:

| stratum | weight | E|Z_psi|^2 | induced scalar |
|---:|---:|---:|---:|
| k=1 | 0.0000740 | 108900 | 240 |
| k=3 | 0.0006461 | 21780 | 48 |
| k=4 | 0.0015506 | 9438 | 20.8 |
| k=5 | 0.0348879 | 8765.7778 | 19.3185 |
| k=6 | 0.9628414 | 3727.3813 | 8.21461 |

The same stratumwise inequality holds in the exact p=5 and p=7 censuses;
the weakest finite case is p=7, k=4, with induced scalar `7.82222`.

This suggests a profilewise proof of `(QVAR)`, but it is not yet a theorem:
the general distributions in the `k>=4` coefficient spaces remain unnamed.

For the four-point residual, p=11 is not primarily a cancellation miracle:
the dominant k=6 stratum already has `||delta_k6||^2=2.15503`, and the full
mixture has `||delta||^2=2.36856`. By contrast p=5 and p=7 require substantial
cross-stratum cancellation.

## 5. General low-profile theorem

The quartic target is now proved on every `k=1` and `k=3` profile stratum,
for every prime.  Only `k>=4` remains in `(QVAR)`.

When `p=3 mod 4`, the quartic character is trivial on `F_p^*`.  Grouping
ordered pairs of the minus set by their projective `F_p`-direction gives,
pointwise,

```
Z_psi(y) = sum_L psi(g_L) a_L(y),
a_L(y)   = (1/4) sum_s (sigma_L(s)-y_inf)^2 >= 0.
```

The functions `(sigma_L-y_inf)/2` are mean-zero ridge functions on distinct
directions, hence orthogonal.  The profile reconstruction and the fixed number
of plus/minus entries then give the second pointwise identity

```
sum_L a_L(y) = p(p^2-1)/4.                  (ENERGY)
```

For `k=1`, all energy is on one direction.  For `k=3`, each affine profile is
a centered permutation of `F_p`, so all three energies equal one third of
`(ENERGY)`.  The direction character has `(p+1)/4` signs of each kind and every
triple has `(p-1)p^2` lifts.  Consequently, with `S=p(p^2-1)/4` and
`m=(p+1)/2`,

```
E_k1 |Z_psi|^2 = S^2,
E_k3 |Z_psi|^2 = S^2 (m-3)/(3(m-1)).
```

Both exceed `3p^2(p^2-1)/16` for every `p=3 mod 4`, `p>=7`.  They reproduce
the exact values `7056,784` at p=7 and `108900,21780` at p=11.

When `p=1 mod 4`, `eta=psi|F_p^*` is the Legendre character.  A `k=3`
affine profile has three nonzero Fourier lines.  On each line the lifted affine
permutation contributes, up to a Gaussian-unit phase, the common magnitude

```
p S_p,
S_p = sum_{a!=0} eta(a)/|1-exp(2 pi i a/p)|^2
    = p^2 L(2,eta)/(2 pi^2).
```

The Euler product gives
`L(2,eta) >= product_l (1+l^-2)^-1 = zeta(4)/zeta(2)=pi^2/15`, so
`S_p>=p^2/30`.  A sum of three Gaussian units has modulus at least one.
Therefore `|Z_psi|^2>=p^6/900`, which clears `(QVAR)` for `p>=13`; p=5 is
the exact base case `|Z_psi|^2=180>225/2`.  On `k=1`, fixed-size subset
sampling in the signed Paley graph gives the exact average

```
E_k1 |Z_psi|^2 = p^3 (p-1)^2 (p+1)/(8(p-2)),
```

which also clears `(QVAR)`.

Independent diagnostics: direct pair sums agree exactly with the profile
formula on all 11,452 p=7 vectors and on sampled p=11 vectors, including the
pointwise total `(ENERGY)`.  The complete p=13 k=4 array has
`E|Z_psi|^2=8788` and induced scalar `9.9047619`; 869 valid p=13 high-activity
MILP samples give `7514.63` and scalar `8.46957`.  These p=13 values are data,
not a general `k>=4` theorem.

The full p=11 directional-energy covariance gives one more useful negative.
On the complete ensemble the quartic sign vector is the top covariance mode,
with eigenvalue `655.2436162` and `6*655.2436162=3931.4616973`.  But it is the
smallest nonzero mode on k=4 (`1573` versus `2510.75`) and is below the top on
the dominant k=6 stratum (`621.2302` versus `625.6471`).  Therefore no symmetry
or stratum-invariance statement can force the quartic mode to be top on each
piece; a full-mixture top-mode argument would need an additional cross-stratum
theorem.  Reproducible record:
`evidence/maxplus_p11/directional_energy_covariance_p11.{py,json}`.

Reproducible diagnostic: `evidence/quartic_profile_attack.py`.

## 6. Odd-coset shell and spherical benchmark

There is a useful lattice reformulation, with one important trap.  Set

```
L = ker_Z(C-pI) = V_+ intersect Z^n.
```

Fix any `y0` in Max+.  The vectors of `L` having every coordinate odd are
exactly `y0+2L`: if `x` and `y0` are odd then `(x-y0)/2` is integral and is
still in `V_+`.  Every vector in this coset has squared norm at least `n`,
with equality exactly when all its coordinates are `+/-1`.  Consequently the
full antipodal Max+ family is precisely the first shell of the odd coset.

It is **not** the first shell of `L`.  If `A` is an affine `F_p`-line with
square direction, then

```
r_A = 1_{{infinity} union A},       C r_A = p r_A,
||r_A||^2 = p+1 < p^2+1=n.
```

Indeed the character sum along `A` is `p-1` on `A` and `-1` off `A`.
Thus ordinary minimum-shell design theorems address a shorter shell and do
not by themselves control Max+.

The coset formulation nevertheless gives a clean benchmark.  Extend
`K_psi(a,b)=psi(a-b)` by zero in the infinity row and column and put
`A_psi=P K_psi P/4` on `V_+`.  Since `K_psi 1=0`,
`Z_psi(y)=y^T A_psi y`.  Additive and multiplicative character
orthogonality give

```
tr(A_psi)=0,          ||A_psi||_HS^2=q(q-1)/32.
```

For the uniform radius-`sqrt(n)` sphere in `d=n/2` dimensions, the standard
quadratic-form moment is therefore

```
V_sph = n^2 [2||A_psi||_HS^2+|tr A_psi|^2] / [d(d+2)]
      = q(q-1)(q+1) / [4(q+5)].
```

This already exceeds the QVAR threshold, by the exact amount

```
V_sph - 3q(q-1)/16 = q(q-1)(q-11) / [16(q+5)] > 0       (p>=5).
```

The exact remaining statement can now be phrased as a lower bound on the
degree-four harmonic excess of the first shell of `y0+2L`: it may be negative,
but not below the negative of the displayed gap.  Proving that harmonic theta
coefficient is nonnegative would be a convenient stronger sufficient result;
nonnegativity has not been proved.

An independent exact audit with PARI/GP, `python-flint`, and `cypari2` found
rank `(p^2+1)/2` and shorter lattice minimum `p+1` at `p=5,7,11,13`.  The
general shorter-vector construction above, rather than that audit, is what
kills the ordinary-minimum-shell shortcut.

## 7. Routes killed

1. **Floor on every profile stratum is false.** The restricted Phi floors
   include p=7 k=4 `5.68889` and p=11 k=4 `4.79742`, both below 6. Any profile
   argument must use mixture or constituent structure, not PSD of every
   restricted stratum.

2. **Exceptional floor on every PSL orbit is false.** At p=7 there is a
   complete PSL orbit of size 1,176 on which `Z_psi=0`; its orbit-restricted
   exceptional scalar is zero. `(QVAR)` is an ensemble-mixing statement, not
   a pointwise or orbitwise inequality.

3. **Ordinary lattice minimum-shell design is the wrong shell.**  The explicit
   norm-`p+1` vectors above are strictly shorter than Max+'s norm `p^2+1` for
   every `p>=5`.  The replacement target is the degree-four harmonic
   coefficient of the first *odd-coset* shell, with the exact allowed negative
   budget displayed in Section 6.

4. **The quartic direction is not the top covariance mode stratumwise.**  The
   exact p=11 k=4 and k=6 covariance spectra above are counterexamples.  A
   proof may still bound the full quartic mode directly, but cannot obtain it
   by asserting that every invariant profile component puts that mode on top.

## 8. Live targets

Leftover 1 is now reduced to two explicit inequalities:

1. prove the quartic variance `(QVAR)` on the union of profile strata `k>=4`;
2. prove the principal variance room
   `||delta||^2 <= n(n+10)^2/[6(n-6)^2]`.

The other two original leftovers (residual (ii) and Type-I multi-level) are
unchanged.

Reproducible artifacts: `src/e1_gmin_m4_prop15589.py`,
`tests/test_prop15589.py`, and the `evidence/maxplus_p11/stratum_*` scripts and
JSON records.
