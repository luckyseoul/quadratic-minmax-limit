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

## 5. Routes killed

1. **Floor on every profile stratum is false.** The restricted Phi floors
   include p=7 k=4 `5.68889` and p=11 k=4 `4.79742`, both below 6. Any profile
   argument must use mixture or constituent structure, not PSD of every
   restricted stratum.

2. **Exceptional floor on every PSL orbit is false.** At p=7 there is a
   complete PSL orbit of size 1,176 on which `Z_psi=0`; its orbit-restricted
   exceptional scalar is zero. `(QVAR)` is an ensemble-mixing statement, not
   a pointwise or orbitwise inequality.

## 6. Live targets

Leftover 1 is now reduced to two explicit inequalities:

1. prove the quartic variance `(QVAR)`;
2. prove the principal variance room
   `||delta||^2 <= n(n+10)^2/[6(n-6)^2]`.

The other two original leftovers (residual (ii) and Type-I multi-level) are
unchanged.

Reproducible artifacts: `src/e1_gmin_m4_prop15589.py`,
`tests/test_prop15589.py`, and the `evidence/maxplus_p11/stratum_*` scripts and
JSON records.
