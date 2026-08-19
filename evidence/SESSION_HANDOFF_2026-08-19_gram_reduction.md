# Session handoff 2026-08-19 — Max+ Gram reduction and leftover audit

Branch: `prop15586-maxplus-gram-reduction`. Ships 15.586 and 15.587.
**No leftover flag is flipped.** E(1) remains open.

## 1. The reduction (15.586, general p, proved)

`P = (I+C/p)/2`, `Max+ = {y in {+-1}^n : Cy = py}`, `N = |Max+|`,
`Z = {B sym : CB = pB, diag B = 0}`, `Phi(B) = E_y[(y^T B y)^2]`.

Every `y` in Max+ has `Py = y` and `y_i^2 = 1`, so for `pi_i = P e_i`

    <y y^T, pi_i pi_i^T> = (y^T P e_i)^2 = ((Py)_i)^2 = y_i^2 = 1

for every `i` and every `y`. Z-perp inside `{B : PBP = B}` is `span{pi_i pi_i^T}`,
so `proj_Zperp(y y^T) = R` is the SAME element for all y, with
`R = 2 sum_i pi_i pi_i^T` and `||R||^2 = 4 tr(P^2) = 4 tr P = 2n`. Hence

    spec(Phi) = nonzero spec(Ghat / N),   Ghat_ab = <y_a,y_b>^2 - 2n

with `#nonzero = dim Z`. `Ghat = (YY^T)o(YY^T) - 2n J` is an INTEGER matrix and a
genuine Gram matrix (PSD by construction) — unlike `G_{u,disj}`, which package
Caveat 1 forbids treating as a Gram because it has negative eigenvalues. That
obstruction was an artifact of the chosen object, not of the problem.

Closed forms (all general p): `dim Z = n(n-6)/8`, `tr Phi = n(n-2)`,
`tr K = -4n` for `K = 8I - Phi`, and `dim span{1, y_i y_j} = n(n-6)/8 + 1`.
`E[y y^T] = I + C/p = 2P`, so the Wick value is exactly `8||B||^2` and the floor
is equivalent to `lambda_max(K) <= 2` — the package's `<delta,psi> <= 2` as an
explicit operator with known trace.

## 2. Exact spectra (finite, p=5,7)

    p=5: 80/13 (mult 26), 144/13 (26), 176/13 (13)
    p=7: 3072/409 (50), 3360/409 (100), 3648/409 (50), 4032/409 (50), 4320/409 (25)

Bottom multiplicity is exactly n, top exactly n/2, both primes.
`lambda_max(K) = 48/n` EXACTLY at p=5 (= 24/13) and 200/409 < 48/50 at p=7.
So the candidate `lambda_* = 8(n-6)/n` is attained at p=5 only: **the floor binds
at the smallest prime**, where p=5 is already a finite check, and by p=7 there is
4x slack. A crude bound closes it for p >= 7 — no sharp identity needed.

## 3. Type I mu identified (15.587)

`mu = max_{|kappa(S)|=1} |E_{y in Max+}[y_i y_j y_k y_l]|` with
`kappa(S) = C_ij C_kl + C_ik C_jl + C_il C_jk` in {-3,-1,1,3}. Recomputed from an
independent exhaustive Max+ enumeration: 3/65 at p=5 (11700 of 14950 four-sets),
109/2863 at p=7 (176400 of 230300) — reproducing `census_gmin_kappa1` exactly.
Denominators are `p*D = N/4`, the same normalisation as the Phi spectrum:
**Type I and the floor are moments of one tensor.**

Bound targets, verified against `L_abs_gmin`/`T_abs` for p = 5..43:

    L = (p-2)/(2p^2)        T = (p-2)/(p(2p-1))        T > L for all p >= 5

`T > L` is exactly why `|mu| <= |T|` cannot close Type I. Margin to the real
target is stable: `mu/L = 0.7692` at p=5, `0.7462` at p=7 — again ~25% headroom,
so a crude bound suffices here too.

Also: `max|m4|` over ALL four-sets is 21/65 (p=5) and 327/2863 (p=7), both < 1.

## 4. Leftover audit — E(1) is three stubs

Each remaining leftover was traced to its blocking predicate. In every case the
surrounding machinery is fully proved and the gate is a hardcoded False:

| Leftover | Blocking predicate | Surroundings |
|---|---|---|
| 1 floor | `phi_F_ge_6_proved_general` | reduction A-D above, general p |
| 2 residual (ii) | `multilevel_ND_k_ge_4p_proved` | theorems A-L all proved |
| 3 Type I | `type_I_aut_e_3AB_positive_general` | 39/39 lemmas True at p=5..23 |

## 5. Negative results — do not re-run

- **Max+ inner-product classes are NOT an association scheme.** Max within-class
  variance of `A_i A_j` is 113 at p=5 and 5369 at p=7. So the spectrum of Ghat is
  not determined by the inner-product distribution alone.
- **The 15.237 C 0-1 pair-span classification survives an exhaustive test.**
  Criterion used: a 0-1 function lies in the pair-span iff `q_B(y) = y^T B y`
  takes at most two values on Max+. Exhaustively over all 5668650 support-3 sets
  at p=5: 3575 two-valued = 2600 triangles + 975 extras. Every extra spans 6
  vertices (three DISJOINT edges whose degree-6 monomial is identically +1) — a
  genuinely different configuration, but every observed mass (1/10, 1/5, 3/10,
  2/5 and complements) lies in `classified_01_pairspan_masses`. No unclassified
  mass found. At p=7 no non-triangle extras in 598510 sampled support-3 sets, so
  the family looks p=5-specific. Support-2 is settled outright: two edges are
  two-valued only if some `|m4| = 1`, and the global max is far below 1.

## 6. The wall, stated once

All three stubs are statistics of Max+, and Max+ is enumerable only for p <= 7:
at p=11 the nullity is 61, so the sweep is 2^61 ~ 2.3e18, and the structured
families give `n_1d = 2772` and `n_{k=3} = 24200` while `n_full` is exactly the
unclassified family. **Closing E(1) is one problem, not three: get Max+ moments
at general p.** The reduction narrows the floor's share of that to 2-point data
(the Gram Ghat); Type I still needs genuine 4-point moments on |kappa|=1 sets.

## 7. Infrastructure note

CUDA 13 dropped sm_70, so GPU JIT is broken for the V100 on this host:
`numba.cuda` fails with `libnvvm: -arch=compute_70 is an unsupported option`, and
any CuPy path needing NVRTC (axis reductions, `cp.unique`, `.sum()` on a bool
array, RawKernel) fails compiling `cuda_fp4.hpp`. Precompiled cuBLAS, cuSOLVER,
elementwise ops and boolean masking still work. Workaround used throughout:
replace axis reductions with a matvec against a ones-vector.
