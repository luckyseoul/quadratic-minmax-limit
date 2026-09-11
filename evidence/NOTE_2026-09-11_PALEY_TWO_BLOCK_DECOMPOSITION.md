# Paley conference two-block decomposition via the nonsquare multiplier

**Status:** draft session note (2026-09-11). The Lemma in §1 is proved below
(elementary, self-contained). §2–§4 are exact re-runnable certificates.
Independent review per repository rules pending. No all-orders or convergence
claim is made.

## 1. Decomposition (proved)

Let `q = 1 (mod 4)` be an odd prime power, `chi` the quadratic character of
`GF(q)` (`chi(0)=0`, `chi(-1)=+1`), and `C` the Paley conference of order
`n = q+1` on `{inf} u GF(q)`: `C_inf,inf = 0`, `C_inf,u = 1`, `C_uv = chi(u-v)`
for `u != v`. Let `S` be the nonzero squares, pick a nonsquare `nu`, and set

```
T  = {inf} u S,        T' = {0} u nu*S.
```

**Lemma.** For all `i, j in T`: `C_ij = -C_{pi(i) pi(j)}`, where `pi(inf) = 0`
and `pi(u) = nu*u` on `S`. Hence, with the blocks listed in these orders and
rows/columns reordered to `[T; T']`,

```
C = [[A, B], [B^T, -A]],     A = C[T,T],   B = C[T,T'].
```

*Proof.* (i) `i = inf`, `j = u in S`: `C_inf,u = 1`, and
`C_{0,nu*u} = chi(0 - nu*u) = chi(-1) chi(nu) chi(u) = (1)(-1)(1) = -1`.
(ii) `i, j in S`: `C_ij = chi(i-j)` and `C_{nu*i, nu*j} = chi(nu*(i-j)) =
chi(nu) chi(i-j) = -chi(i-j)`. ∎

No switching is required in these listings. (Under arbitrary listings of the
same partition the identity holds up to diagonal switching of one block.)

## 2. Certificates

- `scripts/paley_two_block_decomposition.py` — asserts the entrywise identity
  and the full two-block reordering for `q = 5, 9, 13, 17, 25, 29, 37`; prints
  block `Phi` by exact enumeration. All pass; exit 0.
- `tests/test_paley_two_block_decomposition.py` — pytest checks (identity,
  prime powers, recorded block minima, balanced partitions). All pass.

```sh
python3 scripts/paley_two_block_decomposition.py
python3 -m pytest tests/test_paley_two_block_decomposition.py -q
```

## 3. Consequence at q = 5, 13, 17: an optimal conference is a two-block
signing with a half-order minimizer block

Block `Phi` (exact) versus recorded exact half-order minima
(`m_3=3`, `m_5=4`, `m_7=9`, `m_9=12` per `solution.md`; `m_13=20`, `m_15=27`
per `NOTE_2026-09-02_EXTERNAL_N12_N15_CERTIFICATE_AUDIT.md`):

q  | n  | block order | block Phi | m      | match
---|----|-------------|-----------|--------|------
5  | 6  | 3           | 3         | m_3 3  | yes
9  | 10 | 5           | 6         | m_5 4  | no
13 | 14 | 7           | 9         | m_7 9  | yes
17 | 18 | 9           | 12        | m_9 12 | yes
25 | 26 | 13          | 30        | m_13 20| no
29 | 30 | 15          | 35        | m_15 27| no
37 | 38 | 19          | 47        | —      | —

At `q = 13`: `Phi(C14) = 21 = m_14` (recorded; `evidence/n14_kflip.json` shows
no undercut at Hamming distance <= 4), so by §1 the optimal order-14 signing is
a two-block signing `[[A,B],[B^T,-A]]` with `Phi(A) = 9 = m_7`. Same at `q = 5`
(block 3 = m_3, C6 optimal: `Phi(C6) = 5 = m_6`). At `q = 17` the block is
`12 = m_9` (order-18 optimum unknown). An independent brute-force switching
witness for the `q = 13` instance is `evidence/paley_two_block_c14_witness.json`.

## 4. Related finite facts (same session; certifiers in scratch)

- **All balanced splits decompose at q <= 13** (complete scans up to
  switching+permutation, self-tested canonical form): C6 10/10, C10 126/126,
  C14 1716/1716; sampled hits independently brute-force realized.
  `scripts/paley_two_block_scan.py` (slow, ~8–10 min) reproduces the scans and
  the witness search; results JSON retained in scratch.
- `Phi(C18) = 33`; `Phi(C26) = 65 = (1/2)*26*sqrt(25)`, the spectral bound
  saturated (the broad-campaign record 61 at order 26 beats it, so C26 is not
  optimal). Both exact enumerations.
- Two-block family exactness at doublings `n = 2..7` (search-based; scratch
  scouts, not landed).

## 5. Scope

Structure only. This does not prove the multiplier-2/-3 comparisons
(`NOTE_2026-09-01_ORIGINAL_LIMIT_TWO_RAY.md`,
`NOTE_2026-09-01_TETRAHEDRAL_TRIPLING_FRAME.md`), the cross-order criterion
hypotheses, or convergence. Natural next targets: (a) the same decomposition
for other conference classes; (b) characterize when the character-split block
attains the half-order minimum (true at q = 5, 13, 17; false at q = 9, 25, 29
in the computed range); (c) use the decomposition as a structural handle on the
doubling program.
