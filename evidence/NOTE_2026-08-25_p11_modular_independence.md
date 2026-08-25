# The known p=11 modular data do not determine the R1 coefficient

Date: 2026-08-25. This is Proposition 15.641. It closes one proposed
coefficient-determination route, not R1 itself.

## Exact result

At `p=11`, the relevant Kohnen subspace of weight `69/2` and level `44`
has dimension `66` (inside the full dimension-199 modular-form space; Sturm
bound `200`). The currently justified geometric and shell constraints are:

- infinity coefficients `0..19`;
- the half-cusp gap `0..14`;
- the cusp-0 gap `0..5`;
- the cusp-1/4 gap `0..23`;
- the cusp-1/11 gap `0..4`;
- the complete second dual shell, infinity coefficient `20`.

The constraints before the second shell have exact rank `29`; the second
shell adds one rank, giving rank `30` and a 36-dimensional residual kernel.
On that kernel, the second-shell and half-cusp target rows have joint rank
two. In particular, fixing the second shell does not fix the target.

The stronger certificate is one exact 66-coordinate vector `w` with 21
nonzero coordinates:

```text
A_known w = 0
c_second w = 0
c_target w = 1
```

The largest numerator and denominator exponents in this rational witness
are 534 and 499. The binary certificate is stored at

```text
/mnt/storage/e1work/maxplus_p13/r1_modular_attack_2026-08-25/
  p11_modular_independence_witness.gpbin
```

with SHA-256

```text
5bdf184e653079c361f6ee1a2178dd3f4e9b051d9da6625cc3a4910a93b441e7
```

The full 955 MB exact input/cache set is backed up in the same directory;
sizes and hashes are recorded in `evidence/r1_p11_modular_cache_manifest.json`.

An independent PARI read on Soulkiller verified that the payload has 66
coordinates, zero residual, and target `[1]`.

## Meaning and scope

This proves that the known shell coefficients and geometric cusp gaps do
not determine or sign the R1 target through linear modular-form algebra.
The exact second-shell channel is independent of the target channel. It
therefore prevents the second-shell computation from being misreported as
a modular close.

It does **not** produce a theta series with negative representation counts,
does not use the positivity cone of lattice theta series, and does not
refute R1. A successful modular proof would need genuinely new information:
additional complete shells, additional justified cusp data, or nonlinear
positivity specific to the lattice/coset theta series.

## Reproduction

- build the exact space: `scripts/r1_p11_kohnen_cache.gp`;
- normalize its basis: `scripts/r1_p11_kohnen_reduce.gp`;
- generate the exact half-cusp anchor and rational block:
  `scripts/r1_p11_mf2fix_smoke.gp` and
  `scripts/r1_p11_kohnen_half_batch.gp` (the guarded PARI 2.17.3 fix is
  `patches/pari-2.17.3-mf2gaexpansion-relative-zero.patch`);
- generate the other exact cusp blocks with
  `scripts/r1_p11_cusp_batch.gp` (`P11_CUSP=zero`, `quarter`, and `p`);
- audit second-shell independence: `scripts/r1_p11_second_shell_rank.gp`;
- emit the exact nullspace certificate:
  `scripts/r1_p11_modular_independence_witness.gp`;
- independently read the backed-up payload:
  `scripts/r1_p11_verify_modular_witness.gp`;
- machine-readable metadata: `src/e1_gmin_m4_prop15641.py`.
