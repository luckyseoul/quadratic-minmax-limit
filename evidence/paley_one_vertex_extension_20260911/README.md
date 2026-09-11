# Exact one-vertex extension of the order-14 Paley signing

This finite experiment evaluates the exact identity

\[
\min_{b\in\{\pm1\}^{14}}\max_{x\in\{\pm1\}^{14}}
  \bigl(|Q_C(x)|+|b\mathbin\cdot x|\bigr),
\]

for the `q=13` Paley conference signing `C` of order 14. The source has
`Phi(C)=21`. Global sign symmetry permits fixing the first coordinate of
both `b` and `x`, leaving 8,192 incident signings and 8,192 source states;
the programs enumerate all of those representatives exactly with integer
arithmetic.

Soulkiller ran the producing scan with 88 OpenMP workers and obtained 27.
NUKA independently ran a separate exhaustive reduction with 16 workers and
obtained the same source and extension values. NUKA also directly replayed
the selected `b_mask=0` witness against all 16,384 source states.

| result | SHA-256 |
| --- | --- |
| `soulkiller_result.json` | `ad39f64fb7e19be96829e660b3bf4605498b6fa18e7ff07ccc0cced6043f35ee` |
| producer | `c411988b7bbd7b8dd339ff7cd5e8753e69d836a34025886036cf5f90b8616300` |
| independent exhaustive checker | `02483b7586ce48ba4ef9890921e670086d2c55f330e041158a98d495f276e0c8` |
| selected-witness checker | `c484fe7e138f7c771f44868b5d8930ddeaff3d7258128f741ccc03a048f8db9d` |

## Interpretation boundary

This establishes only the stated order-14, source-specific minimum. It does
not show that `C` minimizes the order-14 problem, determine `m_15`, control
the same expression uniformly in `n`, or prove convergence of `m_n/n^(3/2)`.
It is a direct finite test of the one-vertex mechanism needed by a possible
order-to-order argument.
