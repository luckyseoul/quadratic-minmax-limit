# Conditional-cross spectral bridge: exact order-four diagnostic

Classification: exhaustive finite diagnostic only.  This does **not** prove
an all-orders bound on a conditional cross minimizer's operator norm and does
not prove the original limit.

## Exact scope

For the order-four `Phi=4` source used in
`scripts/diagonal_majorizer_gap_scan.py`, the script
`scripts/conditional_spectral_bridge_probe.py` enumerates all `2^16=65,536`
cross signings `B`.  For each it exactly evaluates

\[
 \max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TBy|\bigr).
\]

The conditional minimum is `10`, attained by `184` cross signings.  Their
operator-norm histogram is:

| `||B||op` | count |
|---:|---:|
| `2` | 16 |
| `2.61312592975` | 64 |
| `2.73205080757` | 96 |
| `2 sqrt(2)` | 8 |

Thus this source has no order-four conditional minimizer in the trivially
large-operator-norm regime.  It leaves open the necessary asymptotic bridge,
for example `||B_n||op=O(sqrt(n))` for conditional minimizers relevant to the
joint-shell mismatch argument.

## Independent replay

The deterministic JSON output has SHA-256
`0b23ae0975aa6063282f097b0221fb95a9feb876703f18d910b7d47ba22e354f` on both
the controller and Soulkiller.  The script SHA-256 is
`7345b3b077c94c29ea1970fc6c0ba9f91acc6a7a7b1fa353481fe9724b1d4eb5`.
