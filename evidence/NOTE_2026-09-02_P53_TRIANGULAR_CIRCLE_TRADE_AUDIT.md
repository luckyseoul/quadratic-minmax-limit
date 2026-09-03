# First-live-prime triangular circle-trade audit

**Date:** 2026-09-02
**Status:** exhaustive finite diagnostic; residual (ii) remains open

This audit asks one narrow question left by the positive-mismatch fibre.  If
the Kiss--Somlai Boolean shadow `x` can be changed to another completion of
the same signed triple by flipping one square Miquelian circle `S` disjoint
from the three spikes, then `C_ij x_i x_j=+1` on every pair in `S`.  In that
case `-x 1_S` is a sparse `+p` eigenvector.  Its intersection with the spike
circle is exactly the mismatch created by the trade.

The shardable exact implementation is
`scripts/residual_three_spike_circle_trade_scan.py`.  Four machines split the
first live prime `p=53` by the 2,809 possible finite circle centres:

| worker | centre range | circles checked | switchable circles |
|---|---:|---:|---:|
| soulkiller | `[0,703)` | 18,278 | 0 |
| NUKA | `[703,1405)` | 18,252 | 0 |
| jellyfin | `[1405,2107)` | 18,252 | 0 |
| orin | `[2107,2809)` | 18,252 | 0 |
| **total** | `[0,2809)` | **73,034** | **0** |

Every candidate first passes a base-row signed-clique test and every hit is
then checked on all pairs.  Regression cases show that this is not a vacuous
detector: at `p=11` it finds four switchable circles, including one producing
circle mismatch two; at `p=13` it finds eight switchable circles, all with
mismatch zero.

The four raw shard SHA-256 hashes, in the table order, are

```
4ff1426bbe86d1ec91e731f4876fd48a4cb77226f95fcc66a7c04a14d7e8cb38
ddc1d18755ba0f36c0c3a13788c21d017a2cff4bbc79d7f0be7c3813bf6b8e9f
c71042fee64835ccbbf7b0bf647c878930421a5bffc77daa5d99dca98536ef0a
4d66234630aa9d465921ce86aa91d4daf1b5e6d27132a3c253311750f779d51b
```

Thus the simplest one-circle trade does **not** populate the positive-
mismatch fibre at `p=53`.  This is not a classification of all Boolean
completions: a completion could differ by a more complicated integral
`+p` eigenvector.  The result therefore closes one concrete construction
mechanism, not residual (ii) or the whole positive-mismatch fibre.
