# Proposition 15.750: all-prime multi-level Type-I close

## Result

The multi-level Type-I bad-case box is empty for every prime `p>=5` and
every distinguished edge. The live predicate
`type_I_multilevel_bad_case_ND_closed()` is True. Residual (ii), E1, and the
quadratic-minmax limit remain open.

## Uniform argument for `p>=11`

For a hypothetical bad case `G`, set `H=G union {e}` and use the edge
multiset `W=G+2e`. An isolated vertex of `H` is transported to infinity by
the signed PSL action. In each square direction, the exact Max+ identity on
the Johnson middle slice makes every signed off-fibre block sum equal.
The `l1` budget forces that common value to zero and exactly three parallel
units. Summing directions pins the positive/negative multiplicities of `W`.

In a nonsquare direction define `T=(-S_H-2)/2`. Parallel-edge averaging
forces some direction with `2p E[T]` equal to `4` or `6`. The product of edge
features makes `T mod 2` an affine fibre parity. The central Krawtchouk
recurrence gives bias at most `1/p`, so such a low mean forces `T` to be
everywhere even. The nonzero integral quadratic `B=T/2` then contradicts
Proposition 15.688's sharp floor `4p E[B]>=p-3`.

This is a uniform theorem, not a finite-prime or graph census.

## Exact bases

The tracked artifacts are:

- `e1_type_i_badcase_farkas_p5.json`, SHA-256
  `f3e8c8a0f85fcaf95bc5b0556eced1a3699735e3bc11b78285bb4d25abd8008a`;
- `e1_type_i_badcase_farkas_p7.json`, SHA-256
  `40a6e5156817a421dcc7debe75a55af5749ffd6cdcb12cda78ead75a3d0cc8db`.

At `p=5`, 231 positive integer rows give exact right side `-144`; the
certificate uses the cardinality side condition but no upper bounds. At
`p=7`, 1,226 rows give a full nonnegative-cone certificate in all 1,225 edge
variables, with no cardinality or upper-bound premise. In both cases the
leaf verifier regenerates the Paley conference matrix, validates every
stored Boolean eigenvector, and checks `A^T lambda=0` and `b^T lambda<0`
using Python integers. It reads no eigenshell cache and imports no optimizer
or SciPy. Signed-PSL 2-transitivity transports any distinguished edge to the
canonical certified edge `(infinity,0)`.

## Replay

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15750.py
PYTHONPATH=src pytest -q \
  tests/test_type_i_small_prime_exact.py tests/test_prop15750.py
```

The focused suite includes artifact corruption, no-cache/no-SciPy,
all-edge normalization, symbolic-gap, evidence replay, and live-predicate
tests.
