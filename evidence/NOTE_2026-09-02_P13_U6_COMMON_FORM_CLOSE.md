# Proposition 15.754: the `p=13,k=60,u=6` endpoint closes

**Date:** 2026-09-02

**Pinned predecessor HEAD:**
`cf32d2137d9b40ce631f21bcdf8b2cb6e72d0c81`

## Result

The last open fifth-shell endpoint

```text
p=13, t=4, k=60, u=6
```

is empty. This is an exhaustive finite aggregate/common-form certificate and
a proved endpoint theorem. It is not a graph, orbit, coefficient-cell,
support, or common-realization census.

Consequently, with Propositions 15.752--15.753,

```text
k=4p+8 is impossible for every prime p>=13.
```

Residual (ii) is still globally open at critical `p=5,7`, at
`p=11,k>=50`, in later p13 layers, and beyond Proposition 15.752's band.
The separately tracked positive `p=7,z=7` boundary systems also remain open.
Thus E1, `L=1/2`, and the original convergence problem remain open.

## Exact normalization

The hard means and quotient budget are

```text
a_L = 12+14 k_L,       sum_L k_L = 5.
```

At least two hard directions are exact `b=2` XNOR rows. The unspecialized
row sums, coefficient-offset congruence, and common signed edge total give

```text
hT=5,       P_L=4+k_L,       |E_h|=33,       |E_-h|=28.
```

Opposite count `Q=3` has mean eight, below both the phase-zero nonzero floor
12 and the integral-lift floor 10. Hence every opposite count is `Q=4`.
The seven hard-excess partitions, exact-root counts, collision minima, and
nonexact Parseval bases are

| partition | roots | `C_min` | base |
|---|---:|---:|---:|
| `(1,1,1,1,1)` | 2 | 0 | 303 |
| `(2,1,1,1)` | 3 | 0 | 298 |
| `(2,2,1)` | 4 | 0 | 293 |
| `(3,1,1)` | 4 | 1 | 289 |
| `(3,2)` | 5 | 1 | 284 |
| `(4,1)` | 5 | 2 | 276 |
| `(5)` | 6 | 3 | 259 |

The required nonexact energy is `base+26C`.

## Low-root certificate

For the first two partitions, exhaustive coefficient evaluation of the
common binary quadratic `U=hM2` and quartic `G=hM4-M2^2` gives:

| partition | fixed models | coefficient pairs | maximum | targets |
|---|---:|---:|---:|---:|
| `(1,1,1,1,1)` | 42 | 43,184,232 | 293 | 303, 329, 355 |
| `(2,1,1,1)` | 280 | 10,221,120 | 290 | 298, 324, 350 |

Both maxima are attained in the original pinned one-worker CP table model.
The strict deficits close every collision count allowed by the raw energy
ledger.

## Four-root certificate

Four exact roots force

```text
G=R4*c,       J6=h(M6-M2^3)=R4*Q2.
```

The sign-safe hard/opposite keys are

```text
hard:     (N2, N4-N2^2, N6-N2^3)
opposite: (-N2, -N4-N2^2, -N6+N2^3).
```

The exact sixth-moment regression is

```text
W6=(1,12,1,1,12,12),
```

computed with ordinary integer modular powers before NumPy conversion. For
each partition the explicit join checks both signs, all 70 signed root sets,
218,320,284 `(U,G,J6)` coefficient triples, 336 opposite-compatible
survivors, and 1,008 hard-assignment checks.

- `(2,2,1)`: the reachable energy maximum is 193, below targets
  293, 319, 345, and 371.
- `(3,1,1)`: no hard-row assignment survives at targets
  315, 341, 367, or 393.

## High-root certificate

For `(3,2)`, `(4,1)`, and `(5)`, at least five exact XNOR roots force the
binary quartic `G` to vanish identically. Exact one-worker six-bin models
then impose all 74 translated cuts, the sign-correct hard/opposite quartic
relations, and the six-positive/seven-negative collision floor. The final
strict gaps are

```text
(3,2): 22
(4,1): 28, 50
(5):   38, 42, 56, 58, 82.
```

Solver-free fixed-sphere and quartic-character checks independently exclude
the two sharp equality boundaries.

## Replay

```bash
cd /home/nick/quadratic-minmax-limit-residual-ii
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
  /home/nick/.venvs/mo-exact/bin/python -m e1_gmin_m4_prop15754
PYTHONDONTWRITEBYTECODE=1 \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_prop15754.py \
  tests/test_p13_u6_cut_equalities.py \
  tests/test_p13_u6_high_root_energy.py
```

The focused replay collected 32 tests and exited zero on 2026-09-02.

To regenerate the three common-form artifacts independently:

```bash
PYTHONDONTWRITEBYTECODE=1 /home/nick/.venvs/mo-exact/bin/python \
  scripts/p13_u6_low_root_ug_bound.py
PYTHONDONTWRITEBYTECODE=1 /home/nick/.venvs/mo-exact/bin/python \
  scripts/p13_u6_four_root_ugj.py --case-index 0
PYTHONDONTWRITEBYTECODE=1 /home/nick/.venvs/mo-exact/bin/python \
  scripts/p13_u6_four_root_ugj.py --case-index 1
```

## SHA-256 manifest

```text
f3c54009c8d12494aae14885d268f54dd89cffdd5a3730a766e3730e3b7e2f63  evidence/e1_gmin_m4_prop15754.json
9d9a3b75a00410706df94dd02086357393b00bb463b3f9286596a08eb5faa0a4  evidence/e1_gmin_m4_prop15754_low_root_ug.json
83cfbe4f684c4406a2d23216e4049f083c438ed07fda23214bcffbe1c0d2449b  evidence/e1_gmin_m4_prop15754_four_root_221.json
d924059b403f2ddb33282ff60ef01c32ac7ce3806a276083157eb0b073fec2bc  evidence/e1_gmin_m4_prop15754_four_root_311.json
811eb1833d5551f3fadce62ec6a4302296e301a2cd001bedd399df723afaaebe  scripts/p13_u6_joint_ug_tables.py
213a435e6c8848879e1fb485fa9fdd82832c72b17b863090183c1cd65e62bfaa  scripts/p13_u6_low_root_ug_bound.py
71754b7e9e0f5df38819620f5f90903614152ca4a1516e544fa6fd86f4c0aef1  scripts/p13_u6_four_root_ugj.py
```
