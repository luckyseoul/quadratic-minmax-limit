# Parity-ceiling scripts

Support for `evidence/NOTE_2026-11-11_CONFERENCE_PARITY_CEILING.md`.
Receipts are in `evidence/parity_ceiling_20261111/`. The original MO limit is
OPEN; nothing here bounds `m_n` from below or makes a cross-order claim.

| file | what it does |
|---|---|
| `paley_lib.py` | Paley conference matrix (prime `q`), quadratic character, `q = a^2 + b^2`, local fields |
| `ceiling.py` | closed-form and integral parity ceilings for every symmetric conference matrix of order `q+1`; admissible values of `Q` with their field multisets |
| `torus.py` | explicit non-split-torus witness, exact check of `x^T C x = (q+1) b` |
| `sym.py` + `qopt.c` | best Paley state fixed by an automorphism of prime order (orbit quotient; exact Gray code up to 34 orbits, annealing above) |
| `defect.c` | exact decision: does `C_{q+1}` have a state with a prescribed field multiset (sparse defect-vector search) |
| `sa.c` | annealing lower bounds for `Phi(C_{q+1})` (witnesses only) |
| `quartic.py` | square-wave states on quartic classes: Paley value `(q-1)a/2+1`, lift value `(m-1)(a+b)+1` (note of 2026-11-12) |
| `windows.py` | rigorous windows `L <= Phi(C_{q+1}) <= U` from small-orbit quotients and the parity ceiling |

Run from the repository root:

```sh
python3 scripts/parity_ceiling/ceiling.py                      # ceiling table (seconds)
python3 scripts/parity_ceiling/torus.py 1000                   # 80 primes, exact
gcc -O2 -o /tmp/qopt   scripts/parity_ceiling/qopt.c   -lm
gcc -O2 -o /tmp/defect scripts/parity_ceiling/defect.c
gcc -O2 -o /tmp/sa     scripts/parity_ceiling/sa.c     -lm
QOPT=/tmp/qopt python3 scripts/parity_ceiling/sym.py 61        # all four quotients exact, ~1 min
/tmp/defect 61 7  2 1  -1 2  1 9                               # pattern {3:2, 7:50, 11:9, 15:1}: 0 solutions, <1 s
SECOND_SIGN=1  /tmp/defect 61 7  -1 3  1 12                    # pattern {3:3, 7:47, 11:12}, half 1
SECOND_SIGN=-1 /tmp/defect 61 7  -1 3  1 12                    # half 2 (minutes each)
/tmp/sa 61 30000 1500 101                                      # q restarts sweeps seed
python3 scripts/parity_ceiling/quartic.py 1000                  # both square-wave identities, exact
QOPT=/tmp/qopt python3 scripts/parity_ceiling/windows.py 5 1000  # windows table, ~25 min on two cores
python3 -m pytest tests/test_parity_ceiling_note.py tests/test_quartic_resonance_note.py -q
```

`defect` arguments: `q c t1 n1 t2 n2 ...` where `c` is the modal field value
and the remaining fields are `c + 4 t_k` with multiplicity `n_k`, listed in
the order the vertices are placed (tightest role first). `FIRST=1` stops at
the first solution.
