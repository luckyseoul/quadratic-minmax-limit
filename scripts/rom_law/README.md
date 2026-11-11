# rom_law

Exact and annealed ground states of Haar-conjugated flat spectra, and exact
`Phi` for Paley lifts and restrictions. See
`evidence/NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md`.

    gcc -O3 -march=native -funroll-loops -o gs2 gs2.c
    python3 paley.py              # builder check: C6, C14, C18 -> 5, 21, 33
    python3 rom.py N DRAWS        # exhaustive Haar row (N even, N <= 30 practical on CPU)
    python3 rom_sa.py             # annealing rows 30..50 (lower estimates)
    python3 lift36.py -1          # K(C18, C18 - I), 2^35 states, two processes
    python3 restrict.py Q d1 d2   # C_{Q+1} minus the first d points

`gs2 file [t]` reads int32 n and n*n float64, enumerates the projective cube
by Gray code (optionally split on x_{n-2}), prints max and min of x^T J x.
