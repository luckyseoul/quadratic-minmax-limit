# Φ(C42) = 129: the random-orthogonal law's first prediction is falsified

2026-11-11. Follow-up to
`NOTE_2026-11-11_RANDOM_ORTHOGONAL_MODEL_LAW.md` (§5, prediction 1) and the
GPU-dispatch suggestion in the 2026-11-11 HANDOFF entry. The original MO
limit is OPEN. Nothing here bounds any `m_n`; this note records one exact
value and one failed prediction of the empirical law.

## Result

For the Paley conference matrix `C_42` (q = 41, order n = 42),

    Phi(C_42) = 129,

exactly, over all `2^41` projective states. The maximizing and minimizing
values of `x^T C_42 x` are `+258` and `-258`; at even n the parity of
`x^T A x` is fixed by `n(n-1) mod 4 = 2`, so 258 is admissible.

This is one parity step above `largest_admissible_below`, and it is a
**stop-rule event** in the source note: §5 states "`Phi(C42) >= 129` ... would
count against the law." The predicted value was 127 or 125.

Context values (`S_n = (1/2) n sqrt(n-1)` for n = 2 mod 4):

    S_42            = 134.46560898608982
    Phi / S_42      = 0.9593531087442945
    alpha_42        = 129 / 42^(3/2) = 0.4739317177407108
    Haar * S_42     = 128.9 (source note, interpolated mean)
    Haar * S_42 - Phi = -0.1

The magnitude of the law is close (0.07% in alpha); the parity-step
prediction is wrong. Per the source note this counts against the law, not
against the open limit.

## Certificate independent of the enumeration

A projective state with bitmask `67778113536` (bit i set means spin i is -1;
`x_41 = +1`, 18 negative spins) has, in exact Python integer arithmetic,

    x^T C_42 x = 258,    so Q = 129 and Phi(C_42) >= 129.

Thus the falsification rests on one exactly checked state; the full-cube
enumeration only establishes that 129 is attained on both signs and that no
larger value exists.

## Enumeration and validation

Sharded Gray-code walker (`scripts/rom_law/gs2s.c`, built from the branch's
`gs2.c`): each shard fixes the 12 low free bits and walks the remaining
`2^29` states; `4096` shards x `2^29` = `2^41` projective states, 88 threads,
**567 s wall**.

End-to-end validation of the harness on the known values:

- C14 -> 21, C18 -> 33 (small, full cube);
- C38 -> 109 over `2^37` states, 256 shards, 37 s (matches the September
  2026-09-12 evaluation);
- the branch's `paley.py` builder is entrywise identical to the repo's
  `scripts/paley_structural_lift_family.py` at q = 13 and q = 41
  (`C C^T = q I` checked).

The full-cube value is attained in 376 of the 4096 slices on each sign,
consistent with a codimension-1 maximum rather than a single accidental
witness.

## Scope

- One exact finite value; no asymptotic claim and no bound on `m_42`.
- The sharded walker was written for this test. It reproduces the branch's
  own `gs2` on all validated cases; it is not a new proof method.
- The source note's other predictions (items 2-5) are untouched.
- `evidence/rom_law_20261111/c42_result.json` is the machine receipt.

## Repro

    cd scripts/rom_law
    gcc -O3 -march=native -funroll-loops -o gs2s gs2s.c
    python3 -c "import sys; sys.path.insert(0,'.'); from paley import paley; \
        import numpy as np, struct; C=paley(41); f=open('C42.bin','wb'); \
        f.write(struct.pack('i',len(C))); f.write(np.ascontiguousarray(C,dtype=np.float64).tobytes())"
    seq 0 4095 | xargs -P 88 -I{} ./gs2s C42.bin 12 {} > shards.txt
    # Phi = max(max col 1, -min col 2) / 2

Witness recheck (no enumeration):

    python3 -c "import sys; sys.path.insert(0,'.'); from paley import paley; \
        import numpy as np; C=paley(41).astype(np.int64); m=67778113536; \
        x=np.array([(-1 if (m>>i)&1 else 1) for i in range(42)],dtype=np.int64); \
        print(int(x@C@x))"
    # 258
