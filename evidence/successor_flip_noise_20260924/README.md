# Prepared remote checks — not yet executed

The two new proof notes have analytic author review only. There is no new
remote result log, checker exit status, or independent proof review.
NUKA's SSH preflight failed with `socket: Operation not permitted`.
No mathematical test was run on the controller during this follow-up.

These checkers verify scalar identities, rational comparisons, and print
60-digit numerical evaluations. The decimal evaluations are not certified
interval arithmetic. They do not formally verify Gaussianization,
probabilistic inequalities, regularization, or convergence. The original
convergence problem remains OPEN.

## Pinned source bytes

| File, relative to repository root | SHA-256 |
| --- | --- |
| `evidence/both_phase_disagreement_20260924/check.py` | `34f10712d650e540ee336482b9fe6189a30e3d171996f56413f71fbee63e3783` |
| `evidence/successor_flip_noise_20260924/check.py` | `d6f3654060c0e0fdced11ad9876c349d4048cfc2366f1e3425decef7bee154bb` |
| `evidence/NOTE_2026-09-24_BOTH_PHASE_DISAGREEMENT.md` | `bdf0130e0e4850feb6a2297bfc43e1c1d34567827d7f70755e5c4747c02831b6` |
| `evidence/NOTE_2026-09-24_SUCCESSOR_FLIP_NOISE.md` | `ad9d92eab491112faba320326ae09714f95a292eaa540e14c2d1fdc075d523eb` |

## Manual remote execution

Use NUKA CPU with one serial process, after placing these new files in a
remote checkout. They have not been pushed automatically. Requires Python
and SymPy already installed. Run from that checkout's root:

```bash
bash -eu - <<'SH'
qmm_stage=$(mktemp -d /tmp/qmm-noise-d6f3654060c0.XXXXXXXX)
cp evidence/both_phase_disagreement_20260924/check.py "$qmm_stage/both_phase_check.py"
cp evidence/successor_flip_noise_20260924/check.py "$qmm_stage/successor_check.py"
test "$(sha256sum "$qmm_stage/both_phase_check.py" | cut -d' ' -f1)" = 34f10712d650e540ee336482b9fe6189a30e3d171996f56413f71fbee63e3783
test "$(sha256sum "$qmm_stage/successor_check.py" | cut -d' ' -f1)" = d6f3654060c0e0fdced11ad9876c349d4048cfc2366f1e3425decef7bee154bb
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
printf 'Retained staging directory: %s\n' "$qmm_stage"
for qmm_name in both_phase successor; do
    qmm_rc=0
    python3 -u "$qmm_stage/${qmm_name}_check.py" >"$qmm_stage/${qmm_name}.log" 2>&1 || qmm_rc=$?
    cat "$qmm_stage/${qmm_name}.log"
    printf '%s exit status: %s\n' "$qmm_name" "$qmm_rc"
    test "$qmm_rc" -eq 0
done
SH
```

Retain the staging directory and both logs. An absent second log means the
first check failed or setup stopped; do not report that both passed.
No GPU, matrix simulation, signing census, or previous unchanged check is
requested. Read the complete new proof notes and scripts before running.

The old Orin scalar check has separate provenance at
`../paired_disagreement_floor_20260924/USER_REPORTED_CHECK.md` and does not
cover these new files.
