# Cross-transfer / C_cross sweep scripts (2026-09-18)

These scripts lived in `/home/nick/mo_*.py` on soulkiller and were
reported missing from the repo in
`HANDOFF_2026-09-18_cross_transfer_assessment.md`. They are measurements,
not proofs. Copied here 2026-09-19 so the n=13 2n-transfer numbers can
be replayed.

| File | Role |
|------|------|
| `mo_xterm_test.py` | A-term at a cross-maximizing pair |
| `mo_cross_const.py` | C_cross, sign field vs n |
| `mo_cross_struct.py` | C_cross, structured Gaussian field |
| `mo_cross_sweep.py` | alpha/rho sweep |
| `mo_cross_fine.py` | fine alpha/rho sweep |
| `mo_aterm_verify.py` | A-term at (alpha=a, rho=0.995) |
| `mo_fullmax.py` | full 2n max at n=13 |

The 3n script mentioned in the handoff (`mo_3n.py`) was not present on
disk. Do not treat the n=13 max 0.898 < 0.9216 as a bound on `m_{2n}`.
The recommended continuation is residue (6.20), not a larger sweep.
