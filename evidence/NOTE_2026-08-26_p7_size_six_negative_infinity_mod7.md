# P7 negative infinity-plus-five closure

Proposition 15.659 closes the `p=7`, `c_H=-1` six-point boundary branch
containing infinity. It does not close the six-finite branch, any `p=5`
size-six branch, larger boundaries, residual (ii), R1, or the limit.

For this branch every direction has phase one. The exact type sum is 32,
and Proposition 15.658's edgewise argument makes all four same-type scaled
slacks congruent modulo eight. The directional floors are 6 for odd-fibre
sizes one and five and 14 for size three. Thus at most one size-three
direction occurs per type, and each type has exactly one mean-14 direction:
it is forced when a size-three direction exists, while otherwise exactly one
of the four mean-six floors is elevated by eight.

The exact odd-fibre Johnson catalogs have sizes

- `b=1,5`, phase one, mean six: 1;
- `b=1,5`, phase one, mean fourteen: 1764;
- `b=3`, phase one, mean fourteen: 36.

Two complete integer floor sweeps—CUDA on the V100 and NumPy on Soulkiller's
CPU—agree on the same 83,496 survivors among all `C(49,5)=1,906,884`
boundaries, including the survivor-list hash
`06d2a7d1ba850347d6c876d551cf3822d01c2e6fc52f839833db8b448c329cd0`.
A 2077-second serial NUKA enumeration and a V100-seeded NUKA quotient agree
exactly on all 1,750 square-semilinear orbit records and their total weight.

For each orbit and each permitted pair of elevated directions, the common
282-by-1225 score system was reduced by its 135-dimensional left kernel over
`F_7`. Affine spans of the two Johnson catalogs reject 2,205 of 2,230 cases.
The remaining 25 cases all have catalog shape `36×36`; direct comparison of
all 32,400 exact catalog pairs leaves zero compatible syndromes. NUKA and
Soulkiller independently reproduce both stages with identical mathematical
records.

The compact audit is
`evidence/p7_size6_negative_infinity_audit.json`. Raw recordings and the
exact scripts are archived at
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p7-size-six-negative-infinity/`.
