# Finite negative two-point exclusions at `p=7,11,13`

Date: 2026-08-25. This is Proposition 15.648. After Proposition 15.647
closes all odd primes `p>=17`, this proposition closes the negative-product
two-point branch at `p=11,13` and four of the five arithmetic profiles at
`p=7`. It leaves exactly the balanced `p=7` profile and the separate `p=5`
case. Other boundaries and every top-level gate remain open.

## 1. Exact `p=13` l1 contradiction

The two Proposition 15.647 profiles are `(x,y)=(0,6),(6,0)`. In either
orientation, a zero-baseline direction has `I=9`, `E=44`, and

\[
 K_{st}=\pm(1-w_s-w_t),\qquad w_s\ge0,\qquad\sum_s w_s=10. \tag{1}
\]

Let `u` be the number of positive entries of `w`; then `1<=u<=10` and
`z=13-u` entries vanish. Separating zero-zero, zero-positive, and
positive-positive pairs gives

\[
 \sum_{s<t}|1-w_s-w_t|
 ={z\choose2}+z(10-u)+10(u-1)-{u\choose2}
 =u^2-25u+198.                                         \tag{2}
\]

This decreases on `1<=u<=10` and has minimum `48` at `u=10`. But a
zero-baseline direction has all `E=44` finite edges transverse, so its
inter-fibre l1 norm is at most `44`. The contradiction closes both `p=13`
orientations without computation.

## 2. Exact finite model

For `p=7,11`, the CP-SAT model uses one Boolean variable per finite edge and
one per infinity-star endpoint. It imposes:

- the exact total, baseline, and exceptional direction counts from 15.647;
- the boundary equation `boundary(F)=S symmetric-difference {v}` via native
  XOR constraints;
- `c_H=-1`;
- every exact baseline identity
  `K_st=epsilon_d(2c-delta_s-delta_t-n_s-n_t)`; and
- the redundant exact l1 edge budget for propagation.

All constraints are integer equalities, XORs, or cardinality inequalities.
Thus `INFEASIBLE` is an exact finite nonexistence certificate for the model;
no floating-point tolerance or objective bound is used.

## 3. Exceptional-pair symmetry

Square multiplications of `F_{p^2}` and Frobenius preserve quadratic
direction type and the complete model. Their action on opposite-type
exceptional pairs has:

- two orbits of size 8 at `p=7`, represented by `(0,1),(0,3)`;
- three orbits of size 12 at `p=11`, represented by `(0,1),(0,2),(0,3)`;
- four orbits at `p=13` (not needed for the analytic close).

The orbit computation is reproduced exactly in `e1_gmin_m4_prop15648.py`.

At `p=11`, all three representatives are infeasible for each orientation
`(x,y)=(0,5),(5,0)`, closing that prime. At `p=7`, direct all-pair sweeps
exclude `(0,3),(0,6),(3,0),(6,0)`. The balanced `(3,3)` profile timed out
and is not claimed empty.

## 4. Evidence archive

The ten raw JSON outputs are archived at

`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-25-negative-two-point/`.

Their SHA-256 hashes are recorded in
`evidence/e1_gmin_m4_prop15648.json`. The generating model is
`scripts/residual_negative_full_cpsat.py`.

## Boundary of the result

The remaining `c_H=-1`, `D={infinity,v}` cases are:

- `p=5`, where the three-exception bound does not guarantee a baseline in
  each three-direction type;
- `p=7`, baseline counts `(x,y)=(3,3)`.

The positive-product cases `p=5,7,11,13`, other nonempty boundary profiles,
residual (ii), Type I, R1, global QVAR, and the limit also remain open.
