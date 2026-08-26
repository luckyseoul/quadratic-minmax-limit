# `p=7` infinity-plus-three-finite boundary exclusion

Proposition 15.653 closes every size-four boundary of the form
`D={infinity,u,v,w}` at `p=7`.  Proposition 15.652 already closes the
negative edge-product sign.  The positive sign is reduced to a unique
minimum-mean directional slack model and 416 exact boundary-orbit
certificates.

At Proposition 15.653, four finite points at `p=7` and all `p=5` size-four
cases remained; they are subsequently closed by 15.654--15.656. Larger
boundaries, residual (ii), R1, and the limit remain open.

## 1. Every direction saturates

For infinity plus three finite points, every direction has either one or
three odd boundary fibres.  At `p=7,c_H=+1`, the parity phase is zero and

```text
F(7,1,0) = F(7,3,0) = 8.
```

Each quadratic direction type has four directions and exact budget 32.
Consequently every directional scaled mean is exactly eight.

If there is one odd fibre `B={j}`, pointwise parity and equality of the mean
force

```text
A(X)=x_j.
```

If there are three odd fibres, put `t=|X cap B|`.  There are 35 points on
`J(7,4)`.  The parity vector has 16 mandatory ones, while scaled mean eight
means total slack mass 20.  Thus

```text
A = parity + 2L,       sum_X L(X)=2.
```

There are exactly

```text
C(35+1,2)=630
```

ways to place those two lift units, allowing both at one point.  The
degree-at-most-two Johnson evaluation space has rank 21 and a 14-dimensional
left kernel.  Exact integer/rational evaluation of all 630 corrections
leaves one survivor: both units lie at the unique four-set disjoint from
`B`.  Therefore

```text
A(X)=(t-2)^2.
```

## 2. Sparse coefficient equations

Write `z_s=2x_s-1`; on `J(7,4)`, `sum_s z_s=1`.  The exact score targets are

```text
b=1: epsilon_d S_H(z) = 4 + z_j,

b=3: epsilon_d S_H(z)
     = 5 - sum_{s in B} z_s + sum_{s<t in B} z_s z_t.       (1)
```

Let `I` be the infinity-edge count, `n_s` its fibre counts in direction
`d`, `P_d` the finite parallel-edge count, and

```text
K_st = sum epsilon_d C_uv
```

over selected finite edges joining fibres `s,t`.  If the target in (1) is
written as

```text
c + sum_s l_s z_s + sum_{s<t} q_st z_s z_t,
```

equality on `sum z_s=1` is equivalent to the existence of one integer
`k_d` with

```text
P_d  = c + sum_s l_s + 3k_d - I,
K_st = q_st + k_d - n_s - n_t + l_s + l_t.                  (2)
```

For `b=1`, `(c,l,q)=(4,delta_j,0)`.  For `b=3`, `c=5`, `l_s=-1_B(s)`, and
`q_st=1` exactly inside `B`.  Summing the parallel counts in (2), together
with the 29-edge total and the odd infinity degree, restricts

```text
I in {5,11,17,23,29}.
```

The fixed-boundary CP-SAT model selects the 29 edges of `H` directly,
includes the distinguished edge, imposes the requested odd-degree boundary
and positive Paley-sign product, and enforces the 176 coefficient equations
in (2) plus the infinity-count restriction.  The 29-edge total, boundary,
and product constraints are imposed separately.  The coefficient equations
reconstruct every affine score pointwise, so no affine inequality is lost.

## 3. Complete orbit exhaustion

All `C(49,3)=18,424` finite triples survive the parity budget.  Square field
multiplication and Frobenius form a 48-element subgroup fixing infinity, the
finite point zero, the distinguished edge, Paley signs, and direction type.
Complete enumeration gives

```text
18,424 triples = 416 boundary orbits.
```

Every one of the 416 fixed-boundary coefficient models is `INFEASIBLE`.
Two initial 60-second unknowns were rerun with independent seeds and longer
limits; both became `INFEASIBLE`.  Final totals are

```text
INFEASIBLE 416, UNKNOWN 0, FEASIBLE 0.
```

The independent audit re-enumerates all triples and all 416 orbits, reruns
the exact 630-correction slack classification, checks every result row
against its orbit representative and model scope, and verifies complete
infeasible coverage.

## 4. Archive

Permanent archive:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/
  2026-08-26-p7-four-point/p7_infinity_three_certificate_2026-08-26.tar.gz
SHA256 f9a125a18d287eef63e579b8416022e0c4d91dee0a82334fabb1167cc9356c17
```

Independent audit SHA-256:

```text
d974d27274aacf8987b954a42500d7390bb9f9b895f8b813b3c7faff323264c5
```

Orbit-source SHA-256:

```text
0b3a928e98e11838eac051c768a3b4aa0f3d5e5d32dd17fde8f8021016eae941
```

Core scripts:

- `scripts/p7_size_four_slack_classify.py`;
- `scripts/residual_size_four_boundary_orbits.py`;
- `scripts/residual_boundary_four_lift_cpsat.py`;
- `scripts/residual_size_four_orbit_batch.py`;
- `scripts/p7_size_four_certificate_audit.py`.
