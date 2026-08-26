# `p=5` four-point exclusion from complete eigenshell syndromes

Proposition 15.656 closes every `p=5` size-four residual boundary. Combined
with Propositions 15.652--15.655, every size-four boundary is now excluded
for every odd prime `p>=5`. This does **not** close residual (ii): boundaries
of size at least six remain, as do R1, global QVAR, Type I, and the limit.

## 1. Complete shell equations

At `p=5` there are 26 vertices and 325 edge indicators `x_e`. A residual
graph has 21 edges and contains the distinguished edge. For an eigenshell
sign `eps`, quotienting the full shell by `y~-y` leaves 130 distinct edge
feature rows

```text
f_y(e) = eps y_a y_b C_ab in {+1,-1}.
```

Every normalized edge column has sum 26. If all normalized scores are at
least three, write

```text
sum_e f_y(e)x_e = 3 + 2 A(y),       A(y)>=0.
```

Summing over the 130 representatives gives the exact shell mass

```text
sum_y A(y) = (21*26 - 130*3)/2 = 78.                (1)
```

The four-point odd-degree boundary `D` and Paley edge-product sign `c_H`
fix the parity vector:

```text
(-1)^A(y) = -eps c_H product_(v in D) y_v.
```

Thus `A=P+2L`, with `0<=L<=4` and

```text
sum_y L(y) = (78-sum_y P(y))/2.                     (2)
```

For one shell, edge count, the distinguished edge, and all 130 bad-edge
counts form a common `132 x 325` zero-one matrix. Exact finite-field
elimination gives

```text
rank_F5 = 67,        left-null dimension = 65.
```

Substituting `9-P-2L` for each bad-edge count converts those dependencies
into 65 bounded congruences in the 130 lift variables. Failure of this
bounded system is a rigorous necessary-condition exclusion of the original
integer edge model. The combined two-shell matrix is `262 x 325`, with
rank 113 and 149 left dependencies over `F_5`; in the completed scan every
exclusion already occurs in one shell.

## 2. Complete orbit exhaustion

Proposition 15.632's exact parity-floor filter leaves the following
square-semilinear orbit cases:

```text
product sign   infinity in D   boundaries   orbits
-1             no                  10925      489
-1             yes                  2300      112
+1             no                  10925      489
+1             yes                  2300      112
                                      ----     ----
                                     26450     1202
```

The direct exact scans use all negative-sign cases and the positive-sign
infinity-present cases:

```text
direct orbit cases                         713
shell-local mod-5 infeasible               712
mod-5 unknown                                1
```

The one mod-five timeout is negative sign, no infinity, orbit 164,
represented by boundary `[2,3,12,13]` and having orbit size 24. In its
positive eigenshell, `sum P=56` and `sum L=11`. Rebuilding the `132 x 325`
matrix independently over `F_7` again gives rank 67 and 65 dependencies.
The bounded lift system is exactly `INFEASIBLE`. Solver `UNKNOWN` is never
counted as an exclusion.

Consequently all 713 direct orbit cases are closed, covering 15,525
boundaries.

## 3. Exact transfer of the remaining sign

Let `alpha` be a nonsquare of `F_25`. Multiplication `u->alpha u`, together
with switching only the infinity coordinate, fixes infinity, finite zero,
and the distinguished edge and satisfies

```text
S C[pi,pi] S = -C.
```

It therefore exchanges the two eigenshells while preserving every
normalized score. For a 21-edge graph whose boundary omits infinity,
the infinity degree is even, so

```text
product_(e in pi(H)) C_e
  = (-1)^21 (-1)^deg_H(infinity) product_(e in H) C_e
  = -product_(e in H) C_e.
```

Hence the negative no-infinity exclusion transfers bijectively to the
positive sign. A fresh orbit reconstruction verifies a 489-orbit,
10,925-boundary bijection. When infinity belongs to the boundary the product
does not flip, which is why both signs of that shape were solved directly.

Together with Proposition 15.632's 3,450 floor exclusions, this covers all
`2*binom(26,4)=29,900` `p=5` boundary/sign cases.

## 4. Audit and scope

The structural audit independently rebuilds all four boundary
classifications, both shell matrices modulo 5 and 7, every recorded
parity/lift mass, finite-solver row scope, the unique exception key, and the
nonsquare orbit bijection. It reports all checks true, 712 mod-five
exclusions, one mod-seven exclusion, zero feasible case, and zero remaining
unknown.

This proves:

```text
p=5 size-four boundary                         CLOSED
all size-four boundaries for odd p>=5         CLOSED
first open non-Walsh boundary size for p>=5       6
residual (ii), R1, global QVAR, Type I, L      OPEN
```

## 5. Literature and OEIS context

Post-finding searches for the exact matrix dimensions/ranks and the
Paley/full-shell bounded-syndrome construction found no matching theorem.
Ghinelli--Key (2011) studies codes from ordinary Paley graph and line-graph
incidence matrices; its matrices and ranks are different from the
`132 x 325` full-eigenshell score system used here.

Individual OEIS searches find `26450`, `15525`, and `10925` in unrelated
partition-count sequences. No search matched the structural tuples
`132,325,67,65`, `262,325,113,149`, or `1202,713,712,1`. These are context
checks only; no sequence novelty or priority claim is made.

Ivanisvili--Stolyarov--Vasyunin--Zatitskii, arXiv:2305.03523, was also read
as requested. Its minimal-locally-concave Bellman construction on planar
non-convex domains is a plausible conceptual guide for a two-moment R1/QVAR
compression, but it supplies no theorem for this finite Paley shell system
and is not used in Proposition 15.656.

## 6. Reproduction and permanent archive

Core programs:

- `scripts/p5_size_four_full_shell_mod5_batch.py`;
- `scripts/p5_size_four_full_shell_mod7_exception.py`;
- `scripts/p5_size_four_full_shell_audit.py`;
- `scripts/residual_size_four_boundary_orbits.py`.

Permanent archive:

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/
  2026-08-26-p5-four-point/
  p5_size_four_full_shell_certificate_2026-08-26.tar.gz
```

Hashes:

```text
archive         d5db5e82389ebb0bfcb23e80da5e2322b1d65e74aa8f3804d25275793b7380da
audit           5cafab9272510dc6871818fbfc395c8f3386ca6615f7a1b1e36f4785cf1d7e4f
mod-7 exception 2b91f83aeb543a29b0b4398243115f25189d233c50d41e16f8fd14b121c61b06
```
