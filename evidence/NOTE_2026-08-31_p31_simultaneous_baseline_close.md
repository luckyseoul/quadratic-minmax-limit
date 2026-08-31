# The p=31 endpoint closes by simultaneous baseline coefficients

**Date:** 2026-08-31
**Result status:** proved symbolic endpoint exclusion
**Scope:** `p=31`, `|H|=125`, all-finite `|D|=32`, outside pair slack `R=10`

**Later status:** Proposition 15.734 transports an isolated vertex to
infinity and applies the same coefficient mechanism without any boundary-size
or slack hypothesis. It subsumes this `p=31` close and excludes every
`|H|=4p+1` residual candidate for every prime `p>=13`. The argument below
remains a valid sharper certificate for the historical endpoint profile.

## Result

The endpoint left open by Propositions 15.727--15.732 is impossible at
`p=31`.  No arc-classification completion and no finite configuration search
is used.  The proof compares the coefficients of the fifteen exact
phase-one mean-30 baseline quadratics forced by Proposition 15.728.

Consequences:

- all fifteen mean-30 directions in the Paley-hard type have `b=2`, not
  merely fourteen;
- a block row with `y` 4-secants has at least `5+y` nonrich hard `b=2`
  directions;
- the first possible positive slack at `p=31` is now `R=11`;
- the first prime whose equality endpoint remains unexcluded is `p=37`,
  where `R=12`.

The full first shell was still open at the 15.733 stage. Proposition 15.734
subsequently closes every boundary at `k=4p` for `p>=13`; residual (ii)
remains open only because even `k>4p` and small-prime remnants are not covered.

## Coefficient lemma

Put `z_s=2x_s-1`, so `sum_s z_s=1` on `J(31,16)`.  Fix a direction and
write

`eps S_H=P+sum_s n_s z_s+eps sum_(s<t) K_st z_s z_t`,       (1)

where:

- `I=sum_s n_s` is the infinity-edge count;
- `P` is the selected finite parallel-edge count;
- `K_st` is the signed selected-edge sum between fibres `s,t`.

If the target is `4+tau z_i z_j`, subtracting it from (1) gives a
multilinear quadratic which vanishes on `sum z=1`.  It therefore has the
form

`(sum z_s-1)(c+sum_s a_s z_s)`

after multilinearization.  Constant and linear comparison gives

`n_s=c-a_s`, `P-4=-c+sum_s a_s`,

hence

`I+P-4=30c`.                                             (2)

Integral pair coefficients and polarization imply that `2c` is integral,
so

`15 divides I+P-4`.                                      (3)

For a linear target `4+sigma z_j`, the identical comparison gives

`15 divides I+P-(4+sigma)`.                              (4)

The four equality targets needed below are:

| type/odd fibres | exact baseline `A` | target `eps S_H` | congruence |
|---|---|---|---|
| phase one, `b=2` | `(1-x_i-x_j)^2` | `4+z_i z_j` | `15 | I+P-4` |
| phase one, `b=30` | `1-x_j` | `4-z_j` | `15 | I+P-3` |
| phase zero, `b=2` | `(x_i-x_j)^2` | `4-z_i z_j` | `15 | I+P-4` |
| phase zero, `b=30` | `x_j` | `4+z_j` | `15 | I+P-5` |

Here `j` is the omitted fibre in the `b=30` rows.

## All fifteen hard baselines have b=2

The hard type has means `{30^15,62}`.  Same-type means satisfy

`a_d=I+32P_d-eps_d T-93`,                               (5)

so the fifteen mean-30 directions have one common parallel count `P`.
Proposition 15.728 gives at least fourteen `b=2` directions.  If the
remaining baseline had `b=30`, (3) and (4) would require the same `I+P` to
be both 4 and 3 modulo 15.  Therefore all fifteen have `b=2`.

## The opposite type has no selected finite edge

Set

`rho=(I+P-4)/15`, `s=rho+P`.                              (6)

The mean-62 hard direction has `P+1` parallel edges by (5).  Thus the number
of selected finite edges of the opposite Paley sign is

`E_opp=125-I-(15P+(P+1))=15(8-s)`.                       (7)

The integer `rho` is nonnegative.  Since infinity is outside the boundary,
`I` is even; (6) then makes `s` even.  Equation (7) leaves
`s in {0,2,4,6,8}`.

For an opposite-type direction with `Q_d` selected parallel edges, (5)
gives

`a_d=30s-208+32Q_d`, `sum_opp Q_d=15(8-s)`.              (8)

If `s<8`, nonnegativity forces `Q_d>=7-s`.  The excess above this common
minimum is

`15(8-s)-16(7-s)=8+s<16`,

so at least one direction has `Q_d=7-s` and mean
`a_d=16-2s in {16,12,8,4}`.  A positive phase-zero `b` has floor at least
32.  At `b=0`, parity gives `A=2C` for a nonzero nonnegative integral
quadratic, and Proposition 15.688 gives `4p E C>=28`.  Either way the small
mean is impossible.  Hence

`s=8`, `E_opp=0`, `I=124-16P`.                            (9)

The value `P=0` would give 124 infinity edges and one finite edge.  But if
`U` is the infinity-neighbour set and `F` the finite edge set, then
`D=U triangle partial(F)`, so `I<=|D|+2|F|=34`.  Thus `1<=P<=7`.

All selected finite edges now have the hard sign.  Every opposite direction
has `Q_d=0`, and (8) gives exact mean 32.

## Final contradiction

At phase zero and mean 32, the exact floor table permits only
`b in {0,2,30}`.  The `b=2` coefficient rule would require
`15 | I-4`; the `b=30` rule would require `15 | I-5`.  Equation (9) gives

`I == 4-P (mod 15)`, `1<=P<=7`,

so neither rule can hold.  All sixteen opposite directions therefore have
`b=0`.

The exact global odd-fibre identity is `sum_d b_d=72`.  The fifteen hard
baselines contribute `15*2=30`, while the opposite type contributes zero.
The one hard mean-62 direction would have to contribute `b=42`.  But `b` is
even and at most 30 among 31 fibres.  Contradiction.

Therefore the `p=31,R=10` residual endpoint does not exist.

## Optional stronger cell audit

The deterministic certificate also retains an unused strengthening.  From
the `b=2` coefficient cells,

`K_st=eps(g-n_s-n_t+1_(st=B_d))`, `g=(I+P-4)/15`,

and the transverse-edge norm bound, odd `g` is impossible: among 31 odd
numbers `2n_s-g`, at most 240 of 465 pairs can cancel, giving norm at least
224 against capacity at most 106.  The even case plus the opposite residue
reduces to

`(I,P)=(28,6),(60,4),(92,2)`.

Writing these as `r=1,2,3`, the infinity-star counts in each hard baseline
direction have profiles

- `r=1`: `0^3,1^28`;
- `r=2`: `0^1,2^30` or `1^2,2^29`;
- `r=3`: `2^1,3^30`;

and the exact transverse cells are
`N_st=d_s+d_t+1_(st=B_d)`, where `d_s=r-n_s>=0` and
`sum_s d_s=4-r`.  None of this extra structure is needed for the close.

## Artifacts

- `src/e1_gmin_m4_prop15733.py`
- `tests/test_prop15733.py`
- `evidence/e1_gmin_m4_prop15733.json`
