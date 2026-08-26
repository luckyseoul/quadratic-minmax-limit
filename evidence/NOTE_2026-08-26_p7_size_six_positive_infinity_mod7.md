# p=7 positive infinity-plus-five boundary exclusion

For a six-point boundary containing infinity at `p=7`, the five finite
points give an odd number `b=1,3,5` of odd fibres in every affine
direction.  In the `c_H=+1` branch the parity phase is zero, and every one
of these three fibre sizes has scaled floor eight.  Four directions of
each quadratic type therefore saturate the exact type budget 32, so every
direction has scaled mean eight.

The complete rank-21 evaluation on `J(7,4)` is unique at that mean:
`A(X)=|X cap B| mod 2` for `b=1,5`, and
`A(X)=(|X cap B|-2)^2` for `b=3`.  Consequently every finite five-set fixes
all 280 affine bad-edge counts of a putative 29-edge graph.  Together with
the edge count and distinguished edge, these give the same `282 x 1225`
integer system used in Proposition 15.655.  Its rank over `F_7` is 147,
leaving 135 exact left-null dependencies.

The V100 integer kernel checked all `C(49,5)=1,906,884` finite boundaries
and found zero dependency-compatible right sides in 2.83 seconds.  A
separate NumPy run on NUKA repeated the complete enumeration in 4.47
seconds.  Both runs report the same direction-mask histogram
`{1:2923536, 3:9507960, 5:2823576}` and zero survivors.  The raw records and
the exact implementation are archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-p7-size-six-positive-infinity/`.

This closes only `p=7`, `c_H=+1`, infinity plus five finite points.  The
negative-product infinity branch, six finite points at `p=7`, all `p=5`
size-six cases, larger boundaries, residual (ii), R1, and the limit remain
open.
