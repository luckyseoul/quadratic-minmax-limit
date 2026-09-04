# Proposition 15.769: the p=29 endpoint closes; the old p=23 shortcut does not

Date: 2026-09-04

Scope: the first layer beyond Proposition 15.752 at `p=29` and `p=23`.
This is not a closure of residual (ii).

## Dimension-free three-quarter cube lemma

Let `g` be a nonnegative integer-valued polynomial of degree at most two on
a Boolean cube and suppose `E[g]=3/4`. Then

```text
max(g) <= 6.
```

The result is sharp: with `s=x1+...+x6`,

```text
g=6-3s+binom(s,2)
```

has layer values `6,3,1,0,0,1,3` and mean `3/4`.

For the proof, assume a minimal-dimensional counterexample with integer
maximum `M>=7` at the origin. Every facet mean is in `(1/4)Z`. A facet
through the origin cannot have mean `0`, `1/4`, or `1/2`, by nonnegativity,
the degree-two support floor, and Proposition 15.751's half-mean
maximum-three theorem. It cannot have mean `3/4` by minimality. Its opposite
facet therefore has mean at most `1/2` and maximum at most three. Hence every
nonorigin value of `g` is at most three.

If there are at least five coordinates, average any five-coordinate subcube
by Hamming layer, writing `q(s)` for the layer mean. Quadratic interpolation
at `s=1,3,5` gives

```text
M=q(0)=(15/8)q(1)-(5/4)q(3)+(3/8)q(5) <= 27/4 < 7.
```

Dimensions at most three have total mass below seven. In dimension four,
the zero fourth difference gives

```text
M <= (layer-one mass)+(layer-three mass) <= 12-M,
```

again yielding `M<=6`.

## The p=29 local mass and shell

At `p=29`, suppose `4p E[B]=p+15=44` and put `H=max B`. If `H>=2`, paired
cubes first give `H>=4`. A half-mean cube has maximum at most three, so every
paired cube through the maximizer has mean at least `3/4`. This forces
`H>=12`. The exact stabilizer bound is `H<=88/7<13`, hence `H=12`.
The paired-cube average is then `23/30`, so one paired cube has mean exactly
`3/4`; the new lemma bounds its maximum by six although it contains the
value twelve. Contradiction.

If `H=1`, the corrected Johnson influence calculation gives largest
zero-class complement

```text
L <= 399168/54665 < 8.
```

Thus `L<=7<14`, every pattern extends to the complementary slice, and cube
influence reduces to four active coordinates. The fixed Proposition 15.751
catalog misses density `11/29`. Therefore the local mass 44 is impossible.

At the first post-band `p=29` layer, `t=11`, `k=138`, `|H|=139`. The two old
branches force the already excluded masses `p+9=38` and `p+7=36`. In the new
complement-triple branch, coefficient offset two and the common row sum
force hard parallel count `P=2`. The opposite parallel total is 109.
`Q=6` has mass 14 below the sharp lift floor 26, while after raising every
opposite direction to `Q=7` the surplus is only four. A `Q=7` direction is
therefore forced and has the excluded local mass 44. Hence

```text
p=29, k=138 is empty for every boundary size.
```

## Why the old p=23 offset-one shortcut fails

At `p=23`, the four-coordinate polynomial
`3-2r+binom(r,2)` has local mass `p+13=36` and signed-target offset one.
That is not the whole sharp half-mean equality mechanism. For a five-set
`R`, the same formula has layer values

```text
3,1,0,0,1,3,
```

slice mean `9/23`, paired-cube mean `1/2`, and maximum three. Its signed
target is

```text
3+4C = 5 + sum_({i,j} subset R) z_i z_j,
```

so its coefficient offset is five. Proposition 15.768's forced opposite
counts are `Q in {4,5,6,7}`; this example is compatible with `Q=5`, namely
the `P=4` complement-literal plus all-equal-triple hard family.

Thus the four-coordinate witness plus the offset congruence alone does not
close `p=23,k=110`.  This failed shortcut is now superseded by the complete
equality-globalization and quartic/octic moment argument in
`NOTE_2026-09-04_P23_POST_BAND_MOMENT_CLOSE.md`, which does close that one
endpoint.  The five-coordinate polynomial is the last local equality form,
not a residual graph construction; the new moment certificate excludes its
compatibility with the global coefficient graph.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_prop15769.py
```
