# p=19 low-slack conic reduction

**Date:** 2026-08-29
**Proposition:** 15.689
**Status:** proved conditional on the published complete-arc spectrum;
reduces but does not close the endpoint

Proposition 15.688's exact residue-zero census has 143 phase-labelled
profiles (75 global shapes):

```text
slack   0   4   8  12  16  20  24  28  32
count  54  37  25  13   7   4   1   1   1
```

The complete-arc spectrum of `PG(2,19)` has complete sizes
`10,11,12,13,14,20`; the size-20 arc is the conic. Hence every arc of size
at least 15 is conic-contained.

- Slack zero gives a 16-arc. Three undetermined infinity points force two
  overlapping conic extensions and then three collinear conic points. With
  one or two undetermined directions, the conic is tangent or secant to the
  line at infinity. Every other direction has at least six retained affine
  secants and therefore `b<=4`; every exact profile has a non-undetermined
  `b>=6`.
- Slack four repairs by one deletion to a 15-arc. The deleted off-conic
  point has at least four retained conic secants, forcing slack at least 16.
- Slack eight or twelve repairs by at most two or three deletions. Every
  profile supplies two undetermined infinity points, producing a conic arc
  of size at least 16 or 15. If `j` deleted points are off the conic, the
  line-slack inequality gives `slack >= 4j(5-j)`, namely `16,24,24` for
  `j=1,2,3`.

Thus all `54+37+25+13=129` low-slack profiles are impossible. Exactly
fourteen arithmetic profiles remain:

```text
{16:7, 20:4, 24:1, 28:1, 32:1}
```

External input: G. Faina, S. Marcugini, A. Milani, and F. Pambianco,
*The spectrum of values k for complete k-arcs in PG(2,q) for q<=23*, Ars
Combinatoria **47** (1997), 3--11. H. Sticker's complete classification
table independently gives the same `PG(2,19)` spectrum.

Reproduction:

- `src/e1_gmin_m4_prop15689.py`
- `evidence/e1_gmin_m4_prop15689.json`
- `tests/test_prop15689.py`
