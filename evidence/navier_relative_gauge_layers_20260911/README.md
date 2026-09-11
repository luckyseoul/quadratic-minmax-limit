# Exact relative-gauge connected-layer probe

**Status: exact finite obstruction to truncated correction states; no
all-orders composition conclusion.**

The Navier-style correction-cycle transfer identifies the uncontrolled
connected Eulerian terms in the relative-gauge conditional density as the
actual residual.  This probe calculates their first layers exactly, rather
than assuming the mixed four-cycle Hamiltonian `H4` is a closed state.

For each random pair of symmetric order-5 signings and an order-5-by-5 cross
signing, `scripts/navier_relative_gauge_layer_probe.py` evaluates all 512
projective row/column/sign gauges.  It enumerates the Boolean partition
function in rational arithmetic and obtains the coefficients of `t^4`,
`t^6`, and `t^8` in `log E cosh(t Q_Y)`.  The three block-local normalizers
in the conditional density are gauge invariant, so unequal coefficients
remain unequal in the conditional log density.

The smallest deterministic witness is already at a 4+4 block (`--seed 1
--order 4`): two gauges have equal `H4=-37/3` but respectively have
`H6=4228/45` and `H6=628/45`.  Thus `H4` alone is not a closed
relative-switching state.

Soulkiller independently ran seeds 1 through 88 at a 5+5 block.  Every
instance contained both:

* a pair of gauges with equal `H4` and different `H6`; and
* a pair with equal `(H4,H6)` and different `H8`.

```text
instances                           88
H4/H6 collisions                    88
(H4,H6)/H8 collisions               88
concatenated receipt SHA-256        a7baeec66aca596833c09b68da09271020e6bea57f1f0665dda66c3b91fa8132
```

This is a precise correction-cycle boundary: a proof based on a fixed
four-cycle residual, or even on the pair `(H4,H6)`, cannot recover the exact
finite conditional density.  It does **not** prove that no finite or compact
state exists, nor that the higher layers cannot be uniformly controlled by a
new all-orders argument.  The remaining valid target is a state including or
bounding the complete connected Eulerian tail with a Dini-summable error.

