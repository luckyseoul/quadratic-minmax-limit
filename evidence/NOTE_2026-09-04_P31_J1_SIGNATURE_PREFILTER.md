# Bounded j=1 adaptive-signature prefilter at the p=31 top endpoint

Date: 2026-09-04

**Status:** this is an exhaustive replay over the frozen 10,000-vertex BFS
checkpoint with SHA256
`196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b`.
Every stored design is excluded in the `f=1,d=1` j=1 correction branch by an
exact adaptive-kernel necessary condition.  The BFS checkpoint still has
9,018 frontier vertices and is not an exhaustive enumeration of its
component.  This result is therefore **bounded component evidence**, not a
global certificate.  It does not address the `f=3,d=0` j=1 branch and does
not close residual (ii).

## Exact necessary condition

For a stored exact-profile design, let `G` be the XOR of the sixteen
32-direction half-signature words from the closed kernel formula in
`e1_gmin_m4_p31_top_mobius_boundary_parity.py`.  In the `f=1,d=1` ledger the
correction consists of:

- one fixed antipodal edge in the forced direction `F=5`, with signature
  `e_F`;
- one unused doubled origin orbit in some direction `N`, with signature
  `e_N`; and
- two nonorigin inversion-orbit cancellations, each having a signature word
  of Hamming weight at most two.

Consequently any physical completion must satisfy

```
G + e_F + e_N = v_collision_1 + v_collision_2,
```

and in particular

```
min_N wt(G + e_F + e_N) <= 4.                         (1)
```

This is only a necessary prefilter.  Passing (1) would not establish physical
orbit or centre realizability.

## Pinned inputs and independent option replay

The scan uses no centre search and generates no new BFS vertices.  Its two
inputs are pinned byte-for-byte:

- bounded auxiliary checkpoint:
  `196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b`;
- serialized closed option catalog:
  `6c2c9dd2ca12d007f865c4499dbe038ef53215688eb43bb3759aef7d39daa599`.

The replay independently rebuilds all 7,260 allowed
`(target,auxiliary,scale,P2-mask,g-mask)` records from the exact parallel and
kernel-signature formulae.  Their serialization-independent record-stream
SHA256 is
`aa33d94449b90b8770ec819044552f8a575465203fc0ed29d9e69dcef37764ff`.
Every one of the 10,000 design keys is also checked against the exact
32-direction `P=2` cover ledger before its signature is used.

## Result

The exact histogram of the left side of (1), minimized over all 32 choices of
`N`, is:

| minimum transformed distance | number of checkpoint designs |
|---:|---:|
| 6 | 18 |
| 8 | 99 |
| 10 | 624 |
| 12 | 1,549 |
| 14 | 2,821 |
| 16 | 2,482 |
| 18 | 1,669 |
| 20 | 616 |
| 22 | 116 |
| 24 | 6 |

Thus all 10,000 stored vertices have transformed distance at least six, while
two cancellation words can supply weight at most four.  None passes the
necessary condition.

An intermediate conjecture that every exact-profile design has `wt(G)>=8`
is false.  Exactly three stored designs already have `wt(G)=6`:

| component id | `G` | support | `min_N wt(G+e_F+e_N)` |
|---:|---:|---|---:|
| 5497 | `0010a052` | `1,4,6,13,15,20` | 6 |
| 8971 | `18200083` | `0,1,7,21,27,28` | 6 |
| 9077 | `20d40400` | `10,18,20,22,23,29` | 6 |

Their full 16-option design keys are preserved in the evidence JSON.  These
three counterexamples retract the proposed raw-weight lower bound, while the
correct transformed-distance prefilter still excludes them.

## Replay

From the repository root, with the two pinned files at the displayed paths:

```bash
python scripts/residual_branch_c_j1_signature_prefilter.py \
  /tmp/resii_auxiliary_component_10000_v1.json \
  /tmp/resii_p31_j1_closed_option_catalog_v1.json \
  --output /tmp/resii_p31_j1_signature_prefilter_replay.json
sha256sum \
  /tmp/resii_auxiliary_component_10000_v1.json \
  /tmp/resii_p31_j1_closed_option_catalog_v1.json
```

The checked-in result is
`evidence/resii_p31_j1_signature_prefilter_v1.json`.  Its classification and
Boolean scope fields are intentionally part of the certificate: the bounded
checkpoint is fully scanned, the underlying BFS component is not exhausted,
the alternate j=1 branch is untouched, and residual (ii) remains open.
