# A joint hard-fixed witness exists, but its full boundary system is inconsistent

Date: 2026-09-04

Scope: the displayed `p=31` localized-Mobius sixteen-half family in the E1
common-graph reduction. This is not a residual-(ii) closure.

## Exact positive result

The deterministic certificate
`src/e1_gmin_m4_p31_hard_fixed_joint_witness.py` replays sixteen labelled
halves which simultaneously satisfy all of the following:

- their targets are the sixteen hard directions and their auxiliary
  directions form an SDR;
- their exact raw parallel profile is
  `(15,15,14,14,15,16,16,14,16,14,14,16,15,16,16,14,14,16,16,16,16,14,14,16,14,16,15,16,14,14,15,14)`;
- their selector signature is `0x01000401`, so removing fixed direction zero
  leaves precisely the required endpoint correction `{10,24}`;
- halves 1 and 5 admit thirty scalar-equivalent opposite physical collisions
  in direction zero with endpoint support `{10,24}`; no other half pair has
  that prescribed collision key;
- one explicit centre tuple has exactly that one cancellation and, after
  adding the fixed antipodal edge, gives a simple 479-edge graph with SHA-256
  `6e924c3acb493799f7951a6ab75e22a2628d452d4654643944c1f3871a75a6a4`.

Thus parallel profile, selector signature, and physical collision are jointly
feasible. They cannot be separated and then cited as a contradiction.

## Exact rejection of this witness

For each of the thirty collision lifts, freeze the two collision centres and
leave all 30 choices for each of the other fourteen centres, together with all
15 fixed antipodal edges. The necessary vertex-boundary equations over GF(2)
have 435 variables and 976 equations. Every scalar lift gives coefficient rank
225 and augmented rank 226, hence is inconsistent. The canonical lift's
augmented column matrix hash is
`c32f2cecfe8c378f290fc3cfae70e2fa2e85046bae7981572c6537c1147f4207`.

An independent row-oriented elimination of the canonical lift produces the
same ranks and an explicit XOR contradiction using 104 vertex equations. The
inconsistency is stronger than failure of an exactly-one search: the target
boundary is outside the span of all allowed centre and fixed-edge columns.

Consequently this exact sixteen-half witness cannot be the required common
graph. Other labelled hard-fixed half families, nonzero-form or unbalanced
families, arbitrary lifts, E1, and residual (ii) remain open.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_hard_fixed_joint_witness.py
```
