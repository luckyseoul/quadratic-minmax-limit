# P31 hard-fixed complement-pair centre obstruction

Classification: rigorous obstruction for one fixed sixteen-half family, not a
global exclusion and not closure of residual (ii).

The complement-pair/profile construction and its exact boundary matching give
the sixteen `(target, auxiliary, scale)` choices recorded in
`resii_p31_hard_fixed_center_obstruction_v1.json`.  They replay the exact raw
hard-fixed profile, with aggregate selector word `0x00800005` and correction
support `{2,23}` for the fixed/cancellation direction `0`.

For each of the 120 unordered half pairs, the checker exhausts all `30^2`
nonzero-centre pairs.  Of the resulting 108,000 cases, 107,700 are disjoint
and 300 share one inversion orbit.  Those 300 shared orbits have spatial
directions `4,9,20,21,29`, exactly 60 times each; direction `0` never occurs.
Every orbit repeated among two or more halves occurs in at least one pairwise
intersection.  Consequently no centre tuple for these fixed sixteen halves
can realize the cancellation in required direction `0`, whether the single
cancellation unit is a two-half opposite pair or a three-or-more-half
collision.

Replay:

```bash
python scripts/residual_branch_c_hard_fixed_center_obstruction.py \
  --output /tmp/resii_p31_hard_fixed_center_obstruction_v1.json
pytest -q tests/test_residual_branch_c_hard_fixed_center_obstruction.py
```

The next legitimate search must impose collision compatibility while choosing
the complement pairs and target labels; changing centres alone cannot repair
this fixed family.
