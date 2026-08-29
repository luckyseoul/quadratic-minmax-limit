# p=17 slack-sixteen orbit-extension exclusion

**Status:** proved conditional on Orbiter build 3361's exhaustive PGL orbit
classification, independently checked against Sticker's published count.

Proposition 15.704 left thirteen phase-labelled slack-sixteen profiles, all
without an undetermined direction. Every realization can be repaired to an
arc by at most four deletions; deleting extra retained points if necessary
gives a twelve-point arc core `A` and four restored points `D`.

If `mu_A(x)` is the number of core secants through `x`, the core secant lines
alone contribute at least `4 sum_{x in D} mu_A(x)` to pair slack. Thus a
slack-sixteen realization must satisfy `sum mu_A(x) <= 4`. This includes
shallower repairs and avoids an index-one equality assumption.

Orbiter enumerates 629 PGL classes of twelve-arcs in `PG(2,17)`. The local
independent collinearity audit splits them into 553 complete and 76 extendible
classes, exactly Sticker's published complete-class count. Across every class,
there are 97,122 four-point extensions satisfying the core-secant charge
bound. Direct occupancy checks leave only 47 extensions, representing ten
distinct projective point sets. All have the pattern `(n3,n4)=(4,0)`; the
other two arithmetically possible patterns never occur.

For every survivor, every projective line disjoint from the boundary was used
as the line at infinity. The 6,345 resulting affine charts produce 317
unlabelled Paley-phase profiles. Testing both global phase labellings yields
zero matches with the thirteen arithmetic targets. Therefore all thirteen
slack-sixteen rows are impossible. The p17 remainder drops from 654 to 641:
two slack-zero rows and 639 rows of slack at least twenty.

Reproduce with:

```bash
python scripts/p17_slack16_orbiter_extension.py
python src/e1_gmin_m4_prop15705.py
python -m pytest tests/test_prop15705.py -q
```

The archived Orbiter representative CSV is SHA-256
`0a57481731e10d55eb16a24158d57ca738240a9b32d3f66b9a39d85a64f16e24`.
