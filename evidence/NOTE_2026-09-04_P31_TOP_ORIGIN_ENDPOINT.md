# The p31 top localized-Mobius endpoint has sixteen origin edges

Date: 2026-09-04

**Status:** proved a necessary endpoint condition for the all-active
localized-Mobius construction at `p=31,t=177`.  Parallel parity forces all
sixteen auxiliary directions to be distinct.  Consequently the sole
cancellation unit is nonorigin, the unique fixed antipodal edge is also
nonorigin, and all sixteen raw origin edges survive.  Their signed degree is
`+10` or `+12` in the `hard-minus-opposite` convention.  This does not
construct a common graph and does not close residual (ii).

## 1. The top parity vector

There are `m=16` localized halves and

\[
 m(p-1)=16\cdot30=480
\]

raw half-edge occurrences.  The top graph has `479` edges.  At the `j=0`
Hamming endpoint the fixed-edge theorem leaves exactly one fixed antipodal
edge, no unused doubled orbit, and `478` used nonfixed orbits.  Thus the raw
halves lose exactly one cancellation unit.

Retain the notation

\[
 v=P\pmod2+{\bf1}_{\rm hard},
 \qquad c=\text{parity vector of the auxiliary directions}.
\]

The hard parallel quotas are `14^14,15^2`; hence `v` is supported on the
fourteen hard directions of quota `14`.  The opposite quotas are
`15^3,16^13`; hence `v` is supported on the three opposite directions of
quota `15`.  Therefore

\[
                         \operatorname{wt}(v)=17.          \tag{1}
\]

There is one fixed direction `F` and no doubled column, so the exact
parallel parity equation is

\[
                              c=v+e_F.                     \tag{2}
\]

The parity vector of sixteen auxiliary occurrences has weight at most
sixteen.  If `F` were outside `supp(v)`, (2) would have weight eighteen.
Thus `F` lies in `supp(v)` and

\[
                              \operatorname{wt}(c)=16.     \tag{3}
\]

A parity vector of weight sixteen realized by exactly sixteen occurrences
has no repeated direction: every supported direction occurs exactly once.
Thus the sixteen projective auxiliary directions `M_i` are distinct and
equal `supp(v)\{F}`.

## 2. A common-origin collision forces equal auxiliary directions

For a localized half aimed at a nonzero center, its `t=0` edge is

\[
                              \{0,u_i\},\qquad
                              0\ne u_i\in\ker M_i.         \tag{4}
\]

If two edges in (4) represent the same central-inversion orbit, then
`u_i=+/-u_k`.  This forces `ker M_i=ker M_k`, hence equality of the
projective auxiliary directions.  The converse is neither needed nor
claimed: the magnitude of `u_i` along a fixed kernel line also depends on
the target, center, and chosen auxiliary scale.  Distinctness from (3)
therefore makes all sixteen actual origin-edge orbits distinct.  No
cancellation can occur there.  The sole cancellation unit forced by the
`480 -> 478` support drop is nonorigin, and every one of the sixteen origin
edges survives in `H`.

The unique fixed source edge is antipodal, `\{x,-x\}` with `x!=0`; over the
odd field it cannot contain the origin.  Consequently

\[
                                  \deg_H(0)=16.            \tag{5}
\]

## 3. The two signed-origin cases

The Paley sign of (4) is the sign of `M_i`.  Indeed a kernel representative
of `M_i=(a,b)` is `(b,-a)`, and the two anisotropic norms are both
`a^2+b^2`.

If `F` is hard, it removes one member of the fourteen-hard/three-opposite
support of `v`.  The auxiliary set then has thirteen hard and three opposite
directions.  If `F` is opposite, it has fourteen hard and two opposite
directions.  Hence

\[
 (\deg_+(0),\deg_-(0),\deg_+(0)-\deg_-(0))=
 \begin{cases}
 (13,3,10),&F\text{ hard},\\
 (14,2,12),&F\text{ opposite}.
 \end{cases}                                              \tag{6}
\]

Equivalently, for `w=deg_- - deg_+`, the forced values are `w(0)=-10` and
`w(0)=-12`.  Any top localized-half candidate with an origin cancellation,
origin degree fourteen, or a different signed origin value violates the
exact top parallel parity ledger.

## 4. Scope and replay

The theorem couples the one-cancellation support endpoint to the fixed edge
and to an actual source vertex.  It does not constrain the remaining 462
nonorigin edges enough to prove or disprove the transverse target.

Replay with

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_top_origin_endpoint.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_p31_top_origin_endpoint.py
```

Residual (ii), E1, and `L=1/2` remain **OPEN**.
