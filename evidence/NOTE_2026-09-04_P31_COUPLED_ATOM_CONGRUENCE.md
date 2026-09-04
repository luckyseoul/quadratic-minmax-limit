# Coupled atom congruences survive at the p31 hard-fixed endpoint

Date: 2026-09-04

**Status:** no mod-2, mod-3, or other pure signed-lattice congruence cut
separates the full hard-fixed ledger.  There is an explicit shared physical
479-edge ledger graph with the exact top parallel profile, one hard fixed
edge, no doubled nonfixed inversion orbit, and one common center tuple for
which all 32 normalized row residuals have signed integer atom decompositions
with the required atom-count coordinates.  These decompositions have
negative compact multiplicities, so this is not a nonnegative atom solution
and not a residual-(ii) witness.  A Chvatal--Gomory cut using atom nonnegativity remains
possible.

## 1. The exact compact-atom lattice

On the complete label graph, write

\[
 K(a,b;c)=e_{ab}-e_{ac}-e_{bc}.
\]

Every compact atom has even signed degree: its degrees at `a,b` are zero and
its degree at `c` is `-2`.  Conversely, let `x` be any integer edge vector
whose signed degree is even at every label.

Modulo two, `K(a,b;c)` is the triangle on `a,b,c`.  Triangles through label
zero generate the cycle space of the complete graph, so a sum of compact
atoms can be subtracted from `x` to leave an edgewise-even vector.  The exact
identity

\[
 K(u,w;v)+K(v,w;u)=-2e_{uv}                         \tag{1}
\]

then corrects every remaining cell independently.  Therefore

\[
 \mathbb Z\langle K(a,b;c)\rangle
 =\{x\in\mathbb Z^{\binom n2}:\partial x\equiv0\pmod2\}. \tag{2}
\]

The compact multiplicity coordinate needs no further congruence: because
each compact atom has edge sum `-1`, every decomposition in (2) has total
signed multiplicity `-sum(x)`.

Modulo three, the relabelled form of (1) gives

\[
 K(a,c;b)+K(b,c;a)=-2e_{ab}\equiv e_{ab}\pmod3,     \tag{3}
\]

so every edge coordinate is generated.  Again the atom-count congruence is
exactly the already-forced edge-sum identity.  In an opposite row, reserve
any six distinct all-positive triangles first.  Their total edge sum is 18
and their degrees are even, so (2) applies to the remaining compact target
with exactly the required compact count.

## 2. One shared graph satisfying every row parity

Freeze all sixteen nonzero hard centers at `j_N=1`.  The physical graph
boundary required by all rows at once is

\[
 B(v)=\bigoplus_{N\in H}1_{N(v)=1},                  \tag{4}
\]

For a hard projection `M`, each fibre meets the other fifteen hard lines
once, while the own center line contributes on the center fibre.  Thus its
projected boundary is `F_31 minus {1}`, exactly the boundary of the literal
star.  Every opposite fibre meets all sixteen hard lines once, so its
projected boundary is zero.  Hence (4) synchronizes all 32 mod-two row
conditions through one physical boundary; the rows are not treated as
independent variables.

The executable constructor realizes (4) together with the exact graph
ledger:

1. insert the fixed antipodal edge `{(1,30),(30,1)}` in hard spatial
   direction 1;
2. pair the remaining boundary vertices with clean nonfixed edges;
3. add a clean collinear triangle in each direction whose edge-count parity
   needs toggling;
4. within each direction, replace one edge by a clean length-three path as
   often as needed.  This preserves boundary and direction while adding two
   edges.  A four-cycle would seed an empty even direction, although the
   pinned p31 replay does not need a broad search.

The result has exactly the frozen profile

```text
(15,14,14,15,16,15,16,14,16,14,14,16,16,15,15,14,
 14,16,16,16,16,14,14,16,14,16,16,16,14,14,14,14)
```

and the following exact diagnostics:

```text
physical edges:                         479
fixed antipodal edges:                    1
doubled nonfixed inversion orbits:        0
vertex-boundary weight:                 452
hard projected boundary checks:       16/16
opposite projected boundary checks:   16/16
graph SHA256:
  36aea8d59a4131042de02a999a1f36070cc9d69150c19f17771543d45e46d116
```

## 3. Exact signed replay and scope

For every row, the code forms the actual Paley-signed projection of this one
graph, normalizes by the row sign, and adds the literal star in hard rows.
It reserves six distinct positive triangles in each opposite row and applies
the constructive proof of (2) to the remainder.  All 32 decompositions replay
coefficient by coefficient and have the exact compact and positive atom
counts.

```text
rows replayed:                          32
total nonzero signed compact entries: 16048
total signed compact l1:               25451
minimum compact coefficient:             -17
maximum compact coefficient:               5
all-row signed-decomposition SHA256:
  42745425e873264598bbec521a83e744db1db81ef0d80a39a32bfcecb144bee7
JSON replay SHA256:
  02ef46093a120ff658da0af2936c54db61457e00ad19a6c7fc6e8df45a3b266a
```

The negative coefficient `-17` is the essential limitation.  This result
rules out a closure based only on mod-2/mod-3 (indeed, any signed-lattice)
congruences, including such congruences coupled through the actual shared
graph variables.  It does **not** place this graph in the nonnegative atom
semigroups, does not refute nonnegativity-sensitive Chvatal--Gomory cuts, and
does not construct a common residual graph.  The next integral target is a
nonnegative semigroup inequality or synchronized atom obstruction, not
another modular census.

## Replay

```bash
PYTHONPATH=src python src/e1_gmin_m4_p31_coupled_atom_congruence.py \
  > /tmp/resii_p31_coupled_atom_congruence.json
sha256sum /tmp/resii_p31_coupled_atom_congruence.json
PYTHONPATH=src pytest -q tests/test_p31_coupled_atom_congruence.py
```
