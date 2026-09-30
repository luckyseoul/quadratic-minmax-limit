# Square-affine orbits at +169 on the order-50 plus-I lift

Let \(B\) be the quadratic-residue Seidel matrix of \(\mathbf F_{25}\) and
\(K=K(B)\) its plus-I lift. The square-affine group
\(x\mapsto ax+b\), with \(a\) a nonzero square, has order 300 and acts by
the same permutation on both blocks. That action preserves \(K\).

Coordinate ascent on \(Q\) finds six orbits of states with \(Q=+169\).
They form three pairs. In each pair the two orbits are global sign flips
of one another, and no square-affine map sends a state to its negative.

| size | orbits | local fields \(z_i(Kz)_i\) |
|------|--------|------------------------------|
| 150 | 2 | \(1\) (\(\times2\)), \(5\) (\(\times28\)), \(9\) (\(\times16\)), \(13\) (\(\times4\)) |
| 150 | 2 | \(5\) (\(\times28\)), \(9\) (\(\times22\)) |
| 50 | 2 | \(1\) (\(\times2\)), \(5\) (\(\times24\)), \(9\) (\(\times24\)) |

Each multiset sums to \(338=2\cdot169\). These orbits contain 700 states
at \(+169\). The block map \((x,y)\mapsto(y,-x)\) sends \(Q\) to \(-Q\)
and does not meet them. Multiplication by a nonsquare does not preserve
\(K\); it carries these states to \(-151\) or \(-167\).

Every one of the three local-field patterns splits all 1225 edges. This
does not change the one-edge or two-edge census, and it does not decide
\(m_{50}\) or the limit.

## States fixed by an involution

Every involution \(x\mapsto -x+b\) in this group is a translate of
\(x\mapsto -x\), and each has one fixed point. A \(\pm1\) state fixed by
\(x\mapsto -x\) takes one sign on each pair \(\{t,-t\}\), in both blocks:
\(2^{26}\) states.

Enumerating that set, the maximum of \(Q\) is 169, attained at exactly 28
states. The same count at \(-169\) is 28. The six orbits above contribute
exactly these states: 6 per size-150 orbit and 2 per size-50 orbit. A
state at \(+169\) outside the six orbits cannot be fixed by any of these
involutions, so its square-affine stabilizer is trivial and its orbit has
size 300.
