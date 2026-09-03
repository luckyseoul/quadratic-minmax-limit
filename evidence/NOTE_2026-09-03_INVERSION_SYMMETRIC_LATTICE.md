# The inversion-symmetric edge-Radon lattice

Date: 2026-09-03

Status: exact unrestricted integral cokernel and mod-two image proved. The
restricted central Boolean box, the coupled symmetric half, and residual
(ii) remain open.

## 1. Symmetric ranks

Let \(p=2h+1\) be odd, let \(d=p+1\), and let \(J(e)=-e\) act on source
edges. Let \(I\) be the corresponding central inversion on edge-Radon
target cells. Write

\[
 E^+=\ker_{\mathbf Z}(1-J),\qquad
 {\cal A}^+={\cal A}\cap\ker_{\mathbf Z}(1-I),
\]

where \({\cal A}\) is the ordinary compatible target lattice of
Proposition 15.760.

There are \(dh=(p^2-1)/2\) edge-difference classes and therefore
\(p^2dh\) source edges. Exactly \(dh\) source edges are fixed by \(J\):
the antipodal edges with midpoint zero. Hence

\[
 \operatorname{rank}E^+
 ={p^2dh+dh\over2}
 ={dh(p^2+1)\over2}.                                      \tag{1}
\]

The full compatible target has rational rank \(dph\). Its antisymmetric
part has rank \(dh^2\), so

\[
 \operatorname{rank}{\cal A}^+=dh(h+1).                   \tag{2}
\]

Equivalently, one direction has \(h+1\) inversion-fixed raw coordinates
and \(h^2\) nonfixed coordinate pairs; the \(d\) ordinary compatibility
equations lie in the symmetric part. Rational equivariant projection gives
surjectivity onto (2). Thus

\[
 \boxed{\operatorname{rank}\ker(R:E^+\to{\cal A}^+)
        =dph^2.}                                          \tag{3}
\]

The plus and minus kernels have the same rank, although their source and
target ranks differ.

## 2. Mod-two surjectivity

The proof separates target coordinates fixed and nonfixed under \(I\).

A fixed antipodal source edge \(\{v,-v\}\) maps only to fixed target
coordinates: \(P_L\) when \(L(v)=0\), and the cell \(\{s,-s\}\) otherwise.
There are \(dh\) such source variables. On the nonzero affine fibres this
is the affine Radon transform of an even point function vanishing at the
origin. If all these fibre sums vanish, Fourier inversion recovers every
nonzero frequency and forces the point function to vanish. The fixed-edge
map is therefore injective.

There are \(d(h+1)\) fixed target coordinates. Their compatible subspace
has the following \(d\) independent equations:

1. all \(d\) directional fixed-coordinate totals are equal, giving
   \(d-1\) equations; and
2. the sum of all parallel coordinates is that common total.

Its dimension is \(d(h+1)-d=dh\), exactly the fixed-edge image rank.

A nonfixed source orbit contributes the pair \(e+Je\). Modulo two it
vanishes on fixed target coordinates and gives equal values on each
nonfixed pair of target cells. This is exactly the characteristic-two
reduction of the antisymmetric block in
NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md. That map is surjective
onto the \(dh^2\) paired nonfixed target coordinates.

The two target coordinate supports are disjoint. Consequently

\[
 \boxed{R:E^+/2E^+\longrightarrow{\cal A}^+/2{\cal A}^+
        \text{ is surjective}.}                            \tag{4}
\]

The executable certificate builds the full plus-map over \(\mathbf F_2\)
for \(p=3,5,7,11\). It independently obtains fixed rank \(dh\), nonfixed
rank \(dh^2\), and total rank \(dh(h+1)\). These are fail-when-wrong
replays of the symbolic proof, not a prime census.

## 3. Exact integral cokernel

Put

\[
 C^+={\cal A}^+/R(E^+),\qquad G={\cal A}/R(E).
\]

The natural map \(C^+\to G\) is injective. If \(y\in{\cal A}^+\) and
\(y=Rz\), then

\[
                 2y=R(z+Jz),\qquad z+Jz\in E^+.
\]

Thus a kernel class is killed by two. Surjectivity (4) says that the finite
group \(C^+\) has no two-primary quotient; multiplication by two is an
automorphism, so the kernel is zero.

The image is exactly the \(+1\) eigenspace \(G^+\). Containment is clear.
Conversely, if \(Ig=g\), then \((1+I)y\) represents \(2g\). Proposition
15.760 proves that \(G\) is an elementary \(p\)-group, and two is invertible
modulo the odd prime \(p\). Hence \(g\) lies in the image.

Degree-\(q\) moment rows have inversion eigenvalue \((-1)^q\). Therefore
\(G^+\) is exactly the even-degree part of Proposition 15.759's moment
basis. Its rank is

\[
 S_+(p)
 ={(h-1)(2h^2+5h+6)\over6}.                               \tag{5}
\]

It follows that

\[
 \boxed{{\cal A}^+/R(E^+)
 \cong(\mathbf Z/p\mathbf Z)^{S_+(p)}.}                    \tag{6}
\]

An integral symmetric compatible target has an unrestricted integral
symmetric source if and only if all its even-degree moment congruences
vanish. There is no extra Smith, parity, \(p^2\), or other-prime
obstruction in the symmetric block.

## 4. Exact box after the Möbius antisymmetric lift

Use ordinary row-untransported coordinates, so the desired target is
\(Y=R(\tau x)\). Let \(x_U\) be the physical graph selecting one edge from
each orbit used by the direction-localized Möbius trades, and put

\[
 q_U=\tau x_U,\qquad z=q_U-Jq_U.
\]

The antisymmetric theorem gives

\[
                  Rz=Y-IY.
\]

The forced symmetric source on the used orbits is

\[
 C_U=q_U+Jq_U
    =\sum_{O\ {\rm used}}\tau_O(e_O+Je_O).
\]

After subtracting the actual source \(q_U\), the exact remaining target is

\[
 T_U=Y-Rq_U={Y+IY-RC_U\over2}.                             \tag{7}
\]

This is integral because it is the difference of two integral targets.
It is central by construction.

Every used nonfixed orbit is now frozen at zero: selecting its other edge
would destroy the already fixed antisymmetric difference. On an unused
nonfixed orbit the remaining physical graph may select neither edge or both,
giving ordinary source coordinate \(0\) or
\(\tau_O(e+Je)\). A fixed antipodal edge independently contributes \(0\)
or \(\tau_f f\). Thus the live condition is exactly

\[
 \boxed{
 (z_0+\ker_{\mathbf Z}R|_{E^+})
 \cap
 \left(
  \prod_{O\ {\rm used}}\{0\}
  \times\prod_{O\ {\rm unused}}\{0,\tau_O\}
  \times\prod_{f\ {\rm fixed}}\{0,\tau_f\}
 \right)\ne\varnothing .}                                 \tag{8}
\]

Here a nonfixed quotient coordinate \(\tau_O\) means the source vector
\(\tau_O(e+Je)\).

Equation (6) supplies an unrestricted signed integral central lift of
\(T_U\), because subtracting an actual source preserves every moment
compatibility. Equation (4) independently supplies a mod-two central lift.
They do not imply that one lift satisfies the coefficient restrictions in
(8). That simultaneous affine-box intersection is the exact remaining
symmetric obstruction.

No simple common graph is constructed, and residual (ii), E1, and
\(L=1/2\) remain open.

## Reproduction

    PYTHONPATH=src python src/e1_gmin_m4_inversion_symmetric_lattice.py
    PYTHONPATH=src pytest -q tests/test_inversion_symmetric_lattice.py
