# Whole-slab exchanges in one unused symmetric difference slice

**Status:** proved exact one-slice kernel and bounded-fibre connectivity after
fixed-edge elimination. This does not prove normality of the full unused
configuration, nonemptiness of the restricted symmetric Boolean fibre, or
closure of residual (ii).

## 1. Input after the fixed antipodal choices

In fixed-cell/nonfixed-pair coordinates, the symmetric Radon map has block
form

\[
 R^+=\begin{pmatrix}A&2B\\0&C\end{pmatrix}.
\]

Modulo two, $A$ is an isomorphism onto the compatible fixed-cell data, so
the fixed source choice is unique. After subtracting that choice and dividing
the even fixed-cell residual by two, the remaining problem is

\[
 \sum_{O\notin U} b_O\widehat B_O=\widehat T,qquad
 b_O\in\{0,1\},qquad
 2\sum_Ob_O=|H|-|U|-|a|.                                  \tag{1}
\]

Here $U$ is the set of nonfixed inversion orbits frozen by the chosen
antisymmetric Mobius lift. This note determines the kernel of (1) supported
on one difference class.

Write source edges in midpoint coordinates

\[
 e=(a,[\delta])=\{a-\delta,a+\delta\},
 \qquad [a],[\delta]\in
 \Delta=(\mathbf F_p^2\setminus\{0\})/\{\pm1\}.
\]

Thus the remaining nonfixed pair columns are indexed by
$\Delta\times\Delta$. The theorem below is for the unsigned reduced columns
and holds for every odd prime. In the branch-C specialization
$p\equiv3\pmod4$, the anisotropic Paley source sign is nonzero, depends only
on $[\delta]$, and is therefore constant throughout a fixed difference
slice; multiplying the slice by that sign does not change its kernel. No
Paley-sign claim is made for the unsigned $p=5$ replay.

## 2. Exact one-slice kernel

Fix $[\delta]$, and let $L_\delta$ be a nonzero functional with
$L_\delta(\delta)=0$. For a nonzero square class
$s\in\mathbf F_p^\times/\{\pm1\}$, put

\[
 \mathcal S_{\delta,s}
 =\{([a],[\delta]):L_\delta(a)^2=s\}.                       \tag{2}
\]

There are $h=(p-1)/2$ such slabs, each containing exactly $p$ nonfixed
midpoint orbits. The zero-label slab has $h$ orbits.

**Theorem.** Let $v\in\ker_\mathbf Z\widehat B$ be supported on the one
difference slice $[\delta]$. Then there are integers $\gamma_s$, indexed
by the $h$ nonzero square classes, such that

\[
 v=\sum_{s\ne0}\gamma_s\mathbf1_{\mathcal S_{\delta,s}},
 \qquad \sum_{s\ne0}\gamma_s=0.                            \tag{3}
\]

In particular $v$ vanishes on the zero-label slab. Conversely, every
vector in (3) lies in the kernel. Hence the one-slice kernel is the root
lattice $A_{h-1}$, and its primitive circuits are exactly

\[
 \boxed{
 \mathbf1_{\mathcal S_{\delta,s}}
 -\mathbf1_{\mathcal S_{\delta,t}},\qquad s\ne t.}          \tag{4}
\]

Each circuit exchanges $p$ pair variables for $p$ pair variables. It
therefore preserves the constant weight in (1), changes $2p$ pair
coordinates, and removes and adds $2p$ physical graph edges. Because (4) is
an actual $\widehat B$-kernel relation, it preserves every reduced target
coordinate, including the complete vector of parallel weights in all $p+1$
directions; it is not merely a scalar constant-weight move.

### Proof

Lift $v$ to an even function $f:\mathbf F_p^2\to\mathbf Z$ by
$f(a)=v_{[a],[\delta]}$ for $a\ne0$ and $f(0)=0$. In a direction
$M\ne L_\delta$, the difference is transverse. Vanishing of every reduced
central-cell coordinate says exactly that every affine line sum

\[
 \sum_{a:M(a)=c}f(a)                                      \tag{5}
\]

is zero. For $c\ne0$, a midpoint orbit has one representative in each of
the fibres $c,-c$; for $c=0$, the reduced coordinate is half the even
line sum, so the conclusion is the same.

The finite Fourier slice identity now makes every Fourier coefficient of
$f$ outside the one-dimensional dual space spanned by $L_\delta$ vanish.
Therefore

\[
                         f(a)=g(L_\delta(a))                 \tag{6}
\]

for some $g:\mathbf F_p\to\mathbf Z$. Since $f(0)=0$, one has
$g(0)=0$. Since $f$ is even, $g(-x)=g(x)$. Finally the parallel row, or
equivalently any transverse line sum in (5), gives

\[
                  0=\sum_{x\in\mathbf F_p}g(x)
                    =2\sum_{s\ne0}\gamma_s.                \tag{7}
\]

This proves (3). The converse is the Type-P ridge identity: a zero-sum
profile constant on $L_\delta$-fibres has zero Radon image in every other
direction and zero total in the parallel row.

The primitive support-minimal nonzero vectors of
$\{\gamma\in\mathbf Z^h:\sum\gamma_s=0\}$ are
$e_s-e_t$, proving (4).

## 3. Deleting the Mobius-used orbits

Extend a move on the unused variables by zero on $U$. Equation (3) shows a
strong rigidity: if $U$ meets one point of a slab
$\mathcal S_{\delta,s}$, then that slab's coefficient must be zero. If
$q_\delta$ nonzero slabs are disjoint from $U$, the unused one-slice
kernel is therefore exactly $A_{q_\delta-1}$ on those clean slabs.

It follows that any two binary solutions of (1) which agree outside the
$[\delta]$-slice are connected inside the binary fibre by whole-slab
exchanges (4). Indeed their difference has coefficients
$\gamma_s\in\{-1,0,1\}$; pair each $+1$ slab with a $-1$ slab and swap
them one pair at a time.

There are

\[
                         |\Delta|={p^2-1\over2}              \tag{8}
\]

difference slices. At most $(p+1)(p-1)/2=|\Delta|$ orbits are used by the
disjoint Mobius construction. Hence some difference slice contains at most
one used orbit. It retains at least $h-1$ clean nonzero slabs and therefore
at least

\[
                         \binom{h-1}{2}                      \tag{9}
\]

unused circuits. For $p=31$, this gives one slice with at least fourteen
clean slabs and ninety-one whole-slab circuits.

## 4. Exact boundary of the result

The theorem classifies only moves whose support lies in one difference
slice. A general kernel move can couple several difference classes, and an
arbitrary binary point need not contain an entirely selected clean slab on
which (4) can act. Thus (4) does not prove global exchange connectivity,
normality of the unused-column semigroup, or existence of a Boolean point in
(1). It also does not characterize the image of the punctured map
$\widehat B\bmod2$, which is the next linear obstruction before any global
normality claim. It does prove that fixed-edge elimination leaves a large,
explicit, weight-preserving circuit family and that no smaller partial-slab
repair is possible inside one difference class.

The symbolic record and exact unsigned $p=5,7$ rank replay are in
`src/e1_gmin_m4_symmetric_slice_exchange.py`; focused tests are in
`tests/test_symmetric_slice_exchange.py`.

## Reproduction

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
      tests/test_symmetric_slice_exchange.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_slice_exchange.py
