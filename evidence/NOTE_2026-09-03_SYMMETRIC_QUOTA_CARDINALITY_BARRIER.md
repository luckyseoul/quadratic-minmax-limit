# A cardinality barrier after punctured halved surjectivity

**Status:** proved for every branch-C prime \(p\ge31\). This is a negative
theorem about the sufficiency of punctured mod-two surjectivity and the scalar
direction quotas. It neither excludes nor constructs the actual branch-C
target, and residual (ii) remains open.

## 1. Exact setup

Put

\[
 p=2h+1,\qquad d=p+1,\qquad
 N=|\Delta|=dh=2h(h+1).
\]

After fixed-edge elimination, the unused nonfixed columns split into the
\(d\) disjoint classes indexed by their unique parallel direction. Write
\(c_L\) for the number of retained columns in class \(L\), and prescribe

\[
             0\le n_L\le c_L,\qquad s=\sum_L n_L.       \tag{1}
\]

The halved binary map \(D_U=(C_U,\Phi_U)\) has codomain dimension

\[
                         R=dh(h+1).                       \tag{2}
\]

Assume the separately proved structured-puncture result that \(D_U\) is
onto. Under the explicit fixed-edge inverse, the \(P_L\) functional pulls
back to the sum of the fixed classes \([v]\subset\ker L\). These \(d\)
nonempty sets partition \(\Delta\), so the parallel functionals are
independent on the transformed codomain. Every nonfixed column in group
\(L\) has parallel word \({\bf e}_L\). Consequently, after fixing the
parallel parity vector \((n_L\bmod2)_L\), there are exactly

\[
                         2^{R-d}                          \tag{3}
\]
possible target syndromes.

## 2. The information deficit

The number of Boolean source words with the exact weights (1) is

\[
                   Q(\mathbf n)=\prod_L {c_L\choose n_L}.\tag{4}
\]

Since the groups partition a subset of the full \(N^2\) columns,
Vandermonde's identity gives

\[
 Q(\mathbf n)
 \le {\sum_Lc_L\choose s}
 \le {N^2\choose s}
 \le (N^2)^s.                                             \tag{5}
\]

If \(s=0\), then \(Q=1<2^{R-d}\) immediately.  Assume \(s\ge1\).
For \(h\ge7\), one has \(N=2h(h+1)<2^h\): it holds at \(h=7\), and the
ratio of consecutive left sides is \((h+2)/h<2\). Along branch C the target
edge count satisfies \(|H|\le N-1\), and fixed-edge elimination gives
\(2s=|H|-|U|-|a|\le N-1\). Hence

\[
 Q(\mathbf n)
 <2^{2hs}
 \le2^{h(N-1)}
 <2^{R-d},                                                \tag{6}
\]

where the last inequality follows directly from

\[
 R-d-h(N-1)=d(h-1)+h>0.                                  \tag{7}
\]

The branch under study has \(p\ge31\), hence \(h\ge15\), so (6) applies
uniformly.

## 3. Consequence and exact limitation

Equations (3) and (6) prove:

> **Quota-cardinality barrier.** For every feasible exact direction-quota
> vector arising in the branch-C size range, some target syndrome with the
> same parallel-row parity has no Boolean preimage of those direction
> weights, even when \(D_U\) is onto.

Thus

\[
 \boxed{D_U\text{ onto}+0\le n_L\le c_L
 \quad\not\Longrightarrow\quad
 \exists b:\ D_Ub=\widehat T_U,\ |b|_L=n_L}              \tag{8}
\]

as a target-independent theorem. In fact the exact quota slice is
information-theoretically far too small to cover one fixed-parity target
fibre.

This does **not** show that the actual \(\widehat T_U\) is absent. That
target is highly structured, while the proof of (8) only shows that at least
one compatible syndrome is absent. Closure now requires genuinely
target-sensitive input: a kernel/exchange theorem specialized to the
transverse branch-C target, a direct obstruction for that target, or an
explicit construction. Unused-column abundance, punctured surjectivity, and
the scalar quota inequalities alone cannot supply the missing implication.

## Reproduction

The implementation evaluates only the exact formulas and binomial counts; it
does not build a finite Radon matrix or run a prime census.

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_symmetric_quota_cardinality_barrier.py
```
