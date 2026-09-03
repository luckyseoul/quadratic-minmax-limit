# Global Mobius incidence cannot close the fixed-word collision gate alone

Date: 2026-09-03

Status: an all-prime probabilistic theorem proves that, for every branch
prime and every prescribed set of nonzero hard centers, some choice of the
auxiliaries has pairwise-disjoint physical supports while its fixed-word
hyperbolas have collision surplus at least `2m`.  Thus distinct target
directions, signed conic intersections, arbitrary higher intersections, and
physical ternarity do not by themselves contradict the exact demand
`sigma >= kappa_z+m+q`.  A closure has to couple those intersections to the
required physical cancellations or to the target atom/even-moment equations.
This note does not construct the endpoint target and does not close residual
(ii).

## 1. Exact random-auxiliary incidence law

Let `p=2m-1`, `h=m-1`, and fix any `m` distinct projective target
directions `L_i` with arbitrary nonzero centers `j_i`.  For each target,
choose its auxiliary functional uniformly and independently among the
`p(p-1)` functionals independent of `L_i`.

A nonzero affine block is a functional `N` modulo `N ~ -N`.  If its
projective direction differs from `L_i`, the normal form

    N_z=-(z+1)L_i/j_i + z^2 M_i/[j_i(z-1)],  z != 0,1,

solves uniquely for `M_i` at every `z`.  This gives `p-2` auxiliaries for
`N` and another disjoint `p-2` for `-N`; one conic never contains both.
Hence

    rho = P(the i-th half hits [N]) = 2(p-2)/(p(p-1)).       (1)

If `[N]` is parallel to `L_i`, the probability is zero.  Every projective
direction contains `h` block classes.  There are `m` target directions and
`m` other directions.  Writing

    F_n = E floor(Bin(n,rho)/2)
        = (n rho - (1-(1-2rho)^n)/2)/2,

and `d_b` for the number of half-conics through block `b`, the raw block
collision surplus

    S = sum_b floor(d_b/2)                                  (2)

therefore has exact expectation

    E S = hm(F_m+F_(m-1)).                                  (3)

This calculation already includes the signed equality `N_i=+/-N_k` and
all triple and higher block intersections.

## 2. Physical overlaps are rare for a different exact reason

For two distinct targets, normalize their arbitrary auxiliaries by the
coordinates `(q,A)` and `(r,B)` from the symbolic two-half intersection
classification.  Each coordinate pair runs bijectively over
`F_p x F_p^*`, hence over all `p(p-1)` auxiliaries.

For either physical orientation sign, the direct matching has a zero
branch with `(p-1)^2` auxiliary pairs and a nonzero branch with
`(p-2)(p-3)` pairs.  The swapped matching has another `(p-2)(p-3)` pairs.
Thus the expected number of common physical inversion orbits of two random
halves is exactly

    gamma = [6p^2-24p+26]/[p^2(p-1)^2].                    (4)

Let `T` sum this common-orbit count over all pairs of halves.  Then

    E T = binom(m,2) gamma < 3/4.                           (5)

In particular, Markov gives `P(T>0)<=E T`.  The event `T=0` says that the
physical half supports are pairwise disjoint, so their sum is automatically
ternary.  Notice that a common fixed-word block does not imply a common
physical edge: the block only says that both edges join the two affine lines
`N=-1` and `N=+1`.

## 3. Many block collisions coexist with disjoint physical supports

Changing one auxiliary replaces exactly `p-2` distinct block incidences.
After common old/new blocks are deleted, the numbers removed and added are
equal.  Each removal and addition changes one summand of (2) by zero or one,
so the absolute change in `S` is at most `p-2`, not `2(p-2)`.
McDiarmid therefore gives

    P(S<=2m-1)
      <= exp(-2(E S-(2m-1))^2/[m(p-2)^2]).                 (6)

Here `a=1-2rho >= 7/8`.  Also

    a^(m-1) > (1-2/(m-1))^(m-1) >= (13/15)^15 > 1/10,

where the middle sequence is increasing: the derivative of
`x log(1-2/x)` is `y-log(1+y)>0` with `y=2/(x-2)`.  Expanding (3) gives

    E S = m(m-2)/2 + m(m-1)a^(m-1)(1+a)/4
        > m(35m-67)/64.                                   (7)

For `m>=18`, substituting (7) into (6) makes the exponent exceed `3/2`.
After writing `n=m-18`, the cleared difference is

    1225n^4+62262n^3+1060117n^2+6169740n+2097892 > 0.

Thus this bad-event probability is below `exp(-3/2)<1/4`, while (5) is
below `3/4`.

At the least order `m=16`, use the sharper
`a^15>(7/8)^15>1/8`; indeed `7^15>8^14`.  Equation (6) then has exponent
greater than `4/3`.  Moreover

    E T = 10096/14415,
    exp(-4/3) < 81/293 < 1-10096/14415.

The first exponential inequality follows from the first four positive
terms of its power series:
`exp(4/3)>293/81`.  Therefore, at every branch order `m>=16`, the union of
the two bad events has probability strictly below one.  Some auxiliary
choice simultaneously satisfies

    T=0,                    S>=2m.                         (8)

For this choice every raw nonzero occurrence is a distinct surviving
physical orbit.  Hence the collision variable in the fixed-word ledger is
exactly

    sigma=(|U_np|-C)/2=sum_b floor(d_b/2)=S>=2m.           (9)

## 4. Exact consequence and exact limitation

When `kappa_z=0`, the fixed-word demand is

    sigma >= m+q,             0<=q<=m.                    (10)

The ternary family (8)--(9) reaches at least `2m`, so it is compatible with
the full range of (10).  Consequently no upper bound using only the `m`
distinct target hyperbolas, their signed intersections (including triple
intersections), and physical ternarity can contradict (10).

The witness from the probabilistic theorem has no physical cancellations.
It therefore does not meet the separate branch-C size floor
`kappa>=kappa_0`, even at the top endpoint.  That is precisely the missing
information: a viable continuation must couple the large block collision
surplus to the locations and signs of the *required* nonzero physical
cancellations, or impose the even common moments and prescribed nonfixed
target cells.  Generic Bezout/incidence energy, even globally over all
`m` conics, has been exhausted as a standalone route.

## Reproduction

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_mobius_global_incidence_barrier.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_mobius_global_incidence_barrier.py
```

These commands replay closed formulas.  They do not enumerate primes,
targets, auxiliary families, or physical supports.
