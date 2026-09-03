# Direction-block structure of the halved symmetric row code

Date: 2026-09-03

Status: exact row-code normal form, exact affine-block branch number, and
explicit non-rectangle words of weight $|Delta|$ proved. The minimum
distance, the words attaining weight $ph$, the full classification through
weight $|Delta|$, structured Mobius puncture robustness, residual (ii), E1,
and the original limit remain open.

## 1. The affine-block decomposition of point space

Let $p=2h+1$, $d=p+1$, and let

    Delta=(F_p^2 minus {0})/{+1,-1},       N=|Delta|=dh.

For a projective direction $A$, let $B_A$ be the binary span of the $h$
paired non-origin affine-line blocks of direction $A$. These blocks are
disjoint and each contains $p$ point classes.

The paired-affine-line theorem in
NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md proves that the square
point-block incidence matrix $M$ satisfies

    M M^T=M^T M=I.                                           (1)

Grouping its columns by direction gives the orthogonal direct sum

    F_2^Delta = direct_sum_A B_A,       dim B_A=h.            (2)

Write $l_A$ for the $h$ point classes on the vector line $A$, and let
$s_A$ be the sum of the $h$ affine blocks in $B_A$. They partition the
complement of $l_A$, so

    l_A=1+s_A.                                                (3)

## 2. Exact normal form of Row(D)

View a word on $Delta$ times $Delta$ as a matrix whose left coordinate is
the midpoint class $[a]$ and whose right coordinate is the difference
class $[delta]$. In direction $A$, the raw halved rows are:

* the parallel row $1$ tensor $l_A$;
* fixed-transverse rows $l_A$ tensor $b_J$, for $b_J$ in $B_A$;
* paired rows $b_I$ tensor $b_J$, for $b_I,b_J$ in $B_A$.

Equation (3) gives

    l_A tensor b_J
      = 1 tensor b_J + s_A tensor b_J,                       (4)

and conversely

    1 tensor b_J
      = l_A tensor b_J
        + sum_(b_I in the block basis of B_A) b_I tensor b_J. (5)

Because all affine blocks over all directions form a basis of
$F_2^Delta$, (4)-(5) prove

    Row(D)
      = (<1> tensor F_2^Delta)
        direct_sum direct_sum_A (B_A tensor B_A).             (6)

The sum is direct. In the $M$-basis, a vector $1$ has a nonzero component
in every $B_A$. Hence a matrix $1 q^T$ cannot be block-diagonal by
direction unless $q=0$.

The dimensions in (6) are

    dh + d h^2 = d h(h+1),                                   (7)

agreeing with the full halved-map rank.

Equivalently, after the $M$ change of basis a row-code word has the form

    1 q^T + T,                                                (8)

where $T$ is an arbitrary matrix block-diagonal with respect to (2).
This is the exact remaining minimum-support problem.

## 3. Exact branch number of M

For a nonzero point word $x$, put

    s=wt(x),       t=wt(M^T x).

If $S=supp(x)$ and $n_B=|S intersect B|$, then every point lies on $p$
affine blocks. Two distinct noncollinear point classes lie together in
two blocks, while two distinct collinear classes lie together in none.
Writing $c$ for the number of collinear unordered pairs in $S$ gives

    sum_B n_B = ps,
    sum_B n_B^2 = ps + 4 (binom(s,2)-c).                      (9)

For every nonnegative integer $n$,

    1_(n odd) >= 2n-n^2.

Therefore

    t >= s(p-2s+2)+4c.                                       (10)

When $1<=s<=h$, subtracting $p+1-s$ from the right side of (10) leaves

    (s-1)(p-2s+1)+4c >=0.

Thus $s+t>=p+1$ on the small side. If $s>=h+1$ and $t>=h+1$, the same
inequality is immediate. If $t<=h$, apply the small-side argument to
$M^T x$ and use (1) to recover $x$. Hence

    wt(x)+wt(M^T x) >= p+1.                                  (11)

A point has transform weight $p$, and an affine block has transform
weight one, so the branch number is exactly $p+1$.

This one-dimensional branch theorem is useful but does not by itself
prove a $ph$ lower bound for (8): the block-diagonal matrices couple many
rows, and the two weights in (11) may be distributed differently from row
to row.

## 4. Words at and below the Delta threshold

The fixed-transverse coordinates already give the $N$ rectangles

    l_A tensor b_J,          wt(l_A tensor b_J)=hp=N-h.       (12)

These remain the known candidates for minimum words.

There are also two exact families at weight $N$ which are not rectangles
of type (12).

First, the boundary summand in (6) contains every vertical fibre

    1 tensor e_delta,        wt(1 tensor e_delta)=N.          (13)

There are $N$ such words.

Second, for a nonzero scalar class $[c]$ modulo sign, let

    S_c(a,delta)=1 if [a]=[c delta], and 0 otherwise.         (14)

Scalar multiplication preserves every projective direction and permutes
the $h$ affine blocks inside each $B_A$. Its permutation matrix is
therefore block-diagonal in (2), so $S_c$ belongs to the second summand of
(6). It has one nonzero entry in every row and column and hence

    wt(S_c)=N.                                                (15)

There are $h$ distinct scalar graphs.

Consequently, the statement that every nonzero word of weight at most
$|Delta|$ is a fixed-transverse rectangle is false. Any correct
classification through that threshold must include at least (12)-(15).

Nothing here disproves the sharper conjecture

    minimum distance of Row(D) = ph,

or the possibility that the words of minimum weight are exactly (12).
Those assertions remain open.

## 5. The live obstruction

For a puncture $U$, surjectivity fails exactly when $U$ contains the
support of a nonzero row-code word. The next exact task can now be stated
without raw Radon coordinates:

* bound the support of $1 q^T+T$ when $T$ is block-diagonal as in (8);
* classify all such words through weight $N$; and
* compare those supports with the actual structured Mobius set $U$.

The already-known $ph$ rectangles cannot be contained in the structured
union under the imported midpoint-direction bound. The vertical fibres
and scalar graphs likewise expose concrete support shapes that a complete
structured-puncture proof must exclude. No such full exclusion is claimed
here.

## 6. Fail-when-wrong replay

The implementation constructs only $p=3,5,7$ matrices. It checks (1), that
the raw and normal-form generators have the same rank, that the vertical
fibres and scalar graphs lie in the raw row code with weight $N$, the point
and block equality witnesses for (11), and the second-moment identity (9)
on a fixed tiny subset. These are formula checks, not a prime census or
theorem evidence.

Reproduction:

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_symmetric_halved_row_code.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_halved_row_code.py
