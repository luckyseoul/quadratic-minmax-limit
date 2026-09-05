# Proposition 15.775: eventual first-layer closure by bounded integral quadratics

Date: 2026-09-04.

Classification: **proved infinite-family theorem**, with an explicit
superlinear bounded-support corollary. This is not an all-size theorem,
global residual-(ii) closure, an eventual E1 theorem, or a proof of the
original MathOverflow limit. No prime, graph, coefficient, boundary, or
slice census is used in the new proof.

The main conclusions are:

1. For every odd integer `p>=259201`, no nonnegative integer-valued
   quadratic `A` on `J(p,(p+1)/2)` has `2p E[A]=2p+4` or `2p+6`.
   No affine-parity hypothesis is needed.
2. For every prime `p>=259201`, the two signed shell floors three are
   impossible at `|H|=5p+6`. Thus the entire first residual layer left
   by 15.774, `k=5p+5, t_residual=(p-1)/2+3`, is closed eventually,
   for every boundary and both phases. The eventual next frontier is
   `k>=5p+7`; the smaller-prime frontiers are unchanged.
3. For every prime `p>=29`, `r in {3,4,5}`, and integer `h>=0` with
   `h=r (mod 2)`, both signed shell floors `r` are impossible if
   `h<rp`, or if `46656 h^3 <= p^3(p-1)`. Equivalently the latter
   support band is `h<=p(p-1)^(1/3)/36`. Its growing upper endpoint
   does not remove the remaining all-size quantifier.

The new ingredient is a dimension-free height bound, followed by an
elementary invariant-class bound and an exact middle-slice mean identity.
The first-layer application also uses the proved small-mass spectrum
from [15.774](NOTE_2026-09-04_SMALL_MASS_TWO_TYPE_BRIDGE.md).

## 1. A self-contained cube height bound

Let `f` be a nonnegative integer-valued polynomial of degree at most
two on a Boolean cube of any dimension. Write its sign-coordinate
Fourier expansion as `f=sum_S fhat(S) chi_S`.

Its multilinear polynomial in zero-one coordinates has integer
coefficients: each coefficient is a finite difference of integer
values at cube vertices. Changing coordinates to signs shows that
every Fourier coefficient belongs to `(1/4) Z`. Consequently

```text
|fhat(S)| <= 4 fhat(S)^2,
max f <= sum_S |fhat(S)| <= 4 ||f||_2^2.                (1)
```

For completeness, the needed degree-two norm inequality has the
following elementary noise proof. Set `rho^2=1/3`. On a two-point
space,

```text
E_s (a+rho*s*b)^4 = a^4+2a^2b^2+b^4/9
                  <= (a^2+b^2)^2.                    (2)
```

This tensorizes to `||T_rho g||_4<=||g||_2`: write a function as
`g(x,s)=a(x)+s b(x)`, use Cauchy--Schwarz for the mixed fourth-moment
term after applying noise in the other coordinates, and apply the
induction hypothesis separately to `a,b`. The resulting upper bound
is `(||a||_2^2+||b||_2^2)^2`; the fourth-power coefficient `1/9` in
(2) is at most one. The zero-dimensional case starts the induction.

On degree at most two, the inverse noise multiplier has `L2` norm
at most `rho^(-2)=3`. Hence `||f||_4<=3||f||_2`. Holder interpolation
gives

```text
||f||_2 <= ||f||_1^(1/3) ||f||_4^(2/3)
        <= ||f||_1^(1/3) (3||f||_2)^(2/3),
||f||_2 <= 9 ||f||_1 = 9 E[f].
```

The zero function is immediate; otherwise divide by the positive
power of its `L2` norm. Combining with (1),

```text
max f <= 324 (E[f])^2.                                (3)
```

This proof uses neither a classification of integer-valued quadratics
nor a bounded number of values as an assumption.

## 2. Transfer from slice mean to bounded height

Put `m=(p+1)/2`, `q=(p-1)/2`. Let `A>=0` be an integer-valued
quadratic on `J(p,m)`, with mean `mu` and maximum `H_A`, attained
at `X`. The stabilizer moment certificate of 15.642, also used in
15.688, gives

```text
H_A <= p*mu                         if p=3 (mod 4),
H_A <= p*mu*(p+3)/(p-1)             if p=1 (mod 4).
```

The second expression is therefore a uniform bound in both classes.
Indeed the stabilizer average is a nonnegative quadratic in
`z=|Y intersect X|`, with value `H_A` at `z=m`. The exact nonnegative
three-node quadrature has endpoint weight `1/p` in the first class
and `(p-1)/(p(p+3))` in the second. Its first three moments are those
of the hypergeometric intersection distribution. Dividing the mean
by that endpoint weight proves the displayed bound without any
support-size assumption.

Use the paired cubes through `X`: choose one singleton from `X`,
keep it selected, pair the other `m-1` points of `X` bijectively with
the complement, and select one point from every pair. Every such
cube contains `X`. Average its uniform mean over the singleton and
the matching. On every slice quadratic the result is

```text
T A(X) = (A(X)+p E[A])/(p+1).                         (4)
```

This follows on the constant, coordinate, and pair monomials. A
coordinate in/outside `X` has paired probability `(m+1)/(2m)` / `1/2`;
a pair inside `X` has probability `(m+2)/(4m)`; a pair outside or
crossing `X` has probability `1/4`. These are exactly (4) for the
corresponding slice moments. At least one cube through `X` therefore
has mean at most `(H_A+p*mu)/(p+1)`. Its restriction is a
nonnegative integral quadratic with the same maximum `H_A`.

For `2p*mu=2p+s`, `s in {4,6}`, the preliminary stabilizer bound
and (4) give a cube mean at most

```text
(2p+s)/(p-1) <= 16/7,                p>=29.
```

Apply (3) on this selected cube:

```text
H_A <= 324*((2p+s)/(p-1))^2
    <= 82944/49 < 1800,
mu <= 32/29 < 9/8.                                    (5)
```

We use the convenient uniform bound `H_A<=1800` below. Neither
the cube dimension nor the slice order occurs in this bound.

## 3. Integral derivatives and an invariant-class bound

For a transposition `(ij)`, define

```text
I_ij = (1/4) E[(A-A^(ij))^2].                         (6)
```

If `A-A^(ij)` is not identically zero, condition on `i in X, j notin X`.
The difference is then a nonzero integer-valued affine function `ell`
on `J(N,K)`, where `N=p-2`, `K=(p-1)/2`. If this affine function is
constant, its nonzero support is the whole slice. Otherwise choose
two unequal affine coefficients. The probability that a uniform
`K`-set contains exactly one of their positions is

```text
2K(N-K)/(N(N-1)) = (p-1)/(2(p-2)).
```

Swapping those positions changes `ell` by their nonzero coefficient
difference; at least one value in each such two-point orbit is
nonzero. Thus its conditional support has proportion at least
`(p-1)/(4(p-2))>1/4`. Integrality makes the squared difference at
every nonzero point at least one. The event `x_i != x_j` has
probability `(p+1)/(2p)`, and the other orientation has the same
conditional squared distribution. Therefore every relevant pair has

```text
I_ij >= (p^2-1)/(32p(p-2)) > 1/32.                    (7)
```

The Johnson transposition Laplacian has degree-`j` eigenvalue
`j(p-j+1)`. For a degree-two function its nonconstant eigenvalues
are at most `2(p-1)`. With normalization (6), this gives

```text
sum_(i<j) I_ij <= (p-1) Var(A) <= (p-1) H_A*mu.       (8)
```

The first statement can also be read directly on the slice spaces
spanned by degree-zero, degree-one, and degree-two monomials: the
Laplacian acts on their successive quotients with eigenvalues
`0,p,2(p-1)`, respectively. Self-adjointness gives the orthogonal
decomposition and the bound in (8).

Declare `i~j` when `(ij)` leaves `A` invariant, also including `i=j`.
This is an equivalence relation: invariance under `(ij)` and `(jk)`
implies invariance under `(ik)=(ij)(jk)(ij)`. Thus the relevant pairs
form a complete multipartite graph. If its class sizes are `n_b`,
the largest class has size `p-L`, and the number of relevant pairs is

```text
R=(p^2-sum_b n_b^2)/2 >= (p^2-p(p-L))/2 = pL/2.       (9)
```

Here `sum_b n_b^2 <= (max_b n_b) sum_b n_b` proves the inequality.
Combining (7)--(9), when `L>0`,

```text
L < 64 H_A*mu <= 64*1800*(9/8) = 129600.              (10)
```

The same final strict bound is immediate for `L=0`. Consequently
`L<=129599`. This uses real-valued influence bounds and integral
nonzero differences; it does not assume that `A` is Boolean.

## 4. Full-cube extension and its exact mean

Let `K_0` be a largest invariant class, and `J` its complement of
size `L`. Average a degree-two representative of `A` over permutations
of `K_0`. It retains the same slice values and has the form

```text
A_0(x_J) + A_1(x_J) s_K + c binom(s_K,2),
s_K=sum_(i in K_0) x_i,
```

where `A_0` has degree at most two and `A_1` degree at most one.
Substituting `s_K=m-s_J` gives a degree-two polynomial `f` on the
full cube in `J` and the same function on the slice. If `L<=q`,
every kept-coordinate pattern extends to a slice point: for any
`0<=s_J<=L`, the number `m-s_J` lies between zero and `p-L`.
Thus `f` is nonnegative and integer-valued everywhere, and
`0<=f<=H_A` on the full cube. Its cube mean `nu` belongs to `(1/4) Z`.

For `p>=259201`, (10) gives `L<=129599<q`, so this extension exists.
The exact middle-slice moments are

```text
E_slice[x_i]=(p+1)/(2p),
E_slice[x_i*x_j]=(p+1)/(4p),       i!=j.
```

Compare these with the cube moments `1/2,1/4`. Checking the constant
term too gives, for every such quadratic extension,

```text
p*mu=(p+1)*nu-f(0).                                  (11)
```

If `2p*mu=2p+s`, `s=4,6`, then

```text
nu=1+(f(0)+s/2-1)/(p+1).
```

The numerator is at least one and at most `1802`. Since
`p+1>=259202>4*1802`, this says `1<nu<5/4`, contradicting the
quarter-integral cube mean. This proves conclusion 1 for every odd
integer at or above the displayed threshold.

## 5. Exact first-layer quotas force the excluded masses

Now let `p>=37` be prime and suppose both signed shell floors three
hold for an edge set with `h=|H|=5p+6`. There is an isolated vertex,
since `p^2+1-2h>0`. Transport it to infinity, so the incident count
is zero, and set `M=p+1=2m`, `t_shell=(h-3p)/2=p+3=2m+2`.

For each direction of type `epsilon`, the existing exact identities
are

```text
A_d=(T_H^epsilon-3)/2 >= 0,
a_d=2p E[A_d]=M P_d-epsilon*T-3p=2u_epsilon+M k_d,
0<=u_epsilon<m,    k_d>=0,
sum_(d:type epsilon) k_d=t_shell-u_epsilon,
u_+ + u_- = 3 (mod m).                               (12)
```

These local quadratics have affine parity with an even boundary, but
the eventual local contradiction in section 4 does not require it.
Only the preceding quota reduction uses that affine-parity fact via
15.774: below `2p-10`, the union of possible masses is
`{0,p-3,p-1,p+1}`. It implies the necessary quotient floors

```text
k_min=0    at u=0,m-2,m-1;
k_min=2    at 1<=u<=m-7;
k_min=1    at m-6<=u<=m-3.                            (13)
```

The type capacity `t_shell>=u+m*k_min(u)` confines each residue to
`{0,1,2} union {m-6,...,m-1}`. A low-plus-high residue is at most
one or at least `m-6>3`; a high-plus-high residue is at least
`m-12>3`. Thus (12) has only the ordered pairs `(1,2)` and `(2,1)`.
For `u=1`, every quotient is at least two and their sum is `2m+1`:
exactly one is three and the other `m-1` are two. For `u=2`, their
sum is `2m`: every quotient is two. Their low rows therefore have
the respective masses

```text
2+2M=2p+4,             4+2M=2p+6.                    (14)
```

For prime `p>=259201`, either mass contradicts section 4. This
closes the entire layer, without choosing or classifying a boundary.

An important normalization point: (12) does not force the low
parallel count to be five or the transported signed total to be
`+/-1`. If the `u=1` type is `sigma`, with low parallel count `P`,
then its parallel counts are `P` and one `P+1`; the opposite type
has constant parallel count `Q`, with

```text
P+Q=10,       sigma*T=(p+1)(P-5)+1.                   (15)
```

The proof never assumes `P=5` and works for every possible phase
of every local row. The two signed shell types are both retained.

## 6. A quantitative superlinear support corollary

Keep odd `p>=29`, but replace the fixed local masses by any
nonnegative integral slice quadratic of mean at most an exact
positive number `B`. For a nonzero quadratic, sections 1--3 give

```text
paired cube mean <= 2p*B/(p-1) <= 3B,
H_A <= 2916 B^2,
L < 64 H_A*B <= 186624 B^3.                           (16)
```

The zero-mean case is the zero function: it satisfies the final
height bound `H_A<=2916 B^2` and strict junta bound `L<186624 B^3`
directly for every positive cap `B`. Thus, if `186624 B^3<=q`, every such
quadratic has the full-cube extension from section 4.

Let `p>=29` be prime and suppose an edge set of size `h` has both
signed shell floors `r in {3,4,5}`, with the indispensable parity
condition `h=r (mod 2)`. If `h<rp`, the signed frame mean already
contradicts these floors. Otherwise set `B=h/(2p)>=r/2>=3/2` and
suppose

```text
46656 h^3 <= p^3(p-1),
equivalently 373248 B^3 <= p-1.                       (17)
```

This implies `4B<=p-1`, so `2h=4pB<p^2+1` and there is an
isolated vertex. Use the same transported chart as in section 5.
For the baseline-`r` local rows,

```text
A_d=(T_H^epsilon-r)/2,
a_d=2p E[A_d]=M P_d-epsilon*T-rp,
average_(d:type epsilon) E[A_d]=(h-rp)/(2p)<=B.        (18)
```

Integrality follows from the stated parity condition, not from a
choice of phase. Choose one row of mean at most `B` from each type.
By (16)--(17), both chosen rows have full-cube extensions, say
`f_+,f_-`. Put `alpha_epsilon=f_epsilon(0)` and let their cube
means be `nu_epsilon`. Then

```text
0<=alpha_epsilon<=2916 B^2,
4nu_epsilon in Z,
a_epsilon=2M nu_epsilon-2alpha_epsilon.               (19)
```

Add the two signed identities in (18) and eliminate `T` using
(19). Rearrangement gives the exact divisibility identity

```text
4(alpha_++alpha_-+r)
  = M(4nu_++4nu_- -2(P_++P_-) +4r).                  (20)
```

The parenthesis on the right is an integer. The left is positive,
but, since `r<=2B<=2B^2`,

```text
0 < 4(alpha_++alpha_-+r)
  <= 4(5832 B^2+r)
  <= 23336 B^2 <= 23336 B^3
   < 373248 B^3 <= p-1 < M.                          (21)
```

There is no positive multiple of `M` in this interval. This proves
conclusion 3. Both signs enter (20); neither type, parity phase, or
boundary was discarded. The argument is the same bounded-quadratic
proof family, with explicit constants, not a new equality catalog.

## 7. Scope and executable receipts

The main first-layer threshold uses the sharper fixed-mass bound (5);
it does not require the looser support condition (17). For orders
below `259201`, conclusion 2 makes no change to the 15.774 frontier.
The superlinear corollary independently applies precisely in its
stated numerical support band. Neither result proves that all
remaining residual or minimal-gap supports lie inside that band.

The implementation is
[`src/e1_gmin_m4_prop15775.py`](../src/e1_gmin_m4_prop15775.py), with
[`tests/test_prop15775.py`](../tests/test_prop15775.py). The generated
receipt is `evidence/e1_gmin_m4_prop15775.json`. Its two first-layer
records, `p=524287,6700417`, cover both residue classes modulo four;
three corollary records use `p=6700417`, `h=rp+2`, `r=3,4,5`.
These are exact identity/inequality replays, not a census proving the
universal quantifier. The proof above supplies that quantifier.

Only a fixed small receipt at `p=37` is used to validate the imported
15.774 spectrum contract. The large-order implementation checks a
constant number of residue interval endpoints; it does not construct
the old order-`p` boundary or residue tables at the new orders.
The stabilizer dependency checks its three weights and moments
directly, without its old slice enumeration routine.

Replay on an authorized compute worker, not on the controller:

```sh
PYTHONPATH=src python src/e1_gmin_m4_prop15775.py
PYTHONPATH=src python -m pytest -q tests/test_prop15775.py
```

The tests check normalization, both congruence classes, exact
thresholds, parity and support guards, dependency-failure injections,
the saved receipt, and explicit false global-closure flags. This note
does not claim an execution result; the offloaded validation receipt
and canonical handoff record the actual run separately.
