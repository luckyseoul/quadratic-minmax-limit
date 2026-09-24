# Komlós vector balancing gives a new one-vertex extension theorem

2026-09-24. Cross-order analytic result using the September 2026
resolution of the Komlós vector-balancing conjecture by Guo--Fang--Lu.
The original convergence problem remains OPEN. No finite signing census is
used.

The point of this note is that the one-vertex extension problem is itself a
vector-discrepancy problem once the exact stable-state reduction is applied.
The new dimension-free Komlós theorem changes the relevant near-edge
complexity measure from a union-bound probability mass to a reciprocal-square
mass.

## 1. External theorem used

Guo, Fang and Lu, "Vector Balancing via Directional Total Variation",
arXiv:2609.11189 (2026), prove:

> If v_1,...,v_N are finitely many real vectors with ||v_i||_2<=1, then
> there are signs eps_i in {+1,-1} such that
>
>     || sum_i eps_i v_i ||_infty < K,
>
> where
>
>     K = 3 sqrt(2 pi).

Thus K^2=18 pi. Their proof is dimension-free and cardinality-free.
The arXiv abstract states that the proof was discovered by the Odin
Automatic AI Research Agent. The theorem, not that provenance, is what is
used below.

## 2. Compress the exact extension constraints

Let A be any complete symmetric zero-diagonal signing of order n and put

    F=Phi(A).

For sigma in {+1,-1}, call x sigma-stable when

    sigma x_i (Ax)_i >= 0  for every i.

The previously proved stable-state reduction gives, for every incident
sign row a,

    Phi(A extended by a)
      = max_(sigma,x sigma-stable)
          [sigma Q_A(x)+|a.x|].                            (1)

Identify x with -x. If an antipodal class [x] is stable in one or both
orientations, define its best stable oriented energy

    E_A([x])
      = max { sigma Q_A(x) :
              sigma in {+1,-1}, x is sigma-stable },

and its stable deficit

    D_A([x])=F-E_A([x]) >= 0.                              (2)

Let T_A be the finite collection of these distinct antipodal classes.
Then (1) becomes exactly

    Phi(A extended by a)-F
      = max_([x] in T_A) [ |a.x|-D_A([x]) ].               (3)

No multiplicity from the two orientations remains.

## 3. Reciprocal-square extension theorem

Fix a desired increment r>0. A class with

    D_A([x])+r > n

is automatically harmless, because |a.x|<=n for every sign row a.
Let

    T_A(r)={ [x] in T_A : D_A([x])+r<=n }

be the active stable classes, and define

    Xi_A(r)
      = sum_([x] in T_A(r)) 1/[D_A([x])+r]^2.              (4)

Then:

    THEOREM.
    If

        Xi_A(r) <= 1/(18 pi),                              (5)

    there exists an incident sign row a such that

        Phi(A extended by a) < F+r.                        (6)

Proof. If T_A(r) is empty there is nothing to prove. Otherwise set
S=Xi_A(r)>0. For each old vertex i form a vector indexed by active stable
classes,

    v_i([x])
      = x_i / [(D_A([x])+r) sqrt(S)].                       (7)

Since x_i^2=1,

    ||v_i||_2^2
      = (1/S) sum_([x] in T_A(r)) 1/[D_A([x])+r]^2
      =1.                                                   (8)

Apply the Guo--Fang--Lu Komlós theorem to v_1,...,v_n. There are signs
a_i such that, simultaneously for every active [x],

    | sum_i a_i v_i([x]) | < K.

Therefore

    |a.x|
      < K sqrt(S) [D_A([x])+r]
      <= D_A([x])+r,                                      (9)

using K sqrt(S)<=1. For inactive classes the inequality is strict automatically from
|a.x|<=n<D_A([x])+r. Thus every class has
|a.x|-D_A([x])<r. Insert this in (3), proving (6). QED.

This is an ACTUAL sign row, not a fractional or Gaussian row.

## 4. Exact neutral-increment consequence

Now let A_n be an exact order-n minimizer, F=m_n, and define the exact
neutral raw increment

    r_n
      = m_n [ (1+1/n)^(3/2)-1 ].                           (10)

It satisfies

    m_n+r_n = alpha_n (n+1)^(3/2).                         (11)

Hence if an exact minimizer A_n obeys

    Xi_(A_n)(r_n) <= 1/(18 pi),                            (12)

then the theorem gives

    m_(n+1) <= m_n+r_n,

so

    boxed: alpha_(n+1) <= alpha_n.                         (13)

This is a genuine cross-order implication.

More generally, for rho_n>=0, if

    Xi_(A_n)(r_n+rho_n) <= 1/(18 pi)                       (14)

and

    sum_n rho_n/(n+1)^(3/2) < infinity,                    (15)

then

    alpha_(n+1)
      <= alpha_n + rho_n/(n+1)^(3/2).                      (16)

Summable positive rises force alpha_n to converge. Thus (14)--(15) are a
new Dini-summable convergence criterion expressed entirely in the stable
near-edge spectrum of exact minimizers.

## 5. A new obstruction theorem for every upward step

The contrapositive is unconditional and useful.

If, for an exact minimizer A_n,

    alpha_(n+1) >= alpha_n,                                (17)

then necessarily

    boxed: Xi_(A_n)(r_n) > 1/(18 pi).                      (18)

So even a FLAT normalized step is impossible unless the minimizer carries a
definite reciprocal-square mass of one-flip-stable near-edge states.

Since each term in Xi is at most 1/r_n^2, (18) also forces

    |T_(A_n)(r_n)| > r_n^2/(18 pi).                        (19)

Using m_n=alpha_n n^(3/2),

    r_n
      = (3/2) alpha_n sqrt(n) + O(alpha_n/sqrt(n)),

and therefore every nondecreasing step satisfies asymptotically

    |T_(A_n)(r_n)|
      > [alpha_n^2/(8 pi)+o(1)] n.                         (20)

With the already-proved uniform positive liminf, this is a genuine linear
degeneracy requirement. An increasing alpha step cannot be supported by
only O(1), o(n), or even an arbitrarily small fixed number of dangerous
stable antipodal classes.

The active window in (19) is explicit:

    D_A([x]) < n-r_n.                                      (21)

Thus these are states within O(n) raw energy of the Boolean edge
F=Theta(n^(3/2)), not arbitrary local maxima deep in the cube.

## 6. Dyadic shell form

For k>=0 let N_k count active stable antipodal classes with

    k r <= D_A([x]) < (k+1)r.

Then

    Xi_A(r)
      <= (1/r^2) sum_(k>=0) N_k/(k+1)^2.                  (22)

Consequently the concrete shell condition

    sum_(k>=0) N_k/(k+1)^2
      <= r^2/(18 pi)                                      (23)

is sufficient for an increment-r extension.

This makes the gain over the old random-row union criterion transparent.
At deficit zero, the old exact-tail union argument charges a constant
probability per antipodal maximizer at r=Theta(sqrt(n)), so it can certify
only O(1) such constraints. Komlós charges 1/r^2=Theta(1/n) per zero-deficit
constraint and therefore tolerates a LINEAR family before the certificate
fails. Farther shells are discounted quadratically.

## 7. Scope

This does not prove that (12) holds for every minimizer, so convergence is
not yet proved. It does, however, insert a new September-2026 discrepancy
theorem directly into the missing cross-order step and yields:

1. a deterministic one-vertex sign-row theorem;
2. a Dini-summable convergence criterion;
3. an unconditional structural theorem saying every nondecreasing alpha
   step forces linear stable near-edge degeneracy.

In particular, if alpha_n fails to converge, it must have infinitely many
nondecreasing steps (an eventually strictly decreasing bounded sequence
would converge), so the linear-degeneracy obstruction must occur at
infinitely many orders.

This is independent of the same-order repeated-update lower-bound ladder.
It attacks the other missing side: movement between n and n+1.

Reference:
Shengtao Guo, Ethan X. Fang, Junwei Lu,
"Vector Balancing via Directional Total Variation",
arXiv:2609.11189 (submitted 10 September 2026).
