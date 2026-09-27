# Local-lemma balanced-core one-vertex extension

2026-09-27. Direct strengthening of the balanced-core hypergeometric
one-vertex theorem under the scope freeze. The previous theorem used a union
bound over all residual active stable states. Here the independent random
objects are the balanced signings of the core-signature classes, so bad
events are only locally dependent. The Lovasz local lemma replaces the
global residual-state count by an overlap degree.

Let A be a complete symmetric zero-diagonal signing of order n, let
F=Phi(A), fix a desired increment r>0, and let

    H={s : D_s+r<=n}

be the active antipodal stable family.

Choose a core C subset H, form projective core-signature classes, remove one
coordinate from each odd class, and let U be the leftover set with

    u=|U|<r.

For every remaining even class G choose the canonical signs lambda_i exactly
as in NOTE_2026-09-27_BALANCED_CORE_HYPERGEOMETRIC_EXTENSION.md.

For a residual state s in H\C and an even class G, define

    q_(s,G)
      = min( #{i in G: lambda_i x_i^s=+1},
             #{i in G: lambda_i x_i^s=-1} ).

Put

    J_s={G:q_(s,G)>0},
    q_s=sum_G q_(s,G),
    tau_s=D_s+r-u.

Randomly and independently for each even class G choose a uniformly balanced
sign vector b^G, and put a_i=lambda_i b_i^G on G. Fix arbitrary signs on U.

As proved in the balanced-core note, for each residual state

    P(E_s)
      :=P{|Z_s|>=tau_s}
      <= p_s
      :=2 exp[-tau_s^2/(8q_s)]                              (1)

when q_s>0; if q_s=0 then E_s is impossible.

Crucially, E_s depends only on the independent class variables
{b^G:G in J_s}. Therefore E_s is independent of every collection of events
whose impurity supports are disjoint from J_s.

Define the dependency graph on residual states by

    s~t  iff  J_s intersect J_t is nonempty.

## Theorem 1: asymmetric local-lemma certificate

Suppose there are numbers z_s in (0,1) for all residual states with q_s>0
such that

    p_s <= z_s product_(t~s) (1-z_t).                       (2)

Then with positive probability no bad event E_s occurs. Consequently there
exists an actual incident row a in {+-1}^n satisfying

    Phi(A extended by a)<F+r.                               (3)

Proof. The random balanced signing of each even signature class is one
independent product variable. Event E_s is measurable only with respect to
the variables indexed by J_s. The standard asymmetric Lovasz local lemma
applies to the dependency graph above and gives positive probability that no
E_s occurs. Core states are annihilated up to the u leftovers exactly as in
the balanced-core theorem, inactive stable states are automatically harmless,
and the exact stable-state extension identity gives (3). QED.

## Theorem 2: symmetric usable form

Let

    Delta_dep
      = max_s #{t != s : J_s intersect J_t != empty}.       (4)

If

    boxed:
    max_(s:q_s>0)
      2 exp[-tau_s^2/(8q_s)]
      <= 1/[e(Delta_dep+1)],                                (5)

then (3) holds.

Equivalently, it is enough that every residual active state obey

    boxed:
    tau_s^2
      >= 8 q_s log[2e(Delta_dep+1)].                        (6)

This strictly replaces the union-bound denominator log(2M), where M was the
TOTAL number of residual active states, by log[2e(Delta_dep+1)], where only
states sharing an impure signature class interact.

## Theorem 3: class-load form

For an even signature class G define its residual impurity load

    L_G=#{s in H\C : q_(s,G)>0}.

For a residual state put

    h_s=|J_s|.

Then

    deg(s)
      <= sum_(G in J_s)(L_G-1).                             (7)

Hence if

    Lambda
      :=max_s sum_(G in J_s)(L_G-1),                        (8)

the simpler certificate

    boxed:
    tau_s^2
      >= 8q_s log[2e(Lambda+1)]
      for every residual s with q_s>0                       (9)

constructs a valid extension row.

In particular, if h_s<=h for every residual state and L_G<=L for every
class, then Lambda<=h(L-1), and it is enough that

    boxed:
    tau_s^2
      >=8q_s log[2e(h(L-1)+1)]                              (10)

for every residual active state.

## Constructive version

Because the probability space is a product over signature classes, the
Moser--Tardos resampling algorithm applies directly. Start with independent
uniform balanced signings b^G. While some event E_s is violated, resample
all class variables b^G with G in J_s. Under the standard local-lemma
criterion (2), this terminates in finite expected resampling time and returns
an explicit incident row a satisfying (3).

Thus the theorem is not only existential: a successful finite certificate
can be turned into an actual one-vertex extension row without searching the
2^n row space.

## Neutral consequence

For an exact order-n minimizer take

    r_n=m_n[(1+1/n)^(3/2)-1].

If a core C satisfies u<r_n and either (2), (5), or (9), then

    boxed: alpha_(n+1)<alpha_n.                             (11)

More generally use r=r_n+rho_n. If such a certificate exists for all large
n and

    sum_n rho_n/(n+1)^(3/2)<infinity,

then the frozen one-vertex argument proves convergence of alpha_n.

## Why this advances the live target

The previous balanced-core theorem paid for every residual state globally:

    sum_s p_s<1.

That criterion fails whenever there are many residual constraints, even if
they occupy almost disjoint signature classes.

The present theorem pays only for LOCAL impurity overlap. It can therefore
certify extension with arbitrarily many residual active states provided each
one interacts with only a controlled neighborhood in the signature-class
hypergraph.

The new finite discriminator is now the triple

    (q_s, J_s, L_G),

not the total number of active states. A computation need only build the
impurity incidence graph and test (2) or (9); success is a complete
one-vertex certificate.
