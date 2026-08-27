# Complete p=7 six-finite boundary exclusion

Proposition 15.661 closes both product signs of the residual boundary with
six finite points at `p=7`. It does not close boundaries of size at least
eight, residual (ii), Type I, R1, global QVAR, or the limit.

For six finite boundary points, every direction has odd-fibre count in
`{0,2,4,6}`. The exact scaled floor tables are
`(0,8,8,8)` in phase zero and `(14,6,14,6)` in phase one. Each quadratic
type has four directions and exact budget 32. A complete V100 integer sweep
of all `C(49,6)=13,983,816` boundaries leaves 3,856,300 floor survivors.
The square-semilinear stabilizer reduces these to 80,704 orbits.

Of those, 80,519 ordinary orbits have type-floor sums only 24 or 32. A
sum-24 type elevates exactly one direction by eight units, so each type has
at most one non-singleton complete Johnson catalog. Across 160,745 exact
elevation cases, simultaneous catalog signatures in the 120-dimensional
left kernel over `F_3` and the 135-dimensional left kernel over `F_7`
exclude every case. The same catalog row is used in both fields; this is a
single joined necessary condition, not two unrelated survivor counts.

The remaining 185 deep orbits have floor pairs `(32,16)` or `(24,8)`. A
compact exact model represents all 280 slack values directly, imposing their
fixed parity, the 14 primitive integer degree-two equations in each
direction, exact scaled means and type sums, the score bound `0<=A<=13`, and
both prime-field dependency systems. This closes 92 orbits immediately. The
93 timeouts split into all 930 exact directional-mean allocations: 810 are
certified infeasible directly. The remaining 120 leaves all use only the
complete 36-, 1,764-, or 2,233-row low catalogs; exact two- and three-table
modular hash joins leave zero tuples.

NUKA independently repeats the complete NumPy floor sweep, reproducing the
V100 survivor hash and full histogram. Its serialized quotient reproduces
the ordered 80,704-orbit catalog and profile histogram byte-for-byte under
canonical JSON, and its serialized ordinary modular run reproduces all
160,745 rejections. A nonsquare signed anti-isometry fixes infinity, finite
zero, and the distinguished edge while reversing the conference sign, so
the `c_H=-1` exclusion transfers to `c_H=+1`.

Raw records are archived at
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-27-p7-six-finite/`.
The compact global audit has SHA-256
`d0e4de3041fc875090012f4091b0b57a75d0399c2d10550636089138ba50f6cb`.

Targeted arXiv searches found adjacent generalized-Paley graph work but no
duplicate of this residual classification. Exact OEIS searches find 80,704
in A060716, A133751, and A133756, and 160,745 in A254067, all for unrelated
matrix, factorial/gamma, or Collatz-array constructions. The value 3,856,300
has no exact OEIS hit. These are context checks only; no sequence claim is
made.
