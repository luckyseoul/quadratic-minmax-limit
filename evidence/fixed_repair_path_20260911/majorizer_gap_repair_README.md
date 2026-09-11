# Cone-ranked repair on the canonical-majorizer residual

This finite numerical method test applies the imported Navier--Stokes
positive-cone ranking and correction-cycle pattern to the actual residual in
the complete order-four cross-block diagnostic.  It uses every one of the
65,536 blocks from `gap_n4.json.gz`, with residual
`(relative canonical SDP gap, diagonal dispersion)`.  At every stage the
complete stored residual is recomputed after an accepted integral one- or
two-bit change.  Fractional cone coefficients only rank proposals.

Soulkiller ran all 65,536 starts for each of four modes with 88 workers.
The cone mode improved the numerical gap at 65,046 starts, with mean gap
reduction `0.0772161719286687`.  The matched index-pair and seeded
random-pair controls improved 64,754 and 64,776 starts, with reductions
`0.07128432576824603` and `0.07257528503244955`; single edge reached
63,552 and `0.0609467233652979`.  All modes reached the same already-known
best table value `1/42` numerically.

`majorizer_gap_repair_full.json.gz` is the complete 262,144-row output,
SHA-256 `c56c806624189cf0a8ac9d1de3094c688716ecd94fd53b9e4dca78b22d4670b4`.
NUKA independently replayed every row against the original table and every
reported aggregate with `scripts/verify_majorizer_gap_repair.py`.

This measures discrete local repair on an SDP *numerical diagnostic*.  It
does not provide an exact SDP certificate, establish small canonical gap for
any optimizer, or imply an all-orders convergence argument.
