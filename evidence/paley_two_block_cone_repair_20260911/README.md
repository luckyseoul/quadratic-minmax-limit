# Exact two-block cone-repair diagnostic

This is a finite test of the cross-block construction

\[
K(B)=\begin{pmatrix}C&B\\B^T&-C\end{pmatrix},
\]

where `C` is the order-14 (`q=13`) Paley conference signing.  It uses the
imported Navier correction pattern only to rank integral flips of the 196
entries of `B`.  A proposed block is scored by an exhaustive integer
enumeration of all `2^27` projective Boolean states of `K(B)`.

For deterministic seed 731, Soulkiller used 88 OpenMP workers.  At the
initial exact norm 94, its stress-cone ranking considered all 16 ranked
single flips and all 120 ranked pairs.  None strictly lowered the complete
Boolean norm.  NUKA independently compiled the scorer and replayed the
retained block with 16 workers, obtaining the same norm and witness:

```json
{"n":28,"source_phi":21,"phi":94,"energy":94,"xmask":7791,"ymask":2585}
```

| item | SHA-256 |
| --- | --- |
| Soulkiller result | `ef85e6dcd636bd586f736b98ad6355c811a14b461f980bbc9887644c4d90383d` |
| exact scorer | `9704efe6838667e553d8e4e2bf7fa800dd99a5c5daae7f299bec14e195e93a6a` |
| cone repair driver | `3a328419b2ec98c17fe9854d8a98bcd2d2b325281bdb7e8e7c45833bc600df97` |

## Interpretation boundary

The retained block provides only the trivial finite statement `m_28<=94`.
The stalled seeded cone neighborhood neither bounds the optimum over cross
blocks nor rules out another construction, and it has no all-orders
consequence.  It does show that the imported local correction mechanism has
been applied to the actual two-block lift with exact, rather than sampled,
acceptance.
