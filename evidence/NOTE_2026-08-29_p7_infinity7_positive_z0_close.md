# Proposition 15.714: positive p7 infinity-plus-seven z=0 is closed

Let `z` be the number of undetermined directions of the seven-point affine
boundary. In the positive-product branch, `z=0` means every direction has
odd-fibre count `1`, `3`, or `5`. Both labelled types use phase zero and
their exact budget is 32, so all eight means are forced to eight.

At mean eight each odd-fibre mask has a unique complete Johnson-slice slack:

```text
b=1 or 5: A(X)=|X cap B| mod 2,
b=3:      A(X)=(|X cap B|-2)^2.
```

Thus every seven-point boundary fixes all 280 affine right-side cells. The
common 282-by-1,225 edge system has rank 147 modulo seven and 135 audited
left-null dependencies.

`scripts/p7_infinity7_positive_z0_mod7_gpu.py` combinadically unranked every
one of `C(49,7)=85,900,584` boundaries on Soulkiller's V100. It found the
exact actual-boundary histogram

```text
z=0: 79,447,032
z=1:  6,324,528
z=2:    123,480
z=3:      5,488
z=7:         56.
```

All 79,447,032 `z=0` boundaries violate at least one exact mod-seven
dependency; zero survive. A 100,000-boundary CPU prefix separately
reconstructs 88,715 `z=0` cases and also finds zero survivors.
A same-implementation different-grid rerun with 32,768 rather than 65,535
CUDA blocks exactly reproduces every complete count and hash; its separate
200,000-boundary CPU prefix reconstructs 178,533 `z=0` cases and again finds
zero survivors. No independent implementation validation is claimed.

Therefore the positive `z=0` branch is closed. The projected `b`-profile
outer envelope falls from 1,009 to 792 and the actual positive-boundary scope
to 6,453,552. The `z=1,2,3,7` positive branches, the entire negative branch,
residual (ii), and every top-level gate remain open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found the
surrounding Rédei--Szőnyi, Paley-49, and modular-incidence theory but no prior
occurrence of this exact census or mod-seven exclusion. OEIS records the
generic total `C(49,7)` but not the direction histogram. This is attribution
context, not mathematical evidence or a formal priority claim.
