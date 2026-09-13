# Symbolic cutoff receipt

`check_bounds.py` was staged with SHA-256
`2b4da2f3832cff07e1d549627141926d3c9b39b6208d4b73f88b42d098d546f3`
and run once on NUKA in the isolated virtual environment
`/tmp/qml-entropy-proof.KZwfQE4I/venv`.

The command used one CPU thread and completed with exit status 0 in 0.237 s:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\
  /tmp/qml-entropy-proof.KZwfQE4I/venv/bin/python check_bounds.py
```

`z3-solver` 5.1.0 reported `unsat` for both quantified real-algebra
counterexample searches. The raw structured output is `result.json`.
This is scalar-arithmetic corroboration only: it checks the two cutoff
inequalities used in Section 2 of the parent note, not Bernstein's inequality,
the entropy argument, or the all-orders convergence criterion.
