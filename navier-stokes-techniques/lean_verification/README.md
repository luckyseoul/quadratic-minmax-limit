# Independently verifying OpenAI's Navier–Stokes / Euler Lean proofs

This documents the exact steps used to build and independently re-check the
Lean 4 formalization accompanying "Finite time blowup for Navier–Stokes"
(repo: <https://github.com/openai/NavierStokesAndEuler>), on this machine.

Everything below was installed **without root** (no passwordless sudo here).

## Why two kernels

`lake build` already type-checks everything with Lean's own kernel. The extra
value of Comparator is defence in depth:

1. it re-checks the exported proof term with an **independent kernel**
   implementation (`nanoda`, a Rust type checker for Lean 4), so a bug in
   Lean's kernel alone cannot admit a bad proof;
2. it enforces an **axiom allowlist** — for these challenges,
   `propext`, `Quot.sound`, `Classical.choice` only. Notably `sorryAx` is
   *not* permitted, so the proof cannot be vacuous;
3. it **compares the solution's theorem statements against the challenge
   file's statements**, so a solution cannot quietly prove something weaker
   than Clay's alternatives (C) and (D).

Note the challenge file `ComparatorChallenges/NavierStokes.lean` itself
contains `sorry` at the two theorems — that is *by design*: that file only
**states** the problem. The proofs live in `NavierStokes.ComparatorSolution`.

## Toolchain

```sh
# Lean (elan) -- installs to ~/.elan
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y
export PATH="$HOME/.elan/bin:$PATH"
# The repo pins leanprover/lean4:v4.34.0-rc2 via its lean-toolchain file.

# Rust (rustup) -- needed only for nanoda, installs to ~/.cargo
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# Go -- needed only for landrun, user-local install
curl -sL https://go.dev/dl/go1.23.4.linux-amd64.tar.gz -o /tmp/go.tar.gz
tar -C "$HOME/.local" -xzf /tmp/go.tar.gz
export PATH="$HOME/.local/go/bin:$PATH"
```

## Comparator's three external binaries

```sh
# landrun (sandbox), at the commit pinned by lean-eval's instructions
go install github.com/zouuup/landrun/cmd/landrun@5ed4a3db3a4ad930d577215c6b9abaa19df7f99f
export PATH="$HOME/go/bin:$PATH"

# nanoda (independent kernel)
git clone https://github.com/ammkrn/nanoda_lib.git ~/nanoda_lib
cd ~/nanoda_lib && cargo build --release        # -> target/release/nanoda_bin
export PATH="$HOME/nanoda_lib/target/release:$PATH"

# lean4export + comparator are already lake dependencies of the repo
cd ~/NavierStokesAndEuler
lake build lean4export comparator
```

`COMPARATOR_LANDRUN`, `COMPARATOR_LEAN4EXPORT` and `COMPARATOR_NANODA` can
supply explicit paths instead of relying on `PATH`.

## Build and check

```sh
cd ~/NavierStokesAndEuler
lake exe cache get     # mathlib binary cache
lake build             # 11,251 jobs
lake exe comparator ComparatorChallenges/NavierStokes.json
lake exe comparator ComparatorChallenges/Euler.json
```

## The parallelism patch

Stock Comparator runs nanoda **single-threaded**, which is not an inherent
limit — it is a hardcoded omission. Comparator generates nanoda's config in
`Main.lean` (`runExternalKernel`) and does not emit `num_threads`; nanoda's
`Config.num_threads` is `#[serde(default)]` over a `usize`, so it defaults to
`0`, and `check_all_declars` (`nanoda_lib/src/tc.rs:181`) only takes the
parallel branch when `num_threads > 1`.

The patch applied here adds one entry to that JSON object:

```lean
("num_threads", (88 : Nat)),
```

This is safe with respect to what is being verified: nanoda's
`check_all_declars_par` is a shared atomic work queue over **all**
declarations, each running the identical `check_declar` as the serial path.
It changes scheduling only — not which declarations are checked, not how
strictly, and not the axiom allowlist.

Alternatives considered:

* **Run nanoda directly** on a `lean4export` dump with your own config —
  same independent-kernel check, but loses Comparator's sandbox, its
  statement comparison, and its built-in-kernel pass.
* **A `COMPARATOR_NANODA` wrapper script** that rewrites the config — will
  not work: Comparator spawns the kernel under landrun with `envPass := #[]`
  and `executablePaths := #[]`, so a wrapper has no interpreter it is
  permitted to exec.

Note `lean4export` itself is single-threaded and is a serial phase before
the parallel kernel check; that part is inherent to the export step.

## Results on this machine

See `RESULTS.md` in this directory.
