#!/usr/bin/env python3
"""Probe the integral +p eigenspace lattice of the Paley conference matrix.

This is an unnumbered R1 lab tool.  For

    L_p = ker_Z(C_p - p I),   rank(L_p) = (p^2 + 1)/2,

PARI/GP's ``matkerint`` gives a saturated integral basis.  We report the
Gram determinant and Smith invariants, which decide whether the odd-coset
shell in Prop. 15.589 H belongs to a familiar modular-lattice family.
Nothing here asserts R1 or flips a repository predicate.
"""

from __future__ import annotations

import argparse
import ast
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15270 import _make_field


def paley_conference(p: int) -> np.ndarray:
    """Standard symmetric Paley conference matrix on {infinity} union F_{p^2}."""
    field = _make_field(p)
    chi = field["chi"]
    q = p * p
    C = np.ones((q + 1, q + 1), dtype=np.int64)
    C[0, 0] = 0
    for x in range(q):
        for y in range(q):
            if x == y:
                C[x + 1, y + 1] = 0
            else:
                # Coordinates in _make_field are componentwise additive.
                d = ((y % p - x % p) % p) + ((y // p - x // p) % p) * p
                C[x + 1, y + 1] = int(chi(d))
    if not np.array_equal(C, C.T):
        raise ArithmeticError("Paley conference matrix is not symmetric")
    if not np.array_equal(C @ C, q * np.eye(q + 1, dtype=np.int64)):
        raise ArithmeticError("Paley conference identity failed")
    return C


def gp_matrix(A: np.ndarray) -> str:
    rows = [",".join(str(int(x)) for x in row) for row in A]
    return "[" + ";".join(rows) + "]"


def add_field(x: int, y: int, p: int) -> int:
    return ((x % p + y % p) % p) + ((x // p + y // p) % p) * p


def mul_prime_scalar(a: int, x: int, p: int) -> int:
    return ((a * (x % p)) % p) + ((a * (x // p)) % p) * p


def square_line_vectors(p: int) -> tuple[np.ndarray, np.ndarray]:
    """A rational basis of V+ from square affine-line circle words.

    For each square F_p-direction choose p parallel affine lines.  Keep one
    circle and all within-pencil differences.  The result has
    1 + ((p+1)/2)(p-1) = (p^2+1)/2 columns.
    """
    field = _make_field(p)
    mul, chi = field["mul"], field["chi"]
    q = p * p
    unseen = set(range(1, q))
    directions: list[int] = []
    while unseen:
        g = min(unseen)
        line_units = {mul_prime_scalar(a, g, p) for a in range(1, p)}
        unseen.difference_update(line_units)
        if int(chi(g)) == 1:
            directions.append(g)
    if len(directions) != (p + 1) // 2:
        raise ArithmeticError("wrong number of square directions")

    pencils: list[list[np.ndarray]] = []
    for g in directions:
        line = {mul_prime_scalar(a, g, p) for a in range(p)}
        h = next(x for x in range(1, q) if x not in line)
        rows: list[np.ndarray] = []
        for b in range(p):
            shift = mul_prime_scalar(b, h, p)
            points = {add_field(shift, x, p) for x in line}
            v = np.zeros(q + 1, dtype=np.int64)
            v[0] = 1
            v[[x + 1 for x in points]] = 1
            rows.append(v)
        pencils.append(rows)

    cols = [pencils[0][0]]
    for rows in pencils:
        cols.extend(rows[b] - rows[0] for b in range(1, p))
    basis = np.column_stack(cols)
    all_circles = np.column_stack([v for rows in pencils for v in rows])
    return basis, all_circles


def probe(p: int) -> dict:
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    V, Vall = square_line_vectors(p)
    if not np.array_equal(A @ V, np.zeros_like(V)):
        raise ArithmeticError("square-line columns are not +p eigenvectors")
    program = f"""
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
Ginv=G^-1;
DUALDEN=denominator(Ginv);
LEVEL=DUALDEN;
for(i=1,matsize(G)[1],if((LEVEL*Ginv[i,i])%2,LEVEL=2*DUALDEN;break));
V={gp_matrix(V)};
Vall={gp_matrix(Vall)};
GV=V~*V;
X=matsolve(B~*B,B~*Vall);
H=mathnf(X);
print("RANK=",matsize(B)[2]);
print("DET=",matdet(G));
print("SNF=",matsnf(G));
print("CONTENT=",content(G));
print("DUALDEN=",DUALDEN);
print("LEVEL=",LEVEL);
print("LINEDET=",matdet(GV));
print("LINEINDEX=",sqrtint(matdet(GV)/matdet(G)));
print("ALLCIRCLEINDEX=",abs(matdet(H)));
quit;
"""
    proc = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    parsed: dict[str, object] = {"p": p, "n": p * p + 1}
    for line in proc.stdout.splitlines():
        key, value = line.split("=", 1)
        key = key.strip().lower()
        value = value.strip()
        if key == "snf":
            parsed[key] = ast.literal_eval(value.replace("[", "[").replace("]", "]"))
        else:
            parsed[key] = int(value)
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", nargs="*", type=int, default=[3, 5, 7, 11])
    args = parser.parse_args()
    for p in args.p:
        row = probe(p)
        print(
            f"p={p} n={row['n']} rank={row['rank']} det={row['det']} "
            f"content={row['content']} snf={row['snf']} "
            f"dual_den={row['dualden']} level={row['level']} "
            f"line_det={row['linedet']} line_index={row['lineindex']} "
            f"all_circle_index={row['allcircleindex']}"
        )


if __name__ == "__main__":
    main()
