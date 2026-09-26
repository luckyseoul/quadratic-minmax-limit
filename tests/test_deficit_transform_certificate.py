"""Small exact regressions for the new construction, not a signing census.

Run this explicit file on an external compute host. No recorded research
matrix or historical depth result is enumerated by these tests.
"""

from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from deficit_transform_certificate import (  # noqa: E402
    Slab, canonical, construct_row, contains, enumerate_skeleton,
    load_matrix, pivot_margin, simplify, target_for_neutral_extension,
    transform_section, validate_skeleton, write_new_json,
)


def fiber_section_contains(slabs, pivot, point, buffer):
    """Independent interval description, without generating child normals."""
    lows, highs = [], []
    for s in slabs:
        p = s.normal[pivot]
        dot = sum(c * z for i, (c, z) in enumerate(zip(s.normal, point)) if i != pivot)
        b = s.width(buffer)
        if p == 0:
            if abs(dot) >= b:
                return False
        else:
            endpoints = ((-b - dot) / p, (b - dot) / p)
            lows.append(min(endpoints))
            highs.append(max(endpoints))
    return not lows or (max(lows) < 1 and min(highs) > -1
                        and min(highs) - max(lows) > 2)


def energy(matrix, x):
    return sum(matrix[i][j] * x[i] * x[j]
               for i in range(len(x)) for j in range(i))


class TransformCertificateTests(unittest.TestCase):
    def test_primitive_normalization_keeps_all_three_affine_terms(self):
        s = Slab((2, -4), F(2), F(6), F(2), 3)
        c = canonical(s)
        self.assertEqual(c, Slab((1, -2), F(1), F(3), F(1), 3))
        self.assertEqual(c.width(F(3, 2)), s.width(F(3, 2)) / 2)

    def test_retained_slabs_are_needed_for_a_section(self):
        parents = (Slab((1, 2, 2)),)
        child, pairs = transform_section(parents, 0, F(3, 2), prune_cube=False)
        self.assertEqual(pairs, 0)
        self.assertFalse(contains(child, (0, 1, 1), F(3, 2)))
        self.assertTrue(contains(child, (0, 1, -1), F(3, 2)))
        self.assertEqual(child[0].width(F(3, 2)), F(5, 4))

    def test_positive_deficit_and_mixed_phase_survive(self):
        parents = (Slab((1, 1, 1), phases=1),
                   Slab((1, -1, 1), reserve=F(1), phases=2))
        child, _ = transform_section(parents, 0, F(3, 2), prune_cube=False)
        mixed = next(s for s in child if s.normal == (0, 1, 0))
        self.assertEqual((mixed.mass, mixed.reserve, mixed.cost, mixed.phases),
                         (F(1), F(1, 2), F(1), 3))
        self.assertEqual(mixed.width(F(3, 2)), 1)

    def test_exact_sections_match_independent_interval_geometry(self):
        buffer = F(1, 2)
        parents = (Slab((1, 2, -1), reserve=F(2), phases=1),
                   Slab((-2, 1, 1), reserve=F(3), phases=2),
                   Slab((0, 1, -2), reserve=F(1), phases=1))
        checked = 0
        for pivot in range(3):
            if pivot_margin(parents, pivot, buffer) <= 0:
                continue
            child, _ = transform_section(parents, pivot, buffer, prune_cube=False)
            for coordinates in product((F(-2), F(-1), F(-1, 2), F(0),
                                        F(1, 2), F(1), F(2)), repeat=2):
                point = list(coordinates)
                point.insert(pivot, F(0))
                self.assertEqual(contains(child, point, buffer),
                                 fiber_section_contains(parents, pivot, point, buffer))
                checked += 1
        self.assertGreaterEqual(checked, 49)

    def test_strict_boundary_is_not_pruned(self):
        parents = (Slab((0, 1, 0)),)
        self.assertEqual(len(simplify(parents, F(1))), 1)
        self.assertEqual(len(simplify(parents, F(3, 2))), 0)

    def test_stopped_transform_is_not_an_infeasibility_certificate(self):
        parents = (Slab((1, 1, 1)), Slab((1, -1, 1), phases=2))
        buffer = F(3, 2)
        self.assertTrue(contains(parents, (1, 1, -1), buffer))
        result = construct_row(parents, 3, buffer)
        self.assertEqual(result["status"], "TRANSFORM_STOPPED")
        self.assertNotIn("row", result)

    def test_resource_stop_is_explicitly_inconclusive(self):
        parents = (Slab((1, 1, 1)), Slab((1, -1, 1)))
        result = construct_row(parents, 3, F(5, 2), max_facets=1)
        self.assertEqual(result["status"], "RESOURCE_LIMIT")
        result = construct_row(parents, 3, F(5, 2), max_pairs=0)
        self.assertEqual(result["status"], "RESOURCE_LIMIT")

    def test_recovery_satisfies_original_constraints(self):
        parents = (Slab((1, 1, 1), phases=1),
                   Slab((1, -1, 1), phases=2))
        result = construct_row(parents, 3, F(5, 2))
        self.assertEqual(result["status"], "ROW_FOUND")
        self.assertTrue(contains(parents, result["row"], F(5, 2)))
        self.assertEqual(set(map(abs, result["row"])), {1})

    def test_gray_skeleton_matches_direct_definition_on_one_fixture(self):
        matrix = ((0, 1, 1, -1, 1), (1, 0, -1, 1, -1),
                  (1, -1, 0, 1, 1), (-1, 1, 1, 0, -1),
                  (1, -1, 1, -1, 0))
        skeleton = enumerate_skeleton(matrix)
        validate_skeleton(skeleton, matrix)
        expected, norm = {}, 0
        for suffix in product((-1, 1), repeat=4):
            x = (1,) + suffix
            q = energy(matrix, x)
            norm = max(norm, abs(q))
            phases = 0
            for sigma, bit in ((1, 1), (-1, 2)):
                if all(sigma * energy(matrix, x) >= sigma * energy(
                    matrix, x[:i] + (-x[i],) + x[i + 1:]) for i in range(5)):
                    phases |= bit
            if phases:
                expected[x] = (abs(q), phases)
        actual = {tuple(s["x"]): (s["energy"], s["phases"])
                  for s in skeleton["states"]}
        self.assertEqual(actual, expected)
        self.assertEqual(skeleton["phi"], norm)

    def test_recovered_extension_matches_full_cube_on_one_fixture(self):
        matrix = ((0, -1, -1), (-1, 0, -1), (-1, -1, 0))
        skeleton = enumerate_skeleton(matrix)
        phi, target = skeleton["phi"], 4
        buffer = F(2 * (target - phi) + 1, 2)
        slabs = [Slab(tuple(s["x"]), reserve=F(phi - s["energy"]),
                      phases=s["phases"]) for s in skeleton["states"]]
        result = construct_row(slabs, 3, buffer)
        self.assertEqual(result["status"], "ROW_FOUND")
        row = result["row"]
        direct = max(abs(energy(matrix, x) + t * sum(a * b for a, b in zip(row, x)))
                     for x in product((-1, 1), repeat=3) for t in (-1, 1))
        stable = max(s["energy"] + abs(sum(a * b for a, b in zip(row, s["x"])))
                     for s in skeleton["states"])
        self.assertEqual(direct, stable)
        self.assertLessEqual(direct, target)

    def test_exact_neutral_floor_and_parity(self):
        self.assertEqual(target_for_neutral_extension(27, 15, 0), 28)
        self.assertEqual(target_for_neutral_extension(27, 15, 1), 30)
        self.assertEqual(target_for_neutral_extension(3, 3, 0), 4)

    def test_cache_does_not_gain_an_independent_completeness_claim(self):
        matrix = ((0, -1, -1), (-1, 0, -1), (-1, -1, 0))
        skeleton = enumerate_skeleton(matrix)
        skeleton["states"] = [s for s in skeleton["states"] if s["energy"] == 3]
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            write_new_json(folder / "matrix.json", matrix)
            write_new_json(folder / "skeleton.json", skeleton)
            proc = subprocess.run([
                sys.executable, str(ROOT / "scripts/deficit_transform_certificate.py"),
                "--matrix", str(folder / "matrix.json"), "--skeleton",
                str(folder / "skeleton.json"), "--target-phi", "4",
            ], capture_output=True, text=True, check=True)
            result = json.loads(proc.stdout)
            self.assertEqual(result["certificate"],
                             "conditional_on_cached_skeleton_completeness")
            self.assertNotIn("extension_phi", result)

    def test_nonintegral_input_and_output_overwrite_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "matrix.json"
            write_new_json(path, [[0, 1.0], [1.0, 0]])
            with self.assertRaises(ValueError):
                load_matrix(path)
            with self.assertRaises(FileExistsError):
                write_new_json(path, {})


if __name__ == "__main__":
    unittest.main()
