"""Independent exact bridge checker regression tests."""
import hashlib
import json
from pathlib import Path

from scripts.p23_third_post_band_mesh_check import (
    check_contacts, check_cubes, check_ledger, check_small_kernel,
)


def test_rational_contacts_and_independent_floor_lp():
    result = check_contacts()
    assert result['positive_degree_two_quadratures'] == 9
    assert all(row['independent_LP_floor'] == '46' for row in result['rows'])


def test_complete_small_support_kernel_certificate():
    result = check_small_kernel()
    for row in result['boundaries'].values():
        assert row['raw_coefficients'] == 277
        assert row['contact_rank'] == 250
        assert row['explicit_independent_kernel_vectors'] == 27
        assert len(row['independent_contact_masks']) == 250


def test_middle_boundary_covering_cubes_and_walsh_gram():
    result = check_cubes()
    assert [row['dimension'] for row in result['rows']] == [6, 8, 10, 11, 9, 7, 5]
    assert all(row['Walsh_Gram_off_diagonal'] == 0 for row in result['rows'])


def test_phase_zero_floor_and_offset_arithmetic():
    result = check_ledger()
    assert result['mass32_floor_surviving_boundaries'] == [0, 2, 22]
    assert result['nonzero_boundary_lift_excess'] == 8
    assert all(row['forced_rows_at_least'] == 5 for row in result['offset_ledgers'])


def test_saved_four_node_provenance_matches_the_exact_checker():
    root = Path(__file__).resolve().parents[1]
    replay = json.loads((root / 'evidence/p23_third_post_band_mesh_replay.json').read_text())
    digest = hashlib.sha256((root / replay['script']).read_bytes()).hexdigest()
    assert digest == replay['script_sha256']
    assert len(replay['runs']) == 4
    expected = {'contacts': 'soulkiller', 'small-kernel': 'nuka',
                'cubes': 'jellyfin', 'ledger': 'orin'}
    outputs = {run['output']['mode']: run for run in replay['runs']}
    assert set(outputs) == set(expected)
    for mode, host in expected.items():
        run = outputs[mode]
        assert run['exit_code'] == 0 and run['output']['passed']
        assert run['output']['host'].lower() == host
        assert run['output']['script_sha256'] == digest
    assert outputs['ledger']['output']['architecture'] == 'aarch64'
