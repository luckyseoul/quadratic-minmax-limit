"""Independent exact mesh-checker tests; no candidate theorem imports."""
import hashlib
import json
from functools import lru_cache
from pathlib import Path

from scripts.p1_third_post_band_mesh_check import MODES


ROOT = Path(__file__).resolve().parents[1]
REPLAY = ROOT / 'evidence/p1_third_post_band_mesh_replay.json'
REPLAY_SHA256 = '41add1e11eaf0bb041060b6eca973fd44e09194c0ce1e46b46d289a0f20e643a'
EXPECTED_RUNS = {
    'quadrature': ('soulkiller', 'x86_64', '192.168.1.113', '/home/nick/.venvs/mo-exact/bin/python'),
    'contact-kernel': ('NUKA', 'x86_64', '192.168.1.192', '/home/nick/.venvs/qml-mesh/bin/python'),
    'cube-equality': ('jellyfin', 'x86_64', '192.168.1.191', '/home/nick/.venvs/mo-intel/bin/python'),
    'branch-ledger': ('orin', 'aarch64', '192.168.1.135', '/usr/bin/python3'),
}


@lru_cache(maxsize=None)
def live_result(mode):
    # The evidence binding below reuses the existing checks, including rank403.
    return MODES[mode]()


def test_prime_free_quadrature_neighboring_slice_and_shifted_gap():
    result = live_result('quadrature')
    assert result['shifted_gap_numerator_coefficients'] == [128, 30, 1]
    assert len(result['rows']) == 7
    assert result['rows'][0]['single_nonzero_contact_contribution_lower_bound'] == '132/25'
    assert 33 in [row['odd_order'] for row in result['rows']]


def test_p29_general_slice_odd_contact_kernel():
    result = live_result('contact-kernel')
    assert result['raw_coefficients'] == 436
    assert result['contact_rank'] == 403
    assert result['explicit_independent_kernel_vectors'] == 33
    assert len(result['independent_contact_masks']) == 403


def test_punctured_difference_not_assumed_nonnegative():
    result = live_result('cube-equality')
    assert len(result['forms']) == 4
    new_forms = [row for row in result['forms'] if row['delta'] == 4]
    assert len(new_forms) == 3
    assert all(row['difference_truth_table'][0] == -1 for row in new_forms)
    assert all(row['signed_offset'] == 4 for row in new_forms)


def test_all_eight_branch_ledgers_at_four_representative_primes():
    result = live_result('branch-ledger')
    assert [row['p'] for row in result['rows']] == [29, 37, 41, 53]
    for row in result['rows']:
        assert len(row['branches']) == 8
        assert row['original_k'] == 5*row['p']-3
        assert [branch['forced_rows_at_least'] for branch in row['branches']] == [5, 6, 9, 8, 8, 7, 7, 7]


def test_recorded_mesh_runs_bind_exact_checker_and_distinct_capabilities():
    replay_bytes = REPLAY.read_bytes()
    assert hashlib.sha256(replay_bytes).hexdigest() == REPLAY_SHA256
    replay = json.loads(replay_bytes)
    script_hash = hashlib.sha256((ROOT / replay['script']).read_bytes()).hexdigest()
    assert script_hash == replay['script_sha256']
    assert len(replay['runs']) == 4
    outputs = [run['output'] for run in replay['runs']]
    # output holds each actual node's already-parsed stdout JSON, not a log string.
    assert all(isinstance(output, dict) for output in outputs)
    assert {output['mode'] for output in outputs} == set(EXPECTED_RUNS)
    assert len({output['host'] for output in outputs}) == 4
    for run in replay['runs']:
        output = run['output']
        host, architecture, ip, interpreter = EXPECTED_RUNS[output['mode']]
        assert run['exit_code'] == 0
        assert output['passed'] is True
        assert output['host'] == host
        assert output['architecture'] == architecture
        assert output['script_sha256'] == script_hash
        assert isinstance(output['result'], dict)
        assert f'nick@{ip}' in run['command']
        assert interpreter in run['command']
        assert run['command'].endswith(f"p1_third_post_band_mesh_check.py {output['mode']}'")
        assert host in replay['routing']
    assert 'ARM64' in replay['routing']['orin']
    assert replay['gpu_research_jobs_launched'] is False
    assert replay['old_p23_certificate_rerun'] is False


def test_all_recorded_mode_results_match_live_exact_checks():
    replay = json.loads(REPLAY.read_bytes())
    for run in replay['runs']:
        output = run['output']
        # Normalize tuples as JSON arrays before comparing the complete payloads.
        assert output['result'] == json.loads(json.dumps(live_result(output['mode'])))
